module depletion

   use sdata, only: dp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: depletion_run

   ! Critical boron search safeguards (boron in ppm, slope in dk/dppm).
   real(dp), parameter :: bcon_max = 3000._dp           ! Upper boron limit (ppm), as in check_ppm
   real(dp), parameter :: max_boron_step = 500._dp      ! Max. boron change per search iteration (ppm)
   real(dp), parameter :: slope_default = -8.e-5_dp     ! Fallback boron worth (dk/ppm)
   real(dp), parameter :: slope_lo = -3.e-4_dp          ! Most negative accepted secant slope (dk/ppm)
   real(dp), parameter :: slope_hi = -2.e-5_dp          ! Least negative accepted secant slope (dk/ppm)
   integer, parameter :: bound_hits_limit = 2           ! Consecutive non-critical solves at a bound before stopping

contains

   subroutine depletion_solve(negative_boron)
      ! Solve one depletion interval at the current burnup and operating state.
      !
      ! BCSEARCH uses a secant search with safeguards against feedback noise:
      ! nonphysical secant slopes are replaced by a typical boron worth, each step
      ! is limited and clamped to [0, bcon_max], and a bound is only reported when
      ! consecutive solves at that bound stay non-critical.

      use sdata, only: mode, bcon, rbcon, ke, ftem, mtem, cden, bpos, &
         th_niter, ser, fer, serc, ferc
      use io, only: bther
      use xsec, only: xs_updt
      use cmfd, only: outer
      use th, only: th_iter
      use control, only: check_ppm

      implicit none

      real(dp) :: bc1, bc2, ke1, ke2, slope
      integer :: iteration, zero_hits, max_hits
      logical, intent(out) :: negative_boron

      negative_boron = .false.
      zero_hits = 0
      max_hits = 0

      if (mode == 'FORWARD') then
         if (bther == 1) then
            call th_iter(th_niter, 0)
         else
            call xs_updt(bcon, ftem, mtem, cden, bpos)
            call outer(0)
         end if
         return
      end if

      bcon = rbcon
      call solve_at_boron()
      bc1 = bcon
      ke1 = ke
      if (critical_state()) return
      if (stop_at_bound()) return

      if (bther == 1 .and. bcon < 1.e-5_dp) then
         bcon = 500._dp
      else
         bcon = bcon + (ke - 1._dp) * bcon
      end if
      call limit_guess(bc1)
      call check_ppm(2, bcon)
      call solve_at_boron()
      bc2 = bcon
      ke2 = ke
      if (critical_state()) return
      if (stop_at_bound()) return

      iteration = 2
      do
         ! Use the secant slope only when it is a physical boron worth; otherwise
         ! feedback noise dominates the keff difference and would cause an overshoot.
         slope = slope_default
         if (bc1 /= bc2) then
            if ((ke1 - ke2) / (bc1 - bc2) >= slope_lo .and. &
               (ke1 - ke2) / (bc1 - bc2) <= slope_hi) slope = (ke1 - ke2) / (bc1 - bc2)
         end if

         bcon = bc2 + (1._dp - ke2) / slope
         call limit_guess(bc2)
         iteration = iteration + 1
         call check_ppm(iteration, bcon)
         call solve_at_boron()

         bc1 = bc2
         bc2 = bcon
         ke1 = ke2
         ke2 = ke

         if (critical_state()) exit
         if (stop_at_bound()) return
      end do

   contains

      subroutine solve_at_boron()
         if (bther == 1) then
            call th_iter(2, 0)
         else
            call xs_updt(bcon, ftem, mtem, cden, bpos)
            call outer(0)
         end if
      end subroutine solve_at_boron

      logical function critical_state()
         if (bther == 1) then
            critical_state = abs(ke - 1._dp) < 1.e-5_dp .and. ser < serc .and. fer < ferc
         else
            critical_state = abs(ke - 1._dp) < 1.e-5_dp .and. ser < 1.e-5_dp .and. fer < 1.e-5_dp
         end if
      end function critical_state

      subroutine limit_guess(b_ref)
         ! Limit the change from b_ref and keep the guess within [0, bcon_max].
         real(dp), intent(in) :: b_ref

         bcon = b_ref + max(-max_boron_step, min(max_boron_step, bcon - b_ref))
         bcon = max(0._dp, min(bcon_max, bcon))
      end subroutine limit_guess

      logical function stop_at_bound()
         ! A single non-critical solve at a bound can be feedback noise, so stop
         ! only when consecutive solves at the bound stay on the same side.
         stop_at_bound = .false.

         if (bcon <= 0._dp .and. ke < 1._dp) then
            zero_hits = zero_hits + 1
            if (zero_hits >= bound_hits_limit) then
               negative_boron = .true.
               stop_at_bound = .true.
            end if
         else
            zero_hits = 0
         end if

         if (bcon >= bcon_max .and. ke > 1._dp) then
            max_hits = max_hits + 1
            if (max_hits >= bound_hits_limit) &
               call depletion_error('critical boron concentration exceeds the limit (3000 ppm)')
         else
            max_hits = 0
         end if
      end function stop_at_bound

   end subroutine depletion_solve

   !******************************************************************************!

   subroutine depletion_run()
      use sdata, only: mode, nnod, nmat, mat, vdel, bpos, pow, ppow, npow, ke, bcon, &
         ftem, mtem, cden, aprad, apaxi, afrad, n_depl_intervals, depl_days, &
         depl_power_percent, depl_bank_step
      use io, only: ounit, scr, bther, bcrod, bcycl, asmpow, axipow, asmflux, &
         boutp, print_outp, bvtk, print_vtk
      use dlib, only: dlib_get_node_burnup, dlib_update_node_history, &
         dlib_get_node_history, dlib_get_material_name, dlib_get_hm_density, dlib_is_fuel
      use xsec, only: xs_updt
      use control, only: get_power_dist, print_tail
      use hdf5_output, only: hdf5_is_active, hdf5_write_step
      use cycle, only: cycle_write_state

      implicit none

      integer :: interval, node, nbank
      logical :: negative_boron
      real(dp) :: elapsed_days, node_hm_mass, burnup_increment, interval_power_w, interval_power_mw
      real(dp) :: total_hm_mass, burnup_inventory, average_burnup, maximum_burnup
      real(dp), dimension(:), allocatable :: burnup
      real(dp), dimension(:), allocatable :: cycle_hm_density, cycle_hboron, cycle_htmod, cycle_htfuel
      character(len = 64), dimension(:), allocatable :: cycle_material_name
      logical, dimension(:), allocatable :: cycle_fuel_material

      if (n_depl_intervals < 1) call depletion_error('depletion schedule is empty')
      if (size(depl_days) /= n_depl_intervals .or. &
         size(depl_power_percent) /= n_depl_intervals) then
         call depletion_error('depletion schedule dimensions are inconsistent')
      end if
      nbank = size(depl_bank_step, 1)
      if (.not. allocated(npow)) allocate(npow(nnod))
      allocate(burnup(nnod))
      call dlib_get_node_burnup(burnup)

      write(ounit, *)
      write(ounit, '(A)') '=============================================================='
      write(ounit, '(A)') '                 DEPLETION CALCULATION RESULTS'
      write(ounit, '(A)') '=============================================================='
      call write_depletion_table_header(ounit, nbank, mode == 'BCSEARCH')
      if (scr) then
         write(*, *)
         write(*, '(A)') '=============================================================='
         write(*, '(A)') '                 DEPLETION CALCULATION RESULTS'
         write(*, '(A)') '=============================================================='
         call write_depletion_table_header(6, nbank, mode == 'BCSEARCH')
      end if

      elapsed_days = 0._dp
      do interval = 1, n_depl_intervals
         if (bcrod == 1) bpos = depl_bank_step(:, interval)

         ! The schedule overrides %THER's power level; pow remains rated full power.
         ppow = depl_power_percent(interval)
         interval_power_w = pow * ppow * 0.01_dp
         interval_power_mw = interval_power_w * 1.e-6_dp

         call depletion_solve(negative_boron)
         if (negative_boron) then
            ! Keep output routines and the thermal summary on the last valid solve.
            call get_power_dist(npow)
            write(ounit, '(A,I0,A)') &
               '  WARNING: DEPLETION TERMINATED IN INTERVAL ', interval, &
               ' BECAUSE CRITICAL BORON CONCENTRATION WAS LESS THAN ZERO'
            if (scr) write(*, '(A,I0,A)') &
               '  WARNING: DEPLETION TERMINATED IN INTERVAL ', interval, &
               ' BECAUSE CRITICAL BORON CONCENTRATION WAS LESS THAN ZERO'
            exit
         end if
         call get_power_dist(npow)

         ! The solved state corresponds to the burnup at the start of this interval.
         if (hdf5_is_active) then
            if (mode == 'BCSEARCH') then
               call hdf5_write_step(interval - 1, elapsed_days * 86400._dp, &
                  boron_concentration = bcon, power_w = interval_power_w, node_burnup = burnup)
            else
               call hdf5_write_step(interval - 1, elapsed_days * 86400._dp, &
                  power_w = interval_power_w, node_burnup = burnup)
            end if
         end if

         total_hm_mass = 0._dp
         burnup_inventory = 0._dp
         maximum_burnup = 0._dp
         do node = 1, nnod
            if (.not. dlib_is_fuel(mat(node))) cycle
            node_hm_mass = dlib_get_hm_density(mat(node)) * vdel(node) * 1.e-6_dp
            if (node_hm_mass <= 0._dp) call depletion_error('fuel node has no positive HM mass')

            ! MWth times days gives MWd; dividing by node HM mass gives MWd/kgHM.
            burnup_increment = npow(node) * interval_power_mw * &
               depl_days(interval) / node_hm_mass
            if (.not. ieee_is_finite(burnup_increment)) then
               call depletion_error('burnup increment is not finite')
            end if
            burnup(node) = burnup(node) + burnup_increment
            total_hm_mass = total_hm_mass + node_hm_mass
            burnup_inventory = burnup_inventory + node_hm_mass * burnup(node)
            maximum_burnup = max(maximum_burnup, burnup(node))
         end do
         if (total_hm_mass <= 0._dp) call depletion_error('the core contains no fuel HM mass')
         average_burnup = burnup_inventory / total_hm_mass
         call dlib_update_node_history(burnup, bcon, ftem, mtem)

         elapsed_days = elapsed_days + depl_days(interval)
         call write_depletion_table_row(ounit, interval, elapsed_days, interval_power_mw, ke, &
            bcon, average_burnup, maximum_burnup, depl_bank_step(:, interval), mode == 'BCSEARCH')
         if (scr) call write_depletion_table_row(6, interval, elapsed_days, interval_power_mw, ke, &
            bcon, average_burnup, maximum_burnup, depl_bank_step(:, interval), mode == 'BCSEARCH')

         if (bvtk == 1) call print_vtk(interval - 1)
         ! Check the final accumulated state against DLIB bounds without another transport solve.
         if (interval == n_depl_intervals) call xs_updt(bcon, ftem, mtem, cden, bpos)
      end do

      write(ounit, '(A)') ' Note: Burnup is in MWd/kgHM.'
      write(ounit, '(A)') ' Keff, power shape, and bank positions are at interval start; burnup is after each interval.'
      if (scr) write(*, '(A)') ' Note: Burnup is in MWd/kgHM.'
      if (scr) write(*, '(A)') ' Keff, power shape, and bank positions are at interval start; burnup is after each interval.'

      if (aprad == 1) call asmpow(npow)
      if (apaxi == 1) call axipow(npow)
      if (afrad == 1) call asmflux(1._dp)
      if (boutp == 1) call print_outp(npow)
      if (bther == 1) then
         write(ounit, '(A)') ''
         write(ounit, '(A)') ' THERMAL FEEDBACK SUMMARY (LAST INTERVAL SOLVE)'
         if (scr) write(*, '(A)') ''
         if (scr) write(*, '(A)') ' THERMAL FEEDBACK SUMMARY (LAST INTERVAL SOLVE)'
         call print_tail()
      end if

      if (bcycl == 1) then
         allocate(cycle_hm_density(nmat), cycle_material_name(nmat), cycle_fuel_material(nmat))
         allocate(cycle_hboron(nnod), cycle_htmod(nnod), cycle_htfuel(nnod))
         call dlib_get_node_history(cycle_hboron, cycle_htmod, cycle_htfuel)
         do node = 1, nmat
            cycle_fuel_material(node) = dlib_is_fuel(node)
            cycle_hm_density(node) = dlib_get_hm_density(node)
            call dlib_get_material_name(node, cycle_material_name(node))
         end do
         call cycle_write_state(ounit, scr, cycle_fuel_material, cycle_hm_density, &
            cycle_material_name, burnup, cycle_hboron, cycle_htmod, cycle_htfuel)
         deallocate(cycle_hm_density, cycle_hboron, cycle_htmod, cycle_htfuel)
         deallocate(cycle_material_name, cycle_fuel_material)
      end if

      deallocate(burnup)

   end subroutine depletion_run

   !******************************************************************************!

   subroutine write_depletion_table_header(unit, nbank, include_boron)
      implicit none

      integer, intent(in) :: unit, nbank
      logical, intent(in) :: include_boron

      if (include_boron .and. nbank > 0) then
         write(unit, '(A4,3X,A10,3X,A11,3X,A10,3X,A10,3X,A10,3X,A10,3X,A20)') &
            'Step', 'Elapsed(d)', 'Power(MWth)', 'Keff@start', 'Boron(ppm)', &
            'AvgBU@end', 'MaxBU@end', 'Bank position(steps)'
      else if (include_boron) then
         write(unit, '(A4,3X,A10,3X,A11,3X,A10,3X,A10,3X,A10,3X,A10)') &
            'Step', 'Elapsed(d)', 'Power(MWth)', 'Keff@start', 'Boron(ppm)', 'AvgBU@end', 'MaxBU@end'
      else if (nbank > 0) then
         write(unit, '(A4,3X,A10,3X,A11,3X,A10,3X,A10,3X,A10,3X,A20)') &
            'Step', 'Elapsed(d)', 'Power(MWth)', 'Keff@start', 'AvgBU@end', 'MaxBU@end', 'Bank position(steps)'
      else
         write(unit, '(A4,3X,A10,3X,A11,3X,A10,3X,A10,3X,A10)') &
            'Step', 'Elapsed(d)', 'Power(MWth)', 'Keff@start', 'AvgBU@end', 'MaxBU@end'
      end if

   end subroutine write_depletion_table_header

   !******************************************************************************!

   subroutine write_depletion_table_row(unit, interval, elapsed_days, power_mw, keff, &
      boron, average_burnup, maximum_burnup, bank_positions, include_boron)
      implicit none

      integer, intent(in) :: unit, interval
      real(dp), intent(in) :: elapsed_days, power_mw, keff, boron, average_burnup, maximum_burnup
      real(dp), dimension(:), intent(in) :: bank_positions
      logical, intent(in) :: include_boron
      integer :: bank

      if (include_boron) then
         if (size(bank_positions) > 0) then
            write(unit, '(I4,3X,F10.3,3X,F11.3,3X,F10.5,3X,F10.3,3X,F10.4,3X,F10.4,3X)', &
               advance = 'no') interval, elapsed_days, power_mw, keff, boron, average_burnup, maximum_burnup
         else
            write(unit, '(I4,3X,F10.3,3X,F11.3,3X,F10.5,3X,F10.3,3X,F10.4,3X,F10.4)') &
               interval, elapsed_days, power_mw, keff, boron, average_burnup, maximum_burnup
            return
         end if
      else
         if (size(bank_positions) > 0) then
            write(unit, '(I4,3X,F10.3,3X,F11.3,3X,F10.5,3X,F10.4,3X,F10.4,3X)', &
               advance = 'no') interval, elapsed_days, power_mw, keff, average_burnup, maximum_burnup
         else
            write(unit, '(I4,3X,F10.3,3X,F11.3,3X,F10.5,3X,F10.4,3X,F10.4)') &
               interval, elapsed_days, power_mw, keff, average_burnup, maximum_burnup
            return
         end if
      end if
      do bank = 1, size(bank_positions)
         if (bank == 1) then
            write(unit, '(F10.2)', advance = 'no') bank_positions(bank)
         else
            write(unit, '(1X,F10.2)', advance = 'no') bank_positions(bank)
         end if
      end do
      write(unit, '(A)') ''

   end subroutine write_depletion_table_row

   !******************************************************************************!

   subroutine depletion_error(message)
      use io, only: ounit
      character(len = *), intent(in) :: message

      write(ounit, *) '  ERROR: ', trim(message)
      write(ounit, *) '  KOMODO IS STOPPING'
      write(*, *) '  ERROR: ', trim(message)
      write(*, *) '  KOMODO IS STOPPING'
      stop
   end subroutine depletion_error

   !******************************************************************************!

end module depletion
