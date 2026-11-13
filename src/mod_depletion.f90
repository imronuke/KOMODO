module depletion

   use sdata, only: dp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: depletion_run

contains

   subroutine depletion_solve(negative_boron)
      ! Solve one depletion interval at the current burnup and operating state.

      use sdata, only: mode, bcon, rbcon, ke, ftem, mtem, cden, bpos, &
         th_niter, ser, fer, serc, ferc
      use io, only: bther
      use xsec, only: xs_updt
      use cmfd, only: outer
      use th, only: th_iter
      use control, only: check_ppm

      implicit none

      real(dp) :: bc1, bc2, ke1, ke2, denominator
      integer :: iteration
      logical :: converged
      logical, intent(out) :: negative_boron

      negative_boron = .false.

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

      if (bther == 1 .and. bcon < 1.e-5_dp) then
         bcon = 500._dp
      else
         bcon = bcon + (ke - 1._dp) * bcon
      end if
      if (bcon < 0._dp) then
         bcon = bc1
         negative_boron = .true.
         return
      end if
      call check_ppm(2, bcon)
      call solve_at_boron()
      bc2 = bcon
      ke2 = ke
      if (critical_state()) return

      iteration = 2
      do
         denominator = ke1 - ke2
         if (abs(denominator) <= epsilon(denominator) * max(1._dp, abs(ke1), abs(ke2))) then
            call depletion_error('BCSEARCH secant iteration has indistinguishable keff values')
         end if

         bcon = bc2 + (1._dp - ke2) / denominator * (bc1 - bc2)
         iteration = iteration + 1
         if (bcon < 0._dp) then
            bcon = bc2
            negative_boron = .true.
            return
         end if
         call check_ppm(iteration, bcon)
         call solve_at_boron()

         bc1 = bc2
         bc2 = bcon
         ke1 = ke2
         ke2 = ke

         converged = critical_state()
         if (converged) exit
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
      logical :: negative_boron, schedule_completed
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
      schedule_completed = .true.
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
            schedule_completed = .false.
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

      if (schedule_completed .and. bcycl == 1) then
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
