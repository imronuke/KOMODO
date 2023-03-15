module control
    use iso_fortran_env, only: real64, output_unit, error_unit
    use data
    use print
    use fdm
    use xsec
    use time
    use th
    use transient
    use auxiliaries, only: set_aux_pointer, th_iteration, print_fail_converge
    use utilities

    implicit none

    private

    save

    integer, parameter                    :: dp = real64
    type(fdm_rect), public                :: fdm            ! Finite Difference Object
    type(th_type), allocatable, public    :: th             ! Thermal-hydraulics Object
    type(trans_type), allocatable, public :: tr             ! Transient Object
    type(xs_change_type), public          :: xsc            ! object for xs changes due to a parameter change

    public :: forward, adjoint, fixedsrc, cbc_search, rod_eject

    contains

    !===============================================================================================!
    ! Transient / Rod Ejection calculation
    !===============================================================================================!

    subroutine rod_eject()

        real(dp) :: fn(nnod)

        call calculation_init()

        if (bther == NO) then
            call rod_eject_transient(tr, is_th=.false.)
        else
            call rod_eject_transient(tr, is_th=.true.)
        end if

        call power_dist_print(fdm, fn)
        call print_power_map(fn)

    end subroutine

    !===============================================================================================!
    ! forward calculation
    !===============================================================================================!

    subroutine forward()

        integer  :: g, n
        real(dp) :: fn(nnod)
        logical  :: converge
        real(dp) :: ftemp_err

        call calculation_init()
        call print_head()
        
        if (bther == NO) then
            call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bank_pos)
            call outer_iter(fdm, print_iter = .true., is_converge=converge)
        else
            call th_iteration(25, .true., ftemp_err, converge)
            call print_tail(fdm)
        end if
        if (.not. converge) call print_fail_converge()

        call power_dist_print(fdm, fn)
        call print_power_map(fn)

    end subroutine

    !===============================================================================================!
    ! forward calculation
    !===============================================================================================!

    subroutine fixedsrc()

        integer  :: g, n
        real(dp) :: fn(nnod)
        logical  :: converge

        call calculation_init()
        call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bank_pos)

        call print_head()

        call outer_fixed_src(fdm, print_iter = .true., is_converge=converge)
        if (.not. converge) call print_fail_converge()

        call power_dist_print(fdm, fn)
        call print_power_map(fn)

    end subroutine

    !===============================================================================================!
    ! adjoint calculation
    !===============================================================================================!

    subroutine adjoint()

        integer  :: g, n
        real(dp) :: fn(nnod)
        logical  :: converge

        call calculation_init()
        call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bank_pos)

        call print_head()

        call outer_adjoint(fdm, print_iter = .true., is_converge=converge)
        if (.not. converge) call print_fail_converge()

        call power_dist_print(fdm, fn)
        call print_power_map(fn)

    end subroutine

    !===============================================================================================!
    ! critical boron search calculation
    !===============================================================================================!

    subroutine cbc_search()

        real(dp) :: fn(nnod)

        if (bther == YES) then
            call cbc_search_th()
        else
            call cbc_search_no_th()
        end if

        call power_dist_print(fdm, fn)
        call print_power_map(fn)

    end subroutine

    !===============================================================================================!
    ! critical boron search calculation without TH feedback
    !===============================================================================================!

    subroutine cbc_search_no_th()

        integer  :: n, n_iter
        real(dp) :: K1, K2
        real(dp) :: bc1, bc2
        logical  :: converge

        call calculation_init()

        call print_head()

        n_iter = 2 * nodal_interval

        ! First try
        bcon = xsc % bcon % ref
        call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bank_pos)
        call outer_iter(fdm, print_iter = .false., max_iter=n_iter, is_converge=converge)
        bc1 = bcon
        K1  = fdm % Keff
        write(output_unit,101) 1, bcon, K1, fdm % flux_diff, fdm % fsrc_diff

        ! Second try
        bcon = bcon + (K1 - 1.) * bcon   ! Guess next critical boron concentration
        call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bank_pos)
        call outer_iter(fdm, print_iter = .false.,  max_iter=n_iter, is_converge=converge)
        bc2 = bcon
        K2  = fdm % Keff
        write(output_unit,101) 2, bcon, K2, fdm % flux_diff, fdm % fsrc_diff

        n = 2
        do
            n = n + 1
            bcon = bc2 + (1. - K2) / (K1 - K2) * (bc1 - bc2)
            call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bank_pos)
            call outer_iter(fdm, print_iter = .false.,  max_iter=n_iter, is_converge=converge)
            bc1 = bc2
            bc2 = bcon
            K1  = K2
            K2  = fdm % Keff
            write(output_unit,101) n, bcon, K2, fdm % flux_diff, fdm % fsrc_diff
            if (converge .and. abs(K2-K1) < 1.e-5) exit
        end do

        101 format(I3, F10.2, F14.5, ES14.5, ES13.5)

    end subroutine

    !===============================================================================================!
    ! critical boron search calculation with TH feedback
    !===============================================================================================!

    subroutine cbc_search_th()

        integer   :: n
        real(dp)  :: K1, K2
        real(dp)  :: bc1, bc2
        real(dp)  :: ftem_diff
        logical   :: converge

        call calculation_init()

        call print_head()

        ! First try
        if (allocated(xsc % bcon)) then
            bcon = xsc % bcon % ref
        else
            bcon = 0.0              ! if tabular xs used, set initial boron conc. to 0.0 ppm
        end if
        call th_iteration(n_th_iter, .false., ftem_diff, converge)
        bc1 = bcon
        K1  = fdm % Keff
        write(output_unit,102) 1, bcon, K1, fdm % flux_diff, fdm % fsrc_diff, ftem_diff

        ! Second try
        ! Guess next critical boron concentration
        if (allocated(xsc % bcon)) then
            bcon = bcon + (K1 - 1.) * bcon
        else
            bcon = 500.0              ! if tabular xs used, set initial boron conc. to 0.0 ppm
        end if
        call th_iteration(n_th_iter, .false., ftem_diff, converge)
        bc2 = bcon
        K2  = fdm % Keff
        write(output_unit,102) 2, bcon, K2, fdm % flux_diff, fdm % fsrc_diff, ftem_diff

        n = 2
        do
            n = n + 1
            bcon = bc2 + (1. - K2) / (K1 - K2) * (bc1 - bc2)
            call th_iteration(n_th_iter, .false., ftem_diff, converge)
            bc1 = bc2
            bc2 = bcon
            K1  = K2
            K2  = fdm % Keff
            write(output_unit,102) n, bcon, K2, fdm % flux_diff, fdm % fsrc_diff, ftem_diff
            if (converge .and. abs(K2-K1) < 1.e-5) exit
        end do

        call print_tail(fdm)

        102 format(I3, F9.2, F14.5, ES14.5, ES13.5, ES17.5)

    end subroutine

    !===============================================================================================!
    ! Initialize thermal-hydraulics solver
    !===============================================================================================!

    subroutine th_init()

        ! set data in the th object
        call set_th_data(th, ounit, nzz, zdel, rf, tg, tc, ppitch, cflow, cf, t_inlet)

        ! Allocate necessary data
        allocate(enthalpy(nnod))
        allocate(tfm(nnod, th % n_ring+1))
        tfm = 565.0  ! Initial guess (Kelvin)
        allocate(heat_flux(nnod))
        heat_flux = 0.0
        
        if (.not. allocated(ftem)) then
            allocate(ftem(nnod))
            ftem = 565.0  ! Initial guess (Kelvin)
        endif
        if (.not. allocated(mtem)) then
            allocate(mtem(nnod))
            mtem = 565.0  ! Initial guess (Kelvin)
        endif
        if (.not. allocated(cden)) then
            allocate(cden(nnod))
            cden = 0.734  ! Initial guess (g/cm3)
        endif

        if (allocated(tr)) then
            allocate(flow_rate(nnod))
            flow_rate = cflow         ! initialize sub channel flow rate
        end if

        call set_aux_pointer(fdm, th, xsc)

    end subroutine

    !===============================================================================================!
    ! Initialize calculation
    !===============================================================================================!

    subroutine calculation_init()

        integer  :: i_mat
        real(dp), dimension(nmat,nf) :: lambda_tmp, beta_tmp
        real(dp), dimension(nmat,ng) :: neut_velo_tmp

        allocate(flux(nnod, ng))
        allocate(fsrc(nnod))
        if (.not. allocated(exsrc)) then
            allocate(exsrc(nnod, ng))
            exsrc = 0.0
        end if

        if (.not. allocated(dc)) then
            allocate(dc(nnod, ng, 6))
            dc = 1.0
        end if

        call fdm_time % on
        call set_fdm_data(fdm, ounit, kern, ng, nnod, nxx, nyy, nzz, nmat, &
        east, west, north, south, top, bottom, nodal_interval, 1.e-5_dp, 1.e-5_dp, &
        max_outer, max_outer_th, max_inner, extrp_interval)
        call set_fdm_pointer(fdm, flux, fsrc, exsrc, ind_mat, xdel, ydel, zdel, vdel, &
        xstag, ystag, ix, iy, iz, xyz)
        call fdm_time % off

        call xs_time % on
        call set_xs_data(nnod, ng, ind_mat)
        call xsec_setup(fdm % xs, sigtr, siga, nuf, sigf, sigs, chi, dc, bxtab, bcrod)
        call xs_time % off

        if (allocated(th)) then
            call th_init()
        end if

        if (allocated(tr)) then
            allocate(total_beta(nmat))
            allocate(dfis(nnod))
            call set_fdm_transient(fdm, total_beta, dfis)
            call set_transient_data(tr, ounit, nf, nnod, ng, nmat, total_time, time_step_1, time_step_2, &
            time_mid, small_theta, bextr, bxtab)
            if (bxtab == NO) then
                do i_mat = 1, nmat
                    lambda_tmp   (i_mat,:) = lambda
                    beta_tmp     (i_mat,:) = beta
                    neut_velo_tmp(i_mat,:) = neutron_velo
                end do
            else
                do i_mat = 1, nmat
                    lambda_tmp   (i_mat,:) = m(i_mat) % lambda
                    beta_tmp     (i_mat,:) = m(i_mat) % beta
                    neut_velo_tmp(i_mat,:) = m(i_mat) % velo
                end do
            end if

            call set_transient_pointer(tr, fdm, xsc, beta_tmp, lambda_tmp, neut_velo_tmp, total_beta, dfis, ind_mat)
        end if

    end subroutine

end module
    