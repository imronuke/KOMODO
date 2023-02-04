module control
    use iso_fortran_env, only: real64, output_unit, error_unit
    use data
    use print
    use fdm
    use xsec
    use time
    use th
    use transient
    use utilities

    implicit none

    private

    save

    integer, parameter                    :: dp = real64
    type(fdm_rect), public               :: fdm            ! Finite Difference Object
    type(th_type), allocatable, public    :: th             ! Thermal-hydraulics Object
    type(trans_type), allocatable, public :: tr             ! Transient Object
    type(xs_change_type), public          :: xsc            ! object for xs changes due to a parameter change

    real(dp), allocatable      :: mesh_temp(:,:) ! fuel mesh temperature[nzz, n_fuel]

    public :: forward, adjoint, fixedsrc, cbc_search, rod_eject

    contains

    !===============================================================================================!
    ! Transient / Rod Ejection calculation
    !===============================================================================================!

    subroutine rod_eject()

        real(dp) :: fn(nnod)

        call calculation_init()

        if (bther == NO) then
            call rod_eject_no_th(tr)
        ! else
        !     call cbc_search_no_th()
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

        call calculation_init()
        call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bank_pos)
        
        call print_head()

        call outer_iter(fdm, print_iter = .true., is_converge=converge)
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

        integer  :: n
        real(dp) :: K1, K2
        real(dp) :: bc1, bc2
        logical  :: converge

        call calculation_init()

        call print_head()

        ! First try
        bcon = xsc % bcon % ref
        call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bank_pos)
        call outer_iter(fdm, print_iter = .false., is_converge=converge)
        if (.not. converge) call print_fail_converge()
        bc1 = bcon
        K1  = fdm % Keff
        write(output_unit,101) 1, bcon, K1, fdm % flux_diff, fdm % fsrc_diff

        ! Second try
        bcon = bcon + (K1 - 1.) * bcon   ! Guess next critical boron concentration
        call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bank_pos)
        call outer_iter(fdm, print_iter = .false., is_converge=converge)
        if (.not. converge) call print_fail_converge()
        bc2 = bcon
        K2  = fdm % Keff
        write(output_unit,101) 2, bcon, K2, fdm % flux_diff, fdm % fsrc_diff

        n = 2
        do
            n = n + 1
            bcon = bc2 + (1. - K2) / (K1 - K2) * (bc1 - bc2)
            call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bank_pos)
            call outer_iter(fdm, print_iter = .false., is_converge=converge)
            if (.not. converge) call print_fail_converge()
            bc1 = bc2
            bc2 = bcon
            K1  = K2
            K2  = fdm % Keff
            write(output_unit,101) n, bcon, K2, fdm % flux_diff, fdm % fsrc_diff
            if (abs(bc1-bc2) < 0.1) exit
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

        call calculation_init()

        call print_head()

        ! First try
        bcon = xsc % bcon % ref
        call th_iteration(2, ftem_diff)
        bc1 = bcon
        K1  = fdm % Keff
        write(output_unit,102) 1, bcon, K1, fdm % flux_diff, fdm % fsrc_diff, ftem_diff

        ! Second try
        bcon = bcon + (K1 - 1.) * bcon   ! Guess next critical boron concentration
        call th_iteration(2, ftem_diff)
        bc2 = bcon
        K2  = fdm % Keff
        write(output_unit,102) 2, bcon, K2, fdm % flux_diff, fdm % fsrc_diff, ftem_diff

        n = 2
        do
            n = n + 1
            bcon = bc2 + (1. - K2) / (K1 - K2) * (bc1 - bc2)
            call th_iteration(2, ftem_diff)
            bc1 = bc2
            bc2 = bcon
            K1  = K2
            K2  = fdm % Keff
            write(output_unit,102) n, bcon, K2, fdm % flux_diff, fdm % fsrc_diff, ftem_diff
            if (abs(bc1-bc2) < 0.1) exit
        end do

        call print_tail(fdm)

        102 format(I3, F9.2, F14.5, ES14.5, ES13.5, ES17.5)

    end subroutine

    !===============================================================================================!
    ! thermal-hydraulics iteration
    !===============================================================================================!

    subroutine th_iteration(max_iter, ftem_diff)

        integer, intent(in)    :: max_iter
        real(dp), intent(out)  :: ftem_diff
        
        integer  :: i
        real(dp) :: ftem_prev(nnod)

        do i = 1, max_iter
            ftem_prev = ftem
            call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bank_pos)
            call th_outer(fdm)
            call th_solver(flux, ftem, mtem, cden)
            ftem_diff = maxval(abs(ftem - ftem_prev))
        end do

    end subroutine

    !===============================================================================================!
    ! thermal-hydraulics solver
    !===============================================================================================!

    subroutine th_solver(f0, fuel_temp, mod_temp, cool_dens)

        real(dp), intent(in)   :: f0(:,:)  ! Flux
        real(dp), intent(out)  :: fuel_temp(:)
        real(dp), intent(out)  :: mod_temp(:)
        real(dp), intent(out)  :: cool_dens(:)

        real(dp)  :: linear_pow(nnod)          ! node-wise inear Power
        real(dp)  :: lin_power(nzz)            ! sub-channel axial power
        real(dp)  :: heatf(nzz)
        real(dp)  :: ft(nzz)
        real(dp)  :: mt(nzz)
        real(dp)  :: cd(nzz)
        
        integer :: i, j, k, n

        if (.not. allocated(th)) return

        call th_time % on

        call get_lin_power(f0, linear_pow)

        do j = 1, nyy
            do i = ystag(j) % smin, ystag(j) % smax

                do k = 1, nzz
                    n = xyz(i, j, k)
                    lin_power(k)    = linear_pow(n)
                    heatf(k)        = heat_flux(n)
                    mesh_temp(k, :) = tfm(n, :)
                end do

                call th_update(th, lin_power, mesh_temp, heatf, ft, mt, cd)

                do k = 1, nzz
                    n = xyz(i, j, k)
                    heat_flux(n) = heatf(k)
                    tfm(n, :)    = mesh_temp(k, :)
                    fuel_temp(n) = ft(k)
                    mod_temp(n)  = mt(k)
                    cool_dens(n) = cd(k)
                end do

            end do
        end do

        call th_time % off

    end subroutine

    !===============================================================================================!
    ! calculate linear power density
    !===============================================================================================!

    subroutine get_lin_power(f0, lin_power)

        real(dp), intent(in)  :: f0(:,:)  ! Flux
        real(dp), intent(out) :: lin_power(:)

        real(dp) :: pdist(nnod)
        integer  :: n

        call get_power_dist(f0, pdist)

        do n = 1, nnod
            lin_power(n) = pdist(n) * power * percent_pow * 0.01 &
            / (node_nf(ix(n),iy(n)) * zdel(iz(n)))
        end do   

    end subroutine

    !===============================================================================================!
    ! calculate power distribution (normalize to 1.0)
    !===============================================================================================!

    subroutine get_power_dist(f0, pdist)

        real(dp), intent(in) :: f0(:,:)  ! Flux
        real(dp), intent(out) :: pdist(:)

        integer :: g, n
        real(dp) :: pow, powtot
        
        pdist = 0._dp
        do g= 1, ng
            do n= 1, nnod
              pow = f0(n,g) * fdm % xs % sigf(n,g)
              if (pow < 0.) pow = 0.
              pdist(n) = pdist(n) + pow*vdel(n)
            end do
        end do
        
        ! Normalize to 1.0
        powtot = sum(pdist)
        do n = 1, nnod
            pdist(n) = pdist(n) / powtot
        end do

    end subroutine

    !===============================================================================================!
    ! Initialize thermal-hydraulics solver
    !===============================================================================================!

    subroutine th_init(calc_mode)

        character(*), optional :: calc_mode

        if (.not. allocated(th)) return
        
        ! set data in the th object
        call set_th_data(th, ounit, nzz, zdel, rf, tg, tc, ppitch, cflow, cf, t_inlet)

        ! Allocate necessary data
        allocate(mesh_temp(nzz, th % n_ring+1))
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

        if (present(calc_mode)) then
            if (calc_mode == 'transient') then
                allocate(enthalpy(nnod))
                allocate(flow_rate(nnod))
            end if
        end if

    end subroutine

    !===============================================================================================!
    ! Initialize calculation
    !===============================================================================================!

    subroutine calculation_init()

        integer :: nodal_interval, max_outer

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

        nodal_interval = ceiling((nxx + nyy + nzz) / 2.5)
        if (allocated(th)) then
            max_outer = 20
        else
            max_outer = 1000
        end if

        call fdm_time % on
        call set_fdm_data(fdm, ounit, kern, ng, nnod, nxx, nyy, nzz, nmat, &
        east, west, north, south, top, bottom, nodal_interval, 1.e-5_dp, 1.e-5_dp, max_outer, 2, 5)
        call set_fdm_pointer(fdm, flux, fsrc, exsrc, ind_mat, xdel, ydel, zdel, vdel, &
        xstag, ystag, ix, iy, iz, xyz)
        call fdm_time % off

        call xs_time % on
        call set_xs_data(nnod, ng, ind_mat)
        call xsec_setup(fdm % xs, sigtr, siga, nuf, sigf, sigs, chi, dc)
        call xs_time % off

        if (allocated(th)) then
            if (allocated(tr)) then
                call th_init('transient')
            else
                call th_init()
            end if
        end if

        if (allocated(tr)) then
            allocate(total_beta(nmat))
            allocate(dfis(nnod))
            call set_fdm_transient(fdm, total_beta, dfis)
            call set_transient_data(tr, ounit, nf, nnod, ng, nmat, total_time, time_step_1, time_step_2, &
            time_mid, small_theta, bxtab, bextr)
            call set_transient_pointer(tr, fdm, xsc, beta, lambda, neutron_velo,total_beta, dfis, ind_mat)
        end if

    end subroutine

end module
    