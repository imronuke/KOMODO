module control
    use iso_fortran_env, only: output_unit, error_unit
    use data
    use print
    use fdm
    use xsec
    use time
    use th
    use utilities

    implicit none

    private

    save

    real(dp), allocatable      :: mesh_temp(:,:) ! fuel mesh temperature[nzz, n_fuel]

    public :: forward, adjoint, fixedsrc, cbsearch, cbsearch_th

    contains

    !===============================================================================================!
    ! forward calculation
    !===============================================================================================!

    subroutine forward()

        integer  :: g, n
        real(dp) :: fn(nnod)
        logical  :: converge

        call calculation_init()
        call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bpos)
        
        call print_head()
        call outer_iter(fdm, 1.e-5_dp, 1.e-5_dp, max_outer=1000, max_inner=2, &
        extrp = 5, print_iter = .true., is_converge=converge)

        if (.not. converge) call print_fail_converge()

        fn = 0
        do g = 1, ng
            do n = 1, nnod
                fn(n) = fn(n) + fdm % xs % nuf(n,g) * flux(n,g)
            end do
        end do

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
        call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bpos)

        call print_head()
        call outer_fixed_src(fdm, 1.e-5_dp, 1.e-5_dp, max_outer=1000, max_inner=2, &
        extrp = 5, print_iter = .true., is_converge=converge)

        if (.not. converge) call print_fail_converge()

        fn = 0
        do g = 1, ng
            do n = 1, nnod
                fn(n) = fn(n) + fdm % xs % nuf(n,g) * flux(n,g)
            end do
        end do

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
        call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bpos)

        call print_head()
        call outer_adjoint(fdm, 1.e-5_dp, 1.e-5_dp, max_outer=1000, max_inner=2, &
        extrp = 5, print_iter = .true., is_converge=converge)

        if (.not. converge) call print_fail_converge()

        fn = 0
        do g = 1, ng
            do n = 1, nnod
                fn(n) = fn(n) + fdm % xs % nuf(n,g) * flux(n,g)
            end do
        end do

        call print_power_map(fn)

    end subroutine

    !===============================================================================================!
    ! critical boron search calculation
    !===============================================================================================!

    subroutine cbsearch()

        integer  :: n
        real(dp) :: K1, K2
        real(dp) :: bc1, bc2
        logical  :: converge

        call calculation_init()

        call print_head()

        ! First try
        bcon = xsc % bcon % ref
        call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bpos)
        call outer_iter(fdm, 1.e-5_dp, 1.e-5_dp, max_outer=1000, max_inner=2, &
        extrp = 5, print_iter = .false., is_converge=converge)
        if (.not. converge) call print_fail_converge()
        bc1 = bcon
        K1  = fdm % Keff
        write(output_unit,101) 1, bcon, K1, fdm % flux_diff, fdm % fsrc_diff

        ! Second try
        bcon = bcon + (K1 - 1.) * bcon   ! Guess next critical boron concentration
        call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bpos)
        call outer_iter(fdm, 1.e-5_dp, 1.e-5_dp, max_outer=1000, max_inner=2, &
        extrp = 5, print_iter = .false., is_converge=converge)
        if (.not. converge) call print_fail_converge()
        bc2 = bcon
        K2  = fdm % Keff
        write(output_unit,101) 2, bcon, K2, fdm % flux_diff, fdm % fsrc_diff

        n = 2
        do
            n = n + 1
            bcon = bc2 + (1. - K2) / (K1 - K2) * (bc1 - bc2)
            call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bpos)
            call outer_iter(fdm, 1.e-5_dp, 1.e-5_dp, max_outer=1000, max_inner=2, &
            extrp = 5, print_iter = .false., is_converge=converge)
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
    ! critical boron search calculation
    !===============================================================================================!

    subroutine cbsearch_th()

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

        call print_tail()

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
            call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bpos)
            call fixed_outer(fdm, max_outer=20, max_inner=2, extrp=5)
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

        if (calc_mode == 'transient') then
            allocate(enthalpy(nnod))
            allocate(flow_rate(nnod))
        end if

    end subroutine

    !===============================================================================================!
    ! Initialize calculation
    !===============================================================================================!

    subroutine calculation_init()

        call fdm_time % on

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

        call set_rect_data(fdm, kern, ng, nnod, nxx, nyy, nzz, nmat, &
        east, west, north, south, top, bottom, ounit)
        
        call set_rect_pointer(fdm, flux, fsrc, exsrc, ind_mat, xdel, ydel, zdel, vdel, &
        xstag, ystag, ix, iy, iz, xyz)

        call fdm_time % on

        call xs_time % on

        call set_xs_data(nnod, ng, ind_mat)
        call xsec_setup(fdm % xs, sigtr, siga, nuf, sigf, sigs, chi, dc)

        call xs_time % off

        if (allocated(th)) then
            call th_init()
        else
            return
        end if

    end subroutine

    !===============================================================================================!
    ! Initialize calculation
    !===============================================================================================!

    subroutine print_fail_converge()

        write(error_unit,*)
        write(error_unit,*) '  MAXIMUM NUMBER OF OUTER ITERATION IS REACHED IN FORWARD CALCULATION.'
        write(error_unit,*) '  CHECK PROBLEM SPECIFICATION OR CHANGE ITERATION CONTROL (%ITER).'
        write(error_unit,*) '  PERHAPS BY MAKING FISSION SOURCE INTERPOLATION MORE FREQUENT'
        write(error_unit,*) '  KOMODO IS STOPING...'
        stop

    end subroutine

end module
    