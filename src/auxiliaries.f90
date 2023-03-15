module auxiliaries

    ! Module auxiliaries

    use iso_fortran_env, only: real64, error_unit, output_unit
    use data
    use nodal
    use xsec
    use fdm
    use th
    use time

    implicit none

    private

    save

    integer, parameter  :: dp = real64

    type(fdm_rect), public, pointer          :: fdm
    type(th_type), pointer           :: th
    type(xs_change_type), pointer    :: xsc

    public :: adjust_keff, get_leakage, integrate, xs_update, move_crod, get_adjoint_flux
    public :: set_aux_pointer, th_iteration, th_solver_transient
    public :: par_ave, par_ave_out, print_fail_converge

    contains

    !===============================================================================================!
    ! To adjuts the Keff to 1.0 if it is not equal to 1.0
    ! When starting transient, normally reactor at steady state condition (k = 1.0).
    ! In case the k is not equal to 1.0, then it is adjusted to 1.0
    !===============================================================================================!

    subroutine adjust_keff(c, xsc)

        class(fdm_rect)       :: c
        class(xs_change_type) :: xsc
        
        integer :: i
        logical :: converge
        
        write(output_unit,'(A)', advance='no') '  steady state k-eff not equal to one, force it to one ... '
        write(ounit, *)
        write(ounit, '(A46,F9.6)') '  INITIAL MULTIPLICATION EFFECTIVE (K-EFF) = ',  c % Keff
        write(ounit, *) '  WARNING: THE STEADY STATE K-EFF IS NOT EQUAL TO 1.0'
        write(ounit, *) '  AND NOW IT IS FORCED TO 1.0 BY MODIFYING THE nu*sigf CROSS SECTIONS '
        write(ounit, *)
        do i = 1, 10
            ! Adjust material-wise nu*fission and xs changes due to CR
            nuf = nuf / c % Keff
            if (allocated(xsc % crod)) xsc % crod % nuf = xsc % crod % nuf / c % Keff

            call xsec_update(c % xs, xsc, bcon, ftem, mtem, cden, bank_pos)
            call outer_iter(c, print_iter = .false., is_converge=converge)
            if (.not. converge) call print_fail_converge()
            if (abs(c % Keff-1.0) < 1.e-5) exit
        end do
        if (i == 10) then
            write(error_unit,*) "K-EFF STILL NOT EQUAL TO ONE. KOMODO IS STOPPING"
            stop
        end if

        write(output_unit,*) 'done'
        
    end subroutine

    !===============================================================================================!
    ! To calculate adjoint flux
    !===============================================================================================!

    subroutine get_adjoint_flux(xsc, adj_flux)
        
        class(xs_change_type), intent(in) :: xsc
        real(dp), intent(out)             :: adj_flux(:,:)

        type(fdm_rect), allocatable :: fdm_adj
        real(dp)                     :: flux_tmp(nnod, ng), fsrc_tmp(nnod)
        logical                      :: converge

        allocate(fdm_adj)

        call set_fdm_data(fdm_adj, ounit, kern, ng, nnod, nxx, nyy, nzz, nmat, &
        east, west, north, south, top, bottom, nodal_interval, max_flux_error, max_fsrc_error, &
        max_outer, max_outer_th, max_inner, extrp_interval)

        call set_fdm_pointer(fdm_adj, flux_tmp, fsrc_tmp, exsrc, ind_mat, xdel, ydel, zdel, vdel, &
        xstag, ystag, ix, iy, iz, xyz)

        call xsec_setup(fdm_adj % xs, sigtr, siga, nuf, sigf, sigs, chi, dc, bxtab, bcrod)
        call xsec_update(fdm_adj % xs, xsc, bcon, ftem, mtem, cden, bank_pos)

        call outer_adjoint(fdm_adj, print_iter = .true., is_converge=converge)
        if (.not. converge) call print_fail_converge()

        adj_flux = flux_tmp

        call deallocate_fdm(fdm_adj)
 
    end subroutine

    !===============================================================================================!
    ! To calculate leakage in x, y and z directions
    !===============================================================================================!

    subroutine get_leakage(c, n, g, Lx, Ly, Lz)
        
        class(fdm_rect)        :: c
        integer, intent(in)    :: n, g
        real(dp), intent(out)  :: Lx, Ly, Lz

        real(dp) :: jp, jm
        integer  :: p, m
        integer  :: i, j, k

                
        ! set i, j, k
        i = c % ix(n); j = c % iy(n); k = c % iz(n)
        
        ! x-direction zeroth transverse leakage
        if (i /= c % ystag(j) % smax) p = c % xyz(i+1,j,k)
        if (i /= c % ystag(j) % smin) m = c % xyz(i-1,j,k)

        ! west
        if (i == c % ystag(j) % smax) then
            if (east == reflective) then
                jp = 0.0
            else
                jp = c % df(n,g,1) * c % flux(n,g) - c % dn(n,g,1) * c % flux(n,g)
            end if
        else
            jp = - c % df(n,g,1) * (c % flux(p,g) - c % flux(n,g)) &
                 - c % dn(n,g,1) * (c % flux(p,g) + c % flux(n,g))
        end if

        ! east
        if (i == c % ystag(j) % smin) then
            if (west == reflective) then
                jm = 0.0
            else
                jm = - c % df(n,g,2) * c % flux(n,g) - c % dn(n,g,2) * c % flux(n,g)
            end if
        else
            jm = - c % df(n,g,2) * (c % flux(n,g) - c % flux(m,g)) &
                 - c % dn(n,g,2) * (c % flux(n,g) + c % flux(m,g))
        end if
        
        Lx = (jp - jm)  / c % xdel(i)
        
        ! y-direction zeroth transverse leakage
        if (j /= c % xstag(i) % smax) p = c % xyz(i,j+1,k)
        if (j /= c % xstag(i) % smin) m = c % xyz(i,j-1,k)

        ! north
        if (j == c % xstag(i) % smax) then
            if (north == reflective) then
                jp = 0.0
            else
                jp = c % df(n,g,3) * c % flux(n,g) - c % dn(n,g,3) * c % flux(n,g)
            end if
        else
            jp = - c % df(n,g,3) * (c % flux(p,g) - c % flux(n,g)) &
                 - c % dn(n,g,3) * (c % flux(p,g) + c % flux(n,g))
        end if

        ! south
        if (j == c % xstag(i) % smin) then
            if (south == reflective) then
                jm = 0.0
            else
                jm = - c % df(n,g,4) * c % flux(n,g) - c % dn(n,g,4) * c % flux(n,g)
            end if
        else
            jm = - c % df(n,g,4) * (c % flux(n,g) - c % flux(m,g)) &
                 - c % dn(n,g,4) * (c % flux(n,g) + c % flux(m,g))
        end if
        
        Ly = (jp - jm)  / c % ydel(j)
        
        ! z-direction zeroth transverse leakage
        if (k /= c % nzz) p = c % xyz(i,j,k+1)
        if (k /= 1      ) m = c % xyz(i,j,k-1)

        ! top
        if (k == c % nzz) then
            if (top == reflective) then
                jp = 0.0
            else
                jp = c % df(n,g,5) * c % flux(n,g) - c % dn(n,g,5) * c % flux(n,g)
            end if
        else
            jp = - c % df(n,g,5) * (c % flux(p,g) - c % flux(n,g)) &
                 - c % dn(n,g,5) * (c % flux(p,g) + c % flux(n,g))
        end if

        ! bottom
        if (k == 1) then
            if (bottom == reflective) then
                jm = 0.0
            else
                jm = - c % df(n,g,6) * c % flux(n,g) - c % dn(n,g,6) * c % flux(n,g)
            end if
        else
            jm = - c % df(n,g,6) * (c % flux(n,g) - c % flux(m,g)) &
                 - c % dn(n,g,6) * (c % flux(n,g) + c % flux(m,g))
        end if

        Lz = (jp - jm)  / c % zdel(k)
        
 
    end subroutine

    !===============================================================================================!
    ! To update xs
    !===============================================================================================!

    subroutine xs_update(xs, xsc)

        type(xs_rect)         :: xs
        class(xs_change_type) :: xsc
        
        call xsec_update(xs, xsc, bcon, ftem, mtem, cden, bank_pos)
        
    end subroutine

    !===============================================================================================!
    ! To update xs
    !===============================================================================================!

    subroutine move_crod(t_step, t_next)

        real(dp), intent(in) :: t_step, t_next

        integer :: n

        do n = 1, nbank
            if (direction(n) == 1) then   ! If CR moving down
                if (t_next-t_move(n) > 1.e-5_dp .AND. bank_pos_final(n)-bank_pos(n) < 1.e-5_dp) then
                    bank_pos(n) = bank_pos(n) - t_step*bank_speed(n)
                    if (bank_pos(n) < bank_pos_final(n)) bank_pos(n) = bank_pos_final(n)  ! If bpos exceed, set to fbpos
                end if
            else if (direction(n) == 2) then ! If CR moving up
                if (t_next-t_move(n) > 1.e-5_dp .AND. bank_pos_final(n)-bank_pos(n) > 1.e-5_dp) then
                    bank_pos(n) = bank_pos(n) + t_step*bank_speed(n)
                    if (bank_pos(n) > bank_pos_final(n)) bank_pos(n) = bank_pos_final(n)  ! If bpos exceed, set to fbpos
                end if
            end if
         end do
        
    end subroutine

    !===============================================================================================!
    ! do volume integration on variable a
    !===============================================================================================!

    pure function integrate(a) result(res)

        real(dp), intent(in)          :: a(:)
        real(dp)                      :: res

        integer                 :: n

        res = 0.
        do n = 1, nnod
            res = res + a(n) * vdel(n)
        end do
    
    end function

    !===============================================================================================!
    ! set auxiliaries pointer
    !===============================================================================================!

    subroutine set_aux_pointer(f, t, x)

        type(fdm_rect), target, intent(in)       :: f
        type(th_type), target, intent(in)        :: t
        type(xs_change_type), target, intent(in) :: x

        fdm => f
        th => t
        xsc => x

    end subroutine

    !===============================================================================================!
    ! thermal-hydraulics iteration
    !===============================================================================================!

    subroutine th_iteration(max_iter, print_iter, fuel_temp_err, is_converge)

        integer, intent(in)    :: max_iter
        logical, intent(in)    :: print_iter
        real(dp), intent(out)  :: fuel_temp_err
        logical, intent(out)   :: is_converge
        
        integer  :: i
        real(dp) :: ftem_prev(nnod)

        is_converge = .false.
        do i = 1, max_iter
            ftem_prev = ftem
            call xsec_update(fdm % xs, xsc, bcon, ftem, mtem, cden, bank_pos)
            call th_outer(fdm)
            call th_solver(enthalpy, ftem, mtem, cden)
            fuel_temp_err = maxval(abs(ftem - ftem_prev))
            if (print_iter) then
                write(output_unit,'(I5,F13.5,3ES15.5)') i, fdm % Keff, fdm % fsrc_diff, fdm % flux_diff, fuel_temp_err
                write(ounit,'(I5,F13.5,3ES15.5)') i, fdm % Keff, fdm % fsrc_diff, fdm % flux_diff, fuel_temp_err
            end if
            if ((fdm % max_flux_err > fdm % flux_diff) .and. (fdm % max_fiss_err > fdm % fsrc_diff)) then
                is_converge = .true.
                exit
            end if
        end do

    end subroutine


    !===============================================================================================!
    ! thermal-hydraulics solver
    !===============================================================================================!

    subroutine th_solver(cool_ent, fuel_temp, mod_temp, cool_dens)

        real(dp), intent(out)      :: cool_ent(:)
        real(dp), intent(out)      :: fuel_temp(:)
        real(dp), intent(out)      :: mod_temp(:)
        real(dp), intent(out)      :: cool_dens(:)

        real(dp)  :: linear_pow(nnod)          ! node-wise inear Power
        real(dp)  :: lin_power(nzz)            ! sub-channel axial power
        real(dp)  :: heatf(nzz)
        real(dp)  :: ent(nzz)
        real(dp)  :: ft(nzz)
        real(dp)  :: mt(nzz)
        real(dp)  :: cd(nzz)
        real(dp)  :: mesh_temp(nzz, th % n_ring+1) ! fuel mesh temperature[nzz, n_fuel]
        
        integer :: i, j, k, n

        call th_time % on

        call get_lin_power(fdm % flux, percent_pow*0.01, linear_pow)

        do j = 1, nyy
            do i = ystag(j) % smin, ystag(j) % smax

                do k = 1, nzz
                    n = xyz(i, j, k)
                    lin_power(k)    = linear_pow(n)
                    heatf(k)        = heat_flux(n)
                    mesh_temp(k, :) = tfm(n, :)
                end do

                call th_update(th, lin_power, mesh_temp, heatf, ent, ft, mt, cd)

                do k = 1, nzz
                    n = xyz(i, j, k)
                    heat_flux(n) = heatf(k)
                    tfm(n, :)    = mesh_temp(k, :)
                    cool_ent(n)  = ent(k)
                    fuel_temp(n) = ft(k)
                    mod_temp(n)  = mt(k)
                    cool_dens(n) = cd(k)
                end do

            end do
        end do

        call th_time % off

    end subroutine

    !===============================================================================================!
    ! thermal-hydraulics transient solver
    ! This routine updates ftem, mtem and cden implicitly
    !===============================================================================================!

    subroutine th_solver_transient(power_level, t_step)

        real(dp), intent(in)       :: power_level
        real(dp), intent(in)       :: t_step

        real(dp)  :: linear_pow(nnod)          ! node-wise inear Power
        real(dp)  :: lin_power(nzz)            ! sub-channel axial power
        real(dp)  :: heatf(nzz)
        real(dp)  :: ft(nzz)
        real(dp)  :: mt(nzz)
        real(dp)  :: cd(nzz)
        real(dp)  :: ent(nzz)
        real(dp)  :: mesh_temp(nzz, th % n_ring+1) ! fuel mesh temperature[nzz, n_fuel]
        
        integer :: i, j, k, n

        call th_time % on

        call get_lin_power(fdm % flux, power_level, linear_pow)

        do j = 1, nyy
            do i = ystag(j) % smin, ystag(j) % smax

                do k = 1, nzz
                    n = xyz(i, j, k)
                    lin_power(k)    = linear_pow(n)
                    heatf(k)        = heat_flux(n)
                    ent(k)          = enthalpy(n)
                    ft(k)           = ftem(n)
                    mt(k)           = mtem(n)
                    cd(k)           = cden(n)
                    mesh_temp(k, :) = tfm(n, :)
                end do

                call th_transient(th, t_step, lin_power, flow_rate, ent, mesh_temp, heatf, ft, mt, cd)

                do k = 1, nzz
                    n = xyz(i, j, k)
                    heat_flux(n) = heatf(k)
                    tfm(n, :)    = mesh_temp(k, :)
                    ftem(n)      = ft(k)
                    mtem(n)      = mt(k)
                    cden(n)      = cd(k)
                end do

            end do
        end do

        call th_time % off

    end subroutine


    !===============================================================================================!
    ! calculate linear power density
    !===============================================================================================!

    subroutine get_lin_power(f0, power_level, lin_power)

        real(dp), intent(in)       :: f0(:,:)  ! Flux
        real(dp)                   :: power_level
        real(dp), intent(out)      :: lin_power(:)

        real(dp) :: pdist(nnod)
        integer  :: n

        call get_power_dist(f0, pdist)

        do n = 1, nnod
            lin_power(n) = pdist(n) * power * power_level &
            / (node_nf(ix(n),iy(n)) * zdel(iz(n)))
        end do   

    end subroutine

    !===============================================================================================!
    ! calculate power distribution (normalize to 1.0)
    !===============================================================================================!

    subroutine get_power_dist(f0, pdist)

        real(dp), intent(in)       :: f0(:,:)  ! Flux
        real(dp), intent(out)      :: pdist(:)

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
    ! To calculate average fuel temp (only for active core)
    !===============================================================================================!

    subroutine par_ave(c, par, ave)
        
        class(fdm_base)      :: c
        real(dp), intent(in)  :: par(:)
        real(dp), intent(out) :: ave
        
        real(dp) :: dum, dum2
        integer :: n
        
        dum = 0.; dum2 = 0.
        do n = 1, nnod
           if (c % xs % nuf(n,ng) > 0.) then
              dum = dum + par(n) * vdel(n)
              dum2 = dum2 + vdel(n)
           end if
        end do
        
        ave = dum / dum2
        
    end subroutine
        
    !===============================================================================================!
    ! To calculate average outlet coolant temperature
    !===============================================================================================!
        
    subroutine par_ave_out(c, par, ave)
        
        class(fdm_base)      :: c
        real(dp), intent(in)  :: par(:)
        real(dp), intent(out) :: ave
 
        real(dp) :: dum, dum2
        integer, DIMENSION(nxx,nyy) :: zmax
        integer :: n, i, j, k
        
        ! get number of nodex in axial direction from bottom -> fuel
        do j = 1, nyy
          do i = ystag(j)%smin, ystag(j)%smax
            zmax(i,j) = 0
            do k = 1, nzz/2
              if (c % xs % nuf(xyz(i,j,k),ng) < 1.e-5_dp) zmax(i,j) = zmax(i,j) + 1
            end do
          end do
        end do
        ! get number of nodex in axial direction from fuel -> top reflectors
        do j = 1, nyy
          do i = ystag(j)%smin, ystag(j)%smax
            do k = 1, nzz
              if (c % xs % nuf(xyz(i,j,k),ng) > 1.e-5) zmax(i,j) = zmax(i,j) + 1
            end do
          end do
        end do
        
        dum = 0.; dum2 = 0.
        do n = 1, nnod
           IF (iz(n) == zmax(ix(n),iy(n)) .AND. c % xs % nuf(n,ng) > 1.e-5) THEN
              dum = dum + par(n) * vdel(n)
              dum2 = dum2 + vdel(n)
           end IF
        end do
        
        ave = dum / dum2
        
    end subroutine

    !===============================================================================================!
    ! Give message if iterations fail to converge
    !===============================================================================================!

    subroutine print_fail_converge()

        write(error_unit,*)
        write(error_unit,*) '  MAXIMUM NUMBER OF OUTER ITERATION IS REACHED'
        write(error_unit,*) '  CHECK PROBLEM SPECIFICATION OR CHANGE ITERATION CONTROL (%ITER).'
        write(error_unit,*) '   The two-node nonlinear iteration seems not stable. Try these:'
        write(error_unit,*) '   1. Change iteration control using %ITER card, '
        write(error_unit,*) '      Perhaps by making nodal update less frequent,'
        write(error_unit,*) '      or increase number inner iteration per outer iteration.'
        write(error_unit,*) '   2. Ensure the node sizes in all directions are as uniform as possible. '
        write(error_unit,*) '      Also try smaller node size.'
        write(error_unit,*) '   3. For transient problem, try to reduce time step size.'
        write(error_unit,*) 
        write(error_unit,*) '   If this error persists, contact me at makrus.imron@gmail.com'
        write(error_unit,*) '  KOMODO is stoping'
        stop

    end subroutine

end module