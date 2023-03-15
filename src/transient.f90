module transient

    use iso_fortran_env, only: real64, error_unit, output_unit
    use xsec, only: xs_rect
    use auxiliaries
    use fdm
    use xsec, only: xs_change_type
    use print, only: print_fail_converge, print_transient_step, print_head

    implicit none

    private

    save

    integer, parameter :: dp = real64
    integer            :: fid

    type, public :: trans_type
        integer                    :: nf                    ! Number of delayed neutron precusor family
        integer                    :: nnod                  ! Total number of nodes
        integer                    :: nmat                  ! Number of materials
        integer                    :: ng                    ! number of energy group
        real(dp)                   :: theta = 1._dp         ! Small theta (default is fully implicit)
        real(dp)                   :: big_theta = 0._dp     ! big theta for transient using theta method
        real(dp), allocatable      :: omega(:,:)            ! Exponential transformation constant

        real(dp), allocatable      :: beta(:,:)                ! beta (delayed neutron fraction)
        real(dp), allocatable      :: lambda(:,:)              ! lambda (precursor decay constant)
        real(dp)                   :: core_beta                ! Core averaged beta
        real(dp), allocatable      :: neutron_velo(:,:)        ! Neutron velocity
        real(dp), allocatable      :: precursor(:,:)           ! precursor density
        real(dp), allocatable      :: L(:,:)                   ! Total leakages for node n and group g
        real(dp), pointer          :: total_beta(:) => null()  ! total beta for each material
        real(dp), pointer          :: dfis(:) => null()              

        real(dp)                   :: total_time            ! TOTAL SIMULATION TIME
        real(dp)                   :: t_step_1              ! FIRST TIME STEP
        real(dp)                   :: t_step_2              ! SECOND TIME STEP
        real(dp)                   :: time_mid              ! WHEN SECOND TIME STEP START

        real(dp), allocatable      :: prev_flux(:,:)        ! previous time-step nodes flux
        real(dp), allocatable      :: prev_fsrc(:)          ! previous time-step nodes fission source
        real(dp), allocatable      :: adj_flux(:,:)         ! Adjoint flux
        logical                    :: is_exp                ! is exponential transformation used
        logical                    :: is_tabular            ! is tabular xs used

        type(fdm_rect), pointer       :: f                     ! FDM object
        type(xs_rect), pointer        :: xs                    ! node-wise xs
        type(xs_change_type), pointer :: xsc                   ! xs changes
        integer, pointer              :: ind_mat(:)            ! material index
    end type

    public :: set_transient_data, set_transient_pointer
    public :: rod_eject_transient

    contains

    !===============================================================================================!
    ! set transient data
    !===============================================================================================!

    subroutine set_transient_data(tr, fid_inp, nf, nnod, ng, nmat, total_time, time_step_1, time_step_2, &
        time_mid, theta, bextr, bxtab)

        class(trans_type)       :: tr
        integer, intent(in)     :: fid_inp         ! number of energy group
        integer, intent(in)     :: nf              ! number delayed neutron precursor group family
        integer, intent(in)     :: ng              ! number of energy group
        integer, intent(in)     :: nnod            ! Total number of nodes
        integer, intent(in)     :: nmat            ! Number of materials
        real(dp), intent(in)    :: total_time      ! TOTAL SIMULATION TIME
        real(dp), intent(in)    :: time_step_1     ! FIRST TIME STEP
        real(dp), intent(in)    :: time_step_2     ! SECOND TIME STEP
        real(dp), intent(in)    :: time_mid        ! WHEN SECOND TIME STEP START
        real(dp), intent(in)    :: theta           ! Small theta
        integer, intent(in)     :: bextr           ! is exponential transformation used
        integer, intent(in)     :: bxtab           ! is tabular xs is used

        fid       = fid_inp
        tr % nf   = nf
        tr % ng   = ng
        tr % nmat = nmat
        tr % nnod = nnod
        tr % total_time = total_time
        tr % t_step_1 = time_step_1
        tr % t_step_2 = time_step_2
        tr % time_mid    = time_mid
        tr % theta = theta
        tr % big_theta = (1._dp - theta) / theta
        if (bextr == 1) then
            tr % is_exp = .true.
        else
            tr % is_exp = .false.
        end if
        if (bxtab == 1) then
            tr % is_tabular = .true.
        else
            tr % is_tabular = .false.
        end if

        allocate(tr % prev_flux(nnod,ng))
        allocate(tr % prev_fsrc(nnod))
        allocate(tr % omega(nnod,ng))
        allocate(tr % precursor(nnod,nf))
        allocate(tr % L(nnod,ng))
        allocate(tr % adj_flux(nnod,ng))
        allocate(tr % beta(nmat, nf))
        allocate(tr % lambda(nmat, nf))
        allocate(tr % neutron_velo(nmat, ng))

    end subroutine

    !===============================================================================================!
    ! set transient data
    !===============================================================================================!

    subroutine set_transient_pointer(tr, f, xsc, beta, lambda, neutron_velo, total_beta, dfis, ind_mat)

        class(trans_type)                        :: tr
        type(fdm_rect), intent(in), target       :: f                 ! FDM object
        type(xs_change_type), intent(in), target :: xsc
        real(dp), intent(in), target             :: beta(:,:)         ! beta (delayed neutron fraction)
        real(dp), intent(in), target             :: lambda(:,:)       ! lambda (precursor decay constant)
        real(dp), intent(in), target             :: neutron_velo(:,:) ! Neutron velocity
        real(dp), intent(in), target             :: total_beta(:)     ! total beta for each material
        real(dp), intent(in), target             :: dfis(:)         
        integer, intent(in), target              :: ind_mat(:)        ! material index

        tr % f => f
        tr % xs => f % xs
        tr % xsc => xsc
        tr % ind_mat => ind_mat
        tr % total_beta => total_beta
        tr % dfis => dfis

        tr % beta = beta
        tr % lambda = lambda
        tr % neutron_velo = neutron_velo

    end subroutine

    !===============================================================================================!
    ! To perform rod ejection simulation
    !===============================================================================================!

    subroutine rod_eject_transient(tr, is_th)
        
        class(trans_type)   :: tr
        logical, intent(in) :: is_th
        
        real(dp)               :: power_first
        real(dp)               :: rho
        real(dp)               :: t_now, t_next
        integer                :: i_mat, i, p, i_max, i_step
        real(dp)               :: ftemp_err
        real(dp), allocatable  :: tmp_flux(:,:)
        logical                :: converge
        
        ! Update xsec
        call xs_update(tr % f % xs, tr % xsc)
        
        ! Calculate forward flux at t=0 and check if keff=1
        write(output_unit,*)
        write(output_unit,'(A)', advance='no') '  steady state calculation starts ... '
        if (is_th) then
            call th_iteration(100, .false., ftemp_err, is_converge=converge)
        else
            call outer_iter(tr % f, print_iter = .false., is_converge=converge)
        end if
        if (.not. converge) call print_fail_converge()
        write(output_unit,*) 'done'
        
        ! If keff not equal to 1.0, force to 1.0
        if (abs(tr % f % Keff - 1._dp) > 1.e-5_dp) call adjust_keff(tr % f, tr % xsc)
        
        ! ! Calculate Adjoint flux
        ! write(output_unit,'(A)', advance='no') '  adjoint calculation starts ... '
        ! allocate(tmp_flux(tr % nnod, tr % ng))
        ! tmp_flux = tr % f % flux
        ! call outer_adjoint(fdm, print_iter = .false., is_nodal_updated=.false., is_converge=converge)
        ! if (.not. converge) call print_fail_converge()
        ! tr % adj_flux = tr % f % flux
        ! tr % f % flux = tmp_flux
        ! deallocate(tmp_flux)
        ! write(output_unit,*) 'done'
        ! print *, tr % f % keff
        tr % adj_flux = tr % f % flux
        
        
        ! Calculate Initial precursor density
        call precursor_init(tr)
        
        ! Calculate initial power
        call get_power(tr, power_first)
        
        ! Total beta
        tr % total_beta = 0.
        do p = 1, tr % nf
            do i_mat = 1, tr % nmat
                tr % total_beta(i_mat) = tr % total_beta(i_mat) + tr % beta(i_mat,p)
            end do
        end do

        if (.not. tr % is_tabular) then
            tr % core_beta = tr % total_beta(1)
        else
            tr % core_beta = get_core_beta(tr)

            write(fid,*)
            write(fid,1344) tr % core_beta*1.e5
            write(output_unit,*)
            write(output_unit,1344) tr % core_beta*1.e5
            1344 format ('  CORE AVERAGED DELAYED NEUTRON FRACTION: ', F7.2, ' PCM')
        end if
        
        ! Calculate reactivity
        call reactivity(tr, tr % f % xs % sigr, rho)
        
        ! File output
        call print_head()
        call print_transient_step(0, 0._dp, rho / tr % core_beta, 1._dp, converge=.true.)
        
        ! Start transient calculation
        i_step = 0
        t_next = 0.
        i_max = nint(tr % time_mid/ tr % t_step_1)
        
        ! First Time Step Width
        do i = 1, i_max
        
            i_step = i_step + 1
            t_now = t_next
            t_next = real(i) * tr % t_step_1
        
            if (i > 1) THEN
                if (tr % is_exp) then
                    tr % omega = log(tr % f % flux / tr % prev_flux) / tr % t_step_1
                else
                    tr % omega = 0._dp
                end if
            else
                tr % omega = 0._dp
            end if
        
            call trans_calc(tr, tr % t_step_1, power_first, i_step, t_next, is_th)
        
        end do
        
        ! Second Time Step Width
        i_max = nint((tr % total_time - tr % time_mid) / tr % t_step_2)
        
        do i = 1, i_max
        
            i_step = i_step + 1
            t_now = t_next
            t_next = tr % time_mid + real(i) * tr % t_step_2
        
            if (tr % is_exp) then
                tr % omega = log(tr % f % flux / tr % prev_flux) / tr % t_step_2
            else
                tr % omega = 0._dp
            end if
        
            call trans_calc(tr, tr % t_step_2, power_first, i_step, t_next, is_th)
        
        end do
        
    end subroutine

    !===============================================================================================!
    ! to perform transient calculations for given time step
    !===============================================================================================!

    subroutine trans_calc(tr, t_step, power_first, i_step, t_next, is_th)
        
        class(trans_type)                    :: tr
        real(dp), intent(in)                 :: t_step, power_first, t_next
        integer, intent(in)                  :: i_step
        logical, intent(in)                  :: is_th
        
        real(dp) :: rho
        real(dp) :: power_now
        integer  :: n, g
        logical  :: converge
        
        real(dp) :: lin_power(tr % nnod)          ! Linear power density
        real(dp) :: sigrp(tr % nnod, tr % ng)     ! sigr before modified
        real(dp) :: power_mult                    ! Power multiplication
        real(dp) :: tf, tm, mtf, mtm
        integer  :: i_mat
        
        ! Rod bank changes
        call move_crod(t_step, t_next)
        
        ! Calculate xsec after pertubation
        call xs_update(tr % f % xs, tr % xsc)
        
        ! Modify removal xsec
        sigrp = tr % xs % sigr
        do g = 1, tr % ng
            do n = 1, tr % nnod
                i_mat = tr % ind_mat(n)
                tr % xs % sigr(n,g) = sigrp(n,g) &
                + 1._dp / (tr % theta * tr % neutron_velo(i_mat,g) * t_step) &
                + tr % omega(n,g) / tr % neutron_velo(i_mat,g)
            end do
        end do
        
        ! Save the previous fluxes and fission source
        tr % prev_flux = tr % f % flux
        tr % prev_fsrc = tr % f % fsrc

        ! get transient external source
        tr % f % exsrc = get_exsrc(tr, sigrp, t_step)
        
        ! Transient calculation
        call transient_outer(tr % f, converge)
        
        ! Update precursor density
        call precursor_update(tr, t_step)
        
        ! Calculate power
        call get_power(tr, power_now)
        power_mult = power_now/power_first
        
        ! Calculate reactivity
        call reactivity(tr, sigrp, rho)
        
        if (is_th) then
            call th_solver_transient(power_mult, t_step)
        end if

        call print_transient_step(i_step, t_next, rho / tr % core_beta, power_mult, converge)
        
        
    end subroutine trans_calc

    !===============================================================================================!
    ! To calculate paramaters from previous time step combined as an external source
    ! This external source contain the terms that appear in the kinetic
    !    calculations but do not appear in static calculation
    !===============================================================================================!

    function get_exsrc(tr, sigrp, t_step) result(res_exsrc)
      
        class(trans_type)    :: tr
        real(dp), intent(in) :: sigrp(:,:)
        real(dp), intent(in) :: t_step
        real(dp)             :: res_exsrc(tr % nnod, tr % ng)
      
        real(dp) :: dt(tr % nnod), dtp(tr % nnod)
        integer  :: n, p, g, i_mat
        real(dp) :: a1, a2, pxe, pthet
      
        tr % dfis = 0._dp
        dt        = 0._dp
        dtp       = 0._dp
        do p = 1, tr % nf
            do n = 1, tr % nnod
                i_mat        = tr % ind_mat(n)
                pxe          = exp(-tr % lambda(i_mat,p)*t_step)
                a1           = (1._dp - pxe) / (tr % lambda(i_mat,p)*t_step)
                a2           = 1._dp - a1
                a1           = a1 - pxe
                tr % dfis(n) = tr % dfis(n) + tr % beta(i_mat,p) * a2
                dt(n)        = dt(n) + tr % lambda(i_mat,p) * tr % precursor(n,p) * pxe &
                + tr % beta(i_mat,p) * a1 * tr % prev_fsrc(n)
                dtp(n)       = dtp(n) + tr % lambda(i_mat,p) * tr % precursor(n,p)
            end do
        end do

        do g = 1, tr % ng
            do n = 1, tr % nnod
                ! to do: should we use old sigr?
                ! pthet => all terms from previous time step
                i_mat = tr % ind_mat(n)
                pthet = -tr % L(n,g) - sigrp(n,g) * tr % prev_flux(n,g) + tr % f % scat(n,g) &
                + (1._dp - tr % total_beta(i_mat)) * tr % xs % chi(n,g) * tr % prev_fsrc(n) &
                +  tr % xs % chi(n,g) * dtp(n)
                res_exsrc(n,g) = tr % xs % chi(n,g) * dt(n) &
                + exp(tr % omega(n,g) * t_step) * tr % prev_flux(n,g) / (tr % theta * tr % neutron_velo(i_mat, g) * t_step) &
                + tr % big_theta * pthet
            end do
        end do
      
    end function

    !===============================================================================================!
    ! Calculate power
    !===============================================================================================!

    subroutine get_power(tr,power)
        
        class(trans_type) tr
        real(dp), intent(out):: power
        
        real(dp) :: pn(tr % nnod)
        integer  :: g, n
        real(dp) :: pow
        
        pn = 0._dp
        do g = 1, tr % ng
            do n= 1, tr % nnod
              pow = tr % f % flux(n,g) * tr % xs % sigf(n,g) * tr % f % vdel(n)
              if (pow < 0.) pow = 0.
              pn(n) = pn(n) + pow
            end do
        end do
        
        power = 0._dp
        do n = 1, tr % nnod
            power = power + pn(n)
        end do
        
    end subroutine

    !===============================================================================================!
    ! Calculate Initial precursor density
    !===============================================================================================!

    subroutine precursor_init(tr)
        
        class(trans_type)    :: tr
        
        integer :: n, p, i_mat
        real(dp) :: blamb
        do p = 1, tr % nf
            do n = 1, tr % nnod
                i_mat = tr % ind_mat(n)
                blamb = tr % beta(i_mat,p) / tr % lambda(i_mat,p)
                tr % precursor(n,p)  = blamb * tr % f % fsrc(n)
            end do
        end do
        
        do n = 1, tr % nnod
            if (tr % xs % nuf(n,tr % ng) < 1.e-10) tr % precursor(n,:) = 0.0
        end do
        
    end subroutine

    !===============================================================================================!
    ! To update precursor density
    ! Using analytical solution of the delayed neutron precursor density differential equation
    !===============================================================================================!

    subroutine precursor_update(tr, ht)

        class(trans_type)    :: tr
        real(dp), intent(in) :: ht

        real(dp) :: a1, a2, pxe
        integer  :: n, p, i_mat
        
        do p = 1, tr % nf
            do n = 1, tr % nnod
                i_mat = tr % ind_mat(n)
                pxe  = exp(-tr % lambda(i_mat,p)*ht)
                a1   = (1._dp - pxe) / (tr % lambda(i_mat,p)*ht)
                a2   = 1._dp - a1
                a1   = a1 - pxe
                tr % precursor(n,p)  = tr % precursor(n,p)  * pxe &
                + tr % beta(i_mat,p) / tr % lambda(i_mat,p) * (a1*tr % prev_fsrc(n) + a2*tr % f % fsrc(n))
            end do
        end do
        
    end subroutine

    !===============================================================================================!
    ! To calculate dynamic reactivity
    !===============================================================================================!

    subroutine reactivity(tr, sigrp, rho)
        
        class(trans_type)      :: tr
        real(dp), intent(in)   :: sigrp(:,:)
        real(dp), intent(out)  :: rho
        
        integer  :: n, g, h
        real(dp) :: rem, lea, src, fde
        real(dp) :: Lx                  ! Leakages in x direction
        real(dp) :: Ly                  ! Leakages in y direction
        real(dp) :: Lz                  ! Leakages in z direction

        src = 0.; rem = 0.; lea = 0.; fde = 0.
        do g = 1, tr % ng
            do n = 1, tr % nnod
                call get_leakage(tr % f, n, g, Lx, Ly, Lz)
                tr % L(n,g) = Lx + Ly + Lz
                src = src + tr % adj_flux(n,g) * (tr % f % scat(n,g) + tr % xs % chi(n,g) * tr % f % fsrc(n)) * tr % f % vdel(n)  ! total source
                rem = rem + tr % adj_flux(n,g) * sigrp(n,g) * tr % f % flux(n,g) * tr % f % vdel(n)                      ! removal
                lea = lea + tr % adj_flux(n,g) * tr % L(n,g) * tr % f % vdel(n)                                                   ! leakages
                fde = fde + tr % adj_flux(n,g) * tr % xs % chi(n,g) * tr % f % fsrc(n) * tr % f % vdel(n)                         ! fission source
            end do
        end do

        rho = (src - lea - rem) / fde
        
    end subroutine

    !===============================================================================================!
    ! To calculate core-averaged delayed neutron fraction
    ! used only if tabular xsec is used
    !===============================================================================================!

    pure function get_core_beta(tr) result(core_avg_beta)
        
        class(trans_type), intent(in)    :: tr
        real(dp)                         :: core_avg_beta
        
        INTEGER :: n, p, g, i_mat
        real(dp), dimension(tr % nnod) :: vdum, vdum2
        real(dp) :: F
        
        ! Calculate F
        vdum = 0.
        do g = 1, tr % ng
            do n = 1, tr % nnod
                vdum(n) = vdum(n) + tr % xs % nuf(n,g) * tr % f % flux(n,g)
            end do
        end do
        
        vdum2 = 0.
        do g = 1, tr % ng
            do n = 1, tr % nnod
                vdum2(n) = vdum2(n) + tr % xs % chi(n,g) * vdum(n) * tr % adj_flux(n,g)
            end do
        end do
        
        F = integrate(vdum2)
        
        ! Calculate Delayed neutron fraction (beta)
        core_avg_beta = 0.
        do p = 1, tr % nf
            vdum2 = 0.
            do g = 1, tr % ng
                do n = 1, tr % nnod
                    vdum2(n) = vdum2(n) + tr % xs % chi(n,g) * tr % beta(i_mat,p) * vdum(n) * tr % adj_flux(n,g)
                end do
            end do
            core_avg_beta = core_avg_beta + integrate(vdum2) / F
        end do
        
    end function


end module