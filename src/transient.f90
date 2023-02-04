module transient

    use iso_fortran_env, only: real64, error_unit, output_unit
    use xsec, only: xs_rect
    use transient_aux
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

    real(dp), pointer          :: beta(:)                  ! beta (delayed neutron fraction)
    real(dp), pointer          :: lambda(:)                ! lambda (precursor decay constant)
    real(dp)                   :: core_beta                ! Core averaged beta
    real(dp), pointer          :: neutron_velo(:)          ! Neutron velocity
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
    real(dp), allocatable      :: mod_sigr(:,:)         ! Modified sigma removal to take account of terms that appear in transient calc.
    logical                    :: is_tabular            ! is tabular cross section used
    logical                    :: is_exp                ! is exponential transformation is used

    type(fdm_rect), pointer       :: f                     ! FDM object
    type(xs_rect), pointer        :: xs                    ! node-wise xs
    type(xs_change_type), pointer :: xsc                   ! xs changes
    integer, pointer              :: ind_mat(:)            ! material index
    end type

    public :: set_transient_data, set_transient_pointer
    public :: rod_eject_no_th

    contains

    !===============================================================================================!
    ! set transient data
    !===============================================================================================!

    subroutine set_transient_data(tr, fid_inp, nf, nnod, ng, nmat, total_time, time_step_1, time_step_2, &
        time_mid, theta, bxtab, bextr)

        class(trans_type)                        :: tr
        integer, intent(in)                      :: fid_inp         ! number of energy group
        integer, intent(in)                      :: nf              ! number delayed neutron precursor group family
        integer, intent(in)                      :: ng              ! number of energy group
        integer, intent(in)                      :: nnod            ! Total number of nodes
        integer, intent(in)                      :: nmat            ! Number of materials
        real(dp), intent(in)                     :: total_time      ! TOTAL SIMULATION TIME
        real(dp), intent(in)                     :: time_step_1     ! FIRST TIME STEP
        real(dp), intent(in)                     :: time_step_2     ! SECOND TIME STEP
        real(dp), intent(in)                     :: time_mid        ! WHEN SECOND TIME STEP START
        real(dp), intent(in)                     :: theta           ! Small theta
        integer, intent(in)                      :: bxtab           ! is tabular xs used
        integer, intent(in)                      :: bextr           ! is exponential transformation is used

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
        if (bxtab == 1) then
            tr % is_tabular = .true.
        else
            tr % is_tabular = .false.
        end if
        if (bextr == 1) then
            tr % is_exp = .true.
        else
            tr % is_exp = .false.
        end if

        allocate(tr % prev_flux(nnod,ng))
        allocate(tr % prev_fsrc(nnod))
        allocate(tr % omega(nnod,ng))
        allocate(tr % neutron_velo(ng))
        allocate(tr % precursor(nnod,nf))
        allocate(tr % L(nnod,ng))
        allocate(tr % adj_flux(nnod,ng))
        allocate(tr % mod_sigr(nnod,ng))

    end subroutine

    !===============================================================================================!
    ! set transient data
    !===============================================================================================!

    subroutine set_transient_pointer(tr, f, xsc, beta, lambda, neutron_velo, total_beta, dfis, ind_mat)

        class(trans_type)                        :: tr
        type(fdm_rect), intent(in), target      :: f               ! FDM object
        type(xs_change_type), intent(in), target :: xsc
        real(dp), intent(in), target             :: beta(:)         ! beta (delayed neutron fraction)
        real(dp), intent(in), target             :: lambda(:)       ! lambda (precursor decay constant)
        real(dp), intent(in), target             :: neutron_velo(:) ! Neutron velocity
        real(dp), intent(in), target             :: total_beta(:)   ! total beta for each material
        real(dp), intent(in), target             :: dfis(:)         
        integer, intent(in), target              :: ind_mat(:)      ! material index

        tr % f => f
        tr % xs => f % xs
        tr % xsc => xsc
        tr % ind_mat => ind_mat
        tr % beta => beta
        tr % lambda => lambda
        tr % neutron_velo => neutron_velo
        tr % total_beta => total_beta
        tr % dfis => dfis

    end subroutine

    !===============================================================================================!
    ! To perform rod ejection simulation without TH feedback
    !===============================================================================================!

    subroutine rod_eject_no_th(tr)
        
        class(trans_type) :: tr
        
        real(dp) :: power_first
        real(dp) :: rho
        real(dp) :: t_now, t_next
        integer  :: n, i, j, i_max, i_step
        logical  :: converge
        
        ! Update xsec
        call xs_update(tr % f % xs, tr % xsc)
        
        ! Calculate forward flux at t=0 and check if keff=1
        write(output_unit,*)
        write(output_unit,'(A)', advance='no') '  steady state calculation starts ... '
        call outer_iter(tr % f, print_iter = .false., is_converge=converge)
        if (.not. converge) call print_fail_converge()
        write(output_unit,*) 'done'
        
        ! If keff not equal to 1.0, force to 1.0
        if (abs(tr % f % Keff - 1._dp) > 1.e-5_dp) call adjust_keff(tr % f, tr % xsc)
        
        ! Calculate Adjoint flux
        write(output_unit,'(A)', advance='no') '  adjoint calculation starts ... '
        call get_adjoint_flux(tr % xsc, tr % adj_flux)
        write(output_unit,*) 'done'
        
        ! Calculate Initial precursor density
        call precursor_init(tr)
        
        ! Calculate initial power
        call get_power(tr, power_first)
        
        ! Total beta
        tr % total_beta = 0.
        do j = 1, tr % nf
          do n = 1, tr % nmat
            tr % total_beta(n) = tr % total_beta(n) + tr % beta(j)
          end do
        end do
        tr % core_beta = tr % total_beta(1)
        
        ! Calculate reactivity
        call reactivity(tr, rho)
        
        ! File output
        call print_head()
        call print_transient_step(0, 0._dp, rho / tr % core_beta, 1._dp)
        
        ! Start transient calculation
        i_step = 0
        t_next = 0.
        i_max = nint(tr % time_mid/ tr % t_step_1)
        
        ! First Time Step
        write(*,*) sum(tr % xs % sigr)
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
        
            call trans_calc(tr, tr % t_step_1, power_first, i_step, t_next, is_th=.false.)
        
        end do
        
        ! Second Time Step
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
        
            call trans_calc(tr, tr % t_step_2, power_first, i_step, t_next, is_th=.false.)
        
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
        
        real(dp), dimension(tr % nnod) :: lin_power       ! Linear power density
        real(dp) :: power_mult                            ! Power multiplication
        real(dp) :: tf, tm, mtf, mtm
        
        ! Rod bank changes
        call move_crod(t_step, t_next)
        
        ! Calculate xsec after pertubation
        call xs_update(tr % f % xs, tr % xsc)
        write(*,*) sum(tr % xs % sigr)
        
        ! Modify removal xsec
        if (.not. tr % is_tabular) then
            do g = 1, tr % ng
                do n = 1, tr % nnod
                    tr % mod_sigr(n,g) = tr % xs % sigr(n,g) &
                    + 1._dp / (tr % theta * tr % neutron_velo(g) * t_step) &
                    + tr % omega(n,g) / tr % neutron_velo(g)
                end do
            end do
        ! ELSE
        !   DO g = 1, ng
        !     DO n = 1, nnod
        !       sigr(n,g) = sigr(n,g) + 1._DP / (sth * m(mat(n))%velo(g) * t_step) &
        !                 + omeg(n,g) / m(mat(n))%velo(g)
        !     END DO
        !   END DO
        end IF
        
        ! Save the previous fluxes and fission source
        tr % prev_flux = tr % f % flux
        tr % prev_fsrc = tr % f % fsrc

        ! get transient external source
        tr % f % exsrc = get_exsrc(tr, t_step)
        
        ! Transient calculation
        call transient_outer(tr % f, tr % mod_sigr, converge)
        if (.not. converge) call print_fail_converge()
        
        ! Update precursor density
        call precursor_update(tr, t_step)
        
        ! Calculate power
        call get_power(tr, power_now)
        
        ! Calculate reactivity
        call reactivity(tr, rho)
        
        ! if (thc == 1) then
        !   ! Calculate node power distribution
        !   call PowDis(npow)
        
        !   ! Power change
        !   power_mult = ppow * tpow2/power_first * 0.01_DP
        
        !   ! Calculate linear power density for each nodes (W/cm)
        !   DO n = 1, nnod
        !      lin_power(n) = npow(n) * pow * power_mult &
        !      / (node_nf(ix(n),iy(n)) * zdel(iz(n)))     ! Linear power density (W/cm)
        !   END DO
        
        !   ! TH transient
        !   call th_trans(lin_power,t_step)
        
        !   call par_ave(ftem, tf)
        !   call par_max(tfm(:,1), mtf)
        !   call par_ave(mtem, tm)
        !   call par_max(mtem, mtm)
        ! end if

        call print_transient_step(i_step, t_next, rho / tr % core_beta, power_now/power_first)
        
        
    end subroutine trans_calc

    !===============================================================================================!
    ! To calculate paramaters from previous time step combined as an external source
    ! This external source contain the terms that appear in the kinetic
    !    calculations but do not appear in static calculation
    !===============================================================================================!

    function get_exsrc(tr, t_step) result(exsrc)
      
        class(trans_type)    :: tr
        real(dp), intent(in) :: t_step
        real(dp)             :: exsrc(tr % nnod, tr % ng)
      
        real(dp) :: dt, dtp
        integer  :: n, i, g
        real(dp) :: a1, a2, pxe, pthet
      
        if (.not. tr % is_tabular) then
            tr % dfis = 0._dp
            do n = 1, tr % nnod
                dt = 0._dp; dtp = 0._dp
                do i = 1, tr % nf
                    pxe          = exp(-tr % lambda(i)*t_step)
                    a1           = (1._dp - pxe) / (tr % lambda(i)*t_step)
                    a2           = 1._dp - a1
                    a1           = a1 - pxe
                    tr % dfis(n) = tr % dfis(n) + tr % beta(i) * a2
                    dt           = dt   + tr % lambda(i) * tr % precursor(n,i) * pxe &
                    + tr % beta(i) * a1 * tr % prev_fsrc(n)
                    dtp          = dtp + tr % lambda(i) * tr % precursor(n,i)
                end do
  
                do g = 1, tr % ng
                    ! to do: should we use new sigr?
                    ! pther => all terms from previous time step
                    pthet = -tr % L(n,g) - tr % xs % sigr(n,g) * tr % prev_flux(n,g) + tr % f % scat(n,g) &
                    + (1._dp - tr % total_beta(tr % ind_mat(n))) * tr % xs % chi(n,g) * tr % prev_fsrc(n) &
                    +  tr % xs % chi(n,g) * dtp
                    exsrc(n,g) = tr % xs % chi(n,g) * dt &
                    + exp(tr % omega(n,g) * t_step) * tr % prev_flux(n,g) / (tr % theta * tr % neutron_velo(g) * t_step) &
                    + tr % big_theta * pthet
                end do
            end do
        ! else
        !     dfis = 0._dp
        !     do n = 1, nnod
        !       dt = 0._dp; dtp = 0._dp
        !       do i = 1, tr % nf
        !         pxe  = exp(-m(mat(n))%lamb(i)*t_step)
        !         if (nuf(n,ng) > 0.) then
        !           a1   = (1._dp - pxe) / (m(mat(n))%lamb(i)*t_step)
        !         else
        !           a1 = 0._dp
        !         end if
        !         a2   = 1._dp - a1
        !         a1   = a1 - pxe
        !         dfis(n) = dfis(n) + m(mat(n))%iBeta(i) * a2
        !         dt   = dt   + m(mat(n))%lamb(i) * c0(n,i) * pxe &
        !              + m(mat(n))%iBeta(i) * a1 * fst(n)
        !         dtp  = dtp + m(mat(n))%lamb(i) * c0(n,i)
        !       end do
        
        !       do g = 1, ng
        !         pthet = -L(n,g) - sigrp(n,g) * ft(n,g) + s0(n,g) &
        !               + (1._dp - tbeta(mat(n))) * chi(mat(n),g) * fst(n) &
        !               +  chi(mat(n),g) * dtp
        !         exsrc(n,g) = chi(mat(n),g) * dt &
        !                    + exp(omeg(n,g) * t_step) * ft(n,g) &
        !                    / (sth * m(mat(n))%velo(g) * t_step) + bth * pthet
        !       end do
        !     end do
        end if
      
    end function

    !===============================================================================================!
    ! Calculate power
    !===============================================================================================!

    subroutine get_power(tr,power)
        
        class(trans_type) tr
        real(dp), intent(out):: power
        
        real(dp) :: p(tr % nnod)
        integer  :: g, n
        real(dp) :: pow
        
        p = 0._dp
        do g = 1, tr % ng
            do n= 1, tr % nnod
              pow = tr % f % flux(n,g) * tr % xs % sigf(n,g) * tr % f % vdel(n)
              if (pow < 0.) pow = 0.
              p(n) = p(n) + pow
            end do
        end do
        
        power = 0._dp
        do n = 1, tr % nnod
            power = power + p(n)
        end do
        
    end subroutine

    !===============================================================================================!
    ! Calculate Initial precursor density
    !===============================================================================================!

    subroutine precursor_init(tr)
        
        class(trans_type)    :: tr
        
        integer :: n, i
        real(dp) :: blamb
        
        ! if (bxtab == 1) THEN
        !   do n = 1, nnod
        !      do j = 1, tr % nf
        !        if (nuf(n,ng) > 0.) THEN  !If it is fuel
        !          blamb = m(mat(n))%iBeta(j) / m(mat(n))%lamb(j)
        !          c0(n,j)  = blamb * fs0(n)
        !        else
        !          c0(n,j) = 0.
        !        end if
        !      end do
        !   end do
        ! else
          do i = 1, tr % nf
             do n = 1, tr % nnod
                blamb = tr % beta(i) / tr % lambda(i)
                tr % precursor(n,i)  = blamb * tr % f % fsrc(n)
             end do
          end do
        ! end if
        
        
    end subroutine

    !===============================================================================================!
    ! To update precursor density
    ! Using analytical solution of the delayed neutron precursor density differential equation
    !===============================================================================================!

    subroutine precursor_update(tr, ht)

        class(trans_type)    :: tr
        real(dp), intent(in) :: ht

        real(dp) :: a1, a2, pxe
        integer  :: n, i
        
        ! if (bxtab == 1) THEN
        !   do i = 1, tr % nf
        !     do n = 1, nnod
        !       if (nuf(n,ng) > 0.) THEN  !If it is fuel
        !         pxe  = exp(-m(mat(n))%lamb(i)*ht)
        !         a1   = (1._dp - pxe) / (m(mat(n))%lamb(i)*ht)
        !         a2   = 1._dp - a1
        !         a1   = a1 - pxe
        !         c0(n,i)  = c0(n,i)  * pxe + m(mat(n))%iBeta(i) / m(mat(n))%lamb(i) &
        !                  * (a1*fst(n) + a2*fs0(n))
        !       end if
        !     end do
        !   end do
        ! else
          do i = 1, tr % nf
            pxe  = exp(-tr % lambda(i)*ht)
            a1   = (1._dp - pxe) / (tr % lambda(i)*ht)
            a2   = 1._dp - a1
            a1   = a1 - pxe
            do n = 1, tr % nnod
                tr % precursor(n,i)  = tr % precursor(n,i)  * pxe &
                + tr % beta(i) / tr % lambda(i) * (a1*tr % prev_fsrc(n) + a2*tr % f % fsrc(n))
            end do
          end do
        ! end if
        
        
    end subroutine

    !===============================================================================================!
    ! To calculate dynamic reactivity
    !===============================================================================================!

    subroutine reactivity(tr, rho)
        
        class(trans_type)      :: tr
        real(dp), intent(out)  :: rho
        
        integer  :: n, g, h
        real(dp) :: scg(tr % nnod)   ! scattering source
        real(dp) :: rem, lea, src, fde
        real(dp) :: Lx(tr % nnod, tr % ng)      ! Leakages in x direction
        real(dp) :: Ly(tr % nnod, tr % ng)      ! Leakages in y direction
        real(dp) :: Lz(tr % nnod, tr % ng)      ! Leakages in z direction

        call get_leakage(tr % f, Lx, Ly, Lz)
        
        src = 0.; rem = 0.; lea = 0.; fde = 0.
        do g = 1, tr % ng
            scg = 0.
            do h = 1, tr % ng
                do n = 1, tr % nnod
                    if (g .ne. h) scg(n) = scg(n) + tr % xs % sigs(n,h,g) * tr % f % flux(n,h)
                end do
            end do
            do n = 1, tr % nnod
                tr % L(n,g) = Lx(n,g) + Ly(n,g) + Lz(n,g)
                src = src + tr % adj_flux(n,g) * (scg(n) + tr % xs % chi(n,g) * tr % f % fsrc(n)) * tr % f % vdel(n)  ! total source
                rem = rem + tr % adj_flux(n,g) * tr % xs % sigr(n,g) * tr % f % flux(n,g) * tr % f % vdel(n)                     ! removal
                lea = lea + tr % adj_flux(n,g) * (Lx(n,g) + Ly(n,g) + Lz(n,g)) * tr % f % vdel(n)                  ! leakages
                fde = fde + tr % adj_flux(n,g) * tr % xs % chi(n,g) * tr % f % fsrc(n) * tr % f % vdel(n)                             ! fission source
            end do
        end do

        rho = (src - lea - rem) / fde
        
    end subroutine

    !===============================================================================================!
    ! To calculate core-averaged delayed neutron fraction
    ! used only if tabular xsec is used
    !===============================================================================================!

    ! subroutine get_core_beta(tr)
        
    !     class(trans_type)    :: tr
        
    !     INTEGER :: n, i, g
    !     real(dp), dimension(tr % nnod) :: vdum, vdum2
    !     real(dp) :: F
        
    !     ! Calculate F
    !     vdum = 0.
    !     DO g = 1, tr % ng
    !         DO n = 1, tr % nnod
    !             vdum(n) = vdum(n) + tr % xs % nuf(n,g) * tr % flux(n,g)
    !         END DO
    !     END DO
        
    !     vdum2 = 0.
    !     DO g = 1, tr % ng
    !         DO n = 1, tr % nnod
    !             vdum2(n) = vdum2(n) + tr % xs % chi(n,g) * vdum(n) * tr % adj_flux(n,g)
    !         END DO
    !     END DO
        
    !     F = integrate(vdum2)
        
    !     ! Calculate Delayed neutron fraction (beta)
    !     tr % core_beta = 0._dp
    !     DO i = 1, tr % nf
    !         vdum2 = 0.
    !         DO g = 1, tr % ng
    !             DO n = 1, tr % nnod
    !                 vdum2(n) = vdum2(n) + chi(mat(n),g) * m(mat(n))%iBeta(i) * vdum(n) * af(n,g)
    !             END DO
    !         END DO
    !         tr % core_beta =tr % core_beta + integrate(vdum2) / F
    !     END DO
        
    !     write(fid,*)
    !     write(fid,1344) tr % core_beta*1.e5
        
    !     1344 format ('  CORE AVERAGED DELAYED NEUTRON FRACTION: ', F7.2, ' PCM')
        
        
    ! end subroutine


end module