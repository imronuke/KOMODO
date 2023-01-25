module transient

    use iso_fortran_env, only: real64, output_unit
    use fdm
    use nodal
    use xsec

    implicit none

    private

    save

    integer, parameter :: dp = real64
    integer, parameter :: nf = 6                ! Number of delayed neutron precursor family
    integer            :: fid

    type, public :: trans_type
    integer                    :: nnod                  ! Total number of nodes
    integer                    :: ng                    ! number of energy group
    real(dp)                   :: theta = 1._dp         ! Small theta (default is fully implicit)
    real(dp)                   :: big_theta = 0._dp     ! big theta for transient using theta method
    real(dp)                   :: omega                 ! Exponential transformation constant

    real(dp)                   :: beta(nf)              ! beta (delayed neutron fraction)
    real(dp)                   :: lambda(nf)            ! lambda (precursor decay constant)
    real(dp)                   :: core_beta             ! Core averaged beta
    real(dp), allocatable      :: neutron_velo(:)       ! Neutron velocity
    real(dp), allocatable      :: precursor(:,:)         ! precursor density

    real(dp)                   :: total_time            ! TOTAL SIMULATION TIME
    real(dp)                   :: time_step_1           ! FIRST TIME STEP
    real(dp)                   :: time_step_2           ! SECOND TIME STEP
    real(dp)                   :: time_mid              ! WHEN SECOND TIME STEP START

    real(dp), allocatable      :: L(:,:)                ! Total leakages for node n and group g
    real(dp), allocatable      :: adj_flux(:,:)         ! Adjoint flux

    real(dp), pointer          :: flux(:,:)             ! nodes flux
    real(dp), pointer          :: fsrc(:)               ! nodes fission source
    real(dp), pointer          :: prev_flux(:,:)        ! previous time-step nodes flux
    real(dp), pointer          :: prev_fsrc(:)          ! previous time-step nodes fission source
    type(core_rect), pointer   :: f                     ! FDM object
    type(nodal_type), pointer  :: d                     ! Nodal object
    type(xs_rect), pointer     :: xs                    ! node-wise xs

    end type

    contains

    !===============================================================================================!
    ! set transient data
    !===============================================================================================!

    subroutine set_transient(t, fid_inp, f, flux, fsrc, prev_flux, prev_fsrc, ng, nnod, &
        total_time, time_step_1, time_step_2, time_mid, theta)

        class(trans_type)                    :: t
        integer, intent(in)                  :: fid_inp         ! number of energy group
        type(core_rect), intent(in), target  :: f               ! FDM object
        real(dp), intent(in), target         :: flux(:,:)       ! nodes flux
        real(dp), intent(in), target         :: fsrc(:)         ! nodes fission source
        real(dp), intent(in), target         :: prev_flux(:,:)  ! previous time-step nodes flux
        real(dp), intent(in), target         :: prev_fsrc(:)    ! previous time-step nodes fission source
        integer, intent(in)                  :: ng              ! number of energy group
        integer, intent(in)                  :: nnod            ! Total number of nodes
        real(dp), intent(in)                 :: total_time      ! TOTAL SIMULATION TIME
        real(dp), intent(in)                 :: time_step_1     ! FIRST TIME STEP
        real(dp), intent(in)                 :: time_step_2     ! SECOND TIME STEP
        real(dp), intent(in)                 :: time_mid        ! WHEN SECOND TIME STEP START
        real(dp), intent(in), optional       :: theta           ! Small theta

        fid      = fid_inp
        t % ng   = ng
        t % nnod = nnod
        t % total_time = total_time
        t % time_step_1 = time_step_1
        t % time_step_2 = time_step_2
        t % time_mid    = time_mid
        if (present(theta)) then
            t % theta = theta
            t % big_theta = 1. - 1. / theta
        end if

        t % d => f % d
        t % f => f
        t % xs => f % xs

        t % flux => flux
        t % fsrc => fsrc
        t % prev_flux => prev_flux
        t % prev_fsrc => prev_fsrc

        allocate(t % neutron_velo(ng))
        allocate(t % precursor(nf,ng))
        allocate(t % L(nnod,ng))
        allocate(t % adj_flux(nnod,ng))

    end subroutine

    !===============================================================================================!
    ! To adjuts the Keff to 1.0 if it is not equal to 1.0
    !===============================================================================================!

    ! subroutine adjust_keff()
        
    !     integer :: i
        
    !     WRITE(fid, *)
    !     WRITE(fid, '(A46,F9.6)') '  INITIAL MULTIPLICATION EFFECTIVE (K-EFF) = ', Ke
    !     WRITE(fid, *) '  WARNING: THE STEADY STATE K-EFF IS NOT EQUAL TO 1.0'
    !     WRITE(fid, *) '  AND NOW IT IS FORCED TO 1.0 BY MODIFYING THE nu*sigf CROSS SECTIONS '
    !     WRITE(fid, *)
    !     write(output_unit,*)
    !     write(output_unit,*) ' steady state k-eff not equal to one, force it to one ... done'
    !     do i = 1, 10
    !        xnuf = xnuf / t % d % Keff
    !        dnuf = dnuf / t % d % Keff
    !        CALL XS_updt(bcon, ftem, mtem, cden, bpos)
    !        CALL outer(0)
    !        if (ABS(Ke-1._DP) < 1.e-5_DP) EXIT
    !     end do
    !     if (i == 10) STOP "K-EFF STILL NOT EQUAL TO ONE. KOMODO IS STOPPING"
        
        
    ! end subroutine

    !===============================================================================================!
    ! Calculate Initial precursor density
    !===============================================================================================!

    subroutine precursor_init(t)
        
        class(trans_type)    :: t
        
        integer :: n, i
        real(dp) :: blamb
        
        ! if (bxtab == 1) THEN
        !   do n = 1, nnod
        !      do j = 1, nf
        !        if (nuf(n,ng) > 0.) THEN  !If it is fuel
        !          blamb = m(mat(n))%iBeta(j) / m(mat(n))%lamb(j)
        !          c0(n,j)  = blamb * fs0(n)
        !        else
        !          c0(n,j) = 0.
        !        end if
        !      end do
        !   end do
        ! else
          do i = 1, nf
             do n = 1, t % nnod
                blamb = t % beta(i) / t % lambda(i)
                t % precursor(n,i)  = blamb * t % fsrc(n)
             end do
          end do
        ! end if
        
        
    end subroutine

    !===============================================================================================!
    ! To update precursor density
    !===============================================================================================!

    subroutine precursor_update(t, ht)

        class(trans_type)    :: t
        real(dp), intent(in) :: ht

        real(dp) :: a1, a2, pxe
        integer  :: n, i
        
        ! if (bxtab == 1) THEN
        !   do i = 1, nf
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
          do i = 1, nf
            pxe  = exp(-t % lambda(i)*ht)
            a1   = (1._dp - pxe) / (t % lambda(i)*ht)
            a2   = 1._dp - a1
            a1   = a1 - pxe
            do n = 1, t % nnod
                t % precursor(n,i)  = t % precursor(n,i)  * pxe &
                + t % beta(i) / t % lambda(i) * (a1*t % prev_fsrc(n) + a2*t % fsrc(n))
            end do
          end do
        ! end if
        
        
    end subroutine

    !===============================================================================================!
    ! To calculate dynamic reactivity
    !===============================================================================================!

    subroutine reactivity(t, sigrp, rho)
        
        class(trans_type)                    :: t
        real(dp), dimension(:,:), intent(in) :: sigrp
        real(dp), intent(out)                 :: rho
        
        integer  :: n, g, h
        real(dp) :: scg(t % nnod)   ! scattering source
        real(dp) :: rem, lea, src, fde
        real(dp) :: L1(t % nnod, t % ng)      ! Leakages in x direction
        real(dp) :: L2(t % nnod, t % ng)      ! Leakages in y direction
        real(dp) :: L3(t % nnod, t % ng)      ! Leakages in z direction

        call TLUpd0(t % d, L1, L2, L3)
        
        src = 0.; rem = 0.; lea = 0.; fde = 0.
        do g = 1, t % ng
          scg = 0.
          do h = 1, t % ng
             do n = 1, t % nnod
                if (g .ne. h) scg(n) = scg(n) + t % xs % sigs(n,h,g) * t % flux(n,h)
             end do
          end do
          do n = 1, t % nnod
            src = src + t % adj_flux(n,g) * (scg(n) + t % xs % chi(n,g) * t % fsrc(n)) * t % f % vdel(n)  ! total source
            rem = rem + t % adj_flux(n,g) * sigrp(n,g) * t % flux(n,g) * t % f % vdel(n)                     ! removal
            lea = lea + t % adj_flux(n,g) * (L1(n,g) + L2(n,g) + L3(n,g)) * t % f % vdel(n)                  ! leakages
            fde = fde + t % adj_flux(n,g) * t % xs % chi(n,g) * t % fsrc(n) * t % f % vdel(n)                             ! fission source
           end do
        end do
        
        rho = (src - lea - rem) / fde
        
    end subroutine reactivity

    !===============================================================================================!
    ! To calculate core-averaged delayed neutron fraction
    !===============================================================================================!

    ! subroutine get_core_beta(t)
        
    !     class(trans_type)    :: t
        
    !     INTEGER :: n, i, g
    !     real(dp), dimension(t % nnod) :: vdum, vdum2
    !     real(dp) :: F
        
    !     ! Calculate F
    !     vdum = 0.
    !     DO g = 1, t % ng
    !         DO n = 1, t % nnod
    !             vdum(n) = vdum(n) + t % nuf(n,g) * t % flux(n,g)
    !         END DO
    !     END DO
        
    !     vdum2 = 0.
    !     DO g = 1, t % ng
    !         DO n = 1, t % nnod
    !             vdum2(n) = vdum2(n) + t % chi(n,g) * vdum(n) * t % adj_flux(n,g)
    !         END DO
    !     END DO
        
    !     F = integrate(t, vdum2)
        
    !     ! Calculate Delayed neutron fraction (beta)
    !     t % core_beta = 0._dp
    !     DO i = 1, nf
    !         vdum2 = 0.
    !         DO g = 1, t % ng
    !             DO n = 1, t % nnod
    !                 vdum2(n) = vdum2(n) + chi(mat(n),g) * m(mat(n))%iBeta(i) * vdum(n) * af(n,g)
    !             END DO
    !         END DO
    !         t % core_beta =t % core_beta + integrate(t, vdum2) / F
    !     END DO
        
    !     write(fid,*)
    !     write(fid,1344) t % core_beta*1.e5
        
    !     1344 format ('  CORE AVERAGED DELAYED NEUTRON FRACTION: ', F7.2, ' PCM')
        
        
    ! end subroutine

    !===============================================================================================!
    ! do volume integration on variable a
    !===============================================================================================!

    ! pure function integrate(t, a) result(res)

    !     class(trans_type), intent(in) :: t
    !     real(dp), intent(in)          :: a(:)
    !     real(dp)                      :: res

    !     integer                 :: n

    !     res = 0.
    !     do n = 1, t % nnod
    !         res = res + a(n) * t % vdel(n)
    !     end do
    
    ! end function

end module