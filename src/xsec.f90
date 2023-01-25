module xsec

    use iso_fortran_env, only: real64
    use time

    implicit none

    private

    save

    integer, parameter :: dp = real64
    
    ! Node-wise xs base
    type :: xs_base
        real(dp), allocatable :: D(:,:)              ! Diffusion coefficient
        real(dp), allocatable :: sigtr(:,:)          ! Transport macroscopic XSEC
        real(dp), allocatable :: siga (:,:)          ! Absorption macroscopic XSEC
        real(dp), allocatable :: nuf  (:,:)          ! nu*fission macroscopic XSEC
        real(dp), allocatable :: sigf (:,:)          ! fission macroscopic XSEC
        real(dp), allocatable :: sigr (:,:)          ! Removal macroscopic XSEC
        real(dp), allocatable :: sigs (:,:,:)        ! Scattering macroscopic XSEC
        real(dp), allocatable :: chi(:,:)            ! MG fission spectrum
    end type

    type, extends(xs_base), public :: xs_rect
        real(dp), pointer :: dc(:,:,:)
    end type
    
    ! Change of XS per change of other parameters
    type :: xs_changes
        real(dp)              :: ref                 ! Reference parameter
        real(dp), allocatable :: sigtr(:,:)          ! Transport macroscopic XSEC
        real(dp), allocatable :: siga (:,:)          ! Absorption macroscopic XSEC
        real(dp), allocatable :: nuf  (:,:)          ! nu*fission macroscopic XSEC
        real(dp), allocatable :: sigf (:,:)          ! fission macroscopic XSEC
        real(dp), allocatable :: sigs (:,:,:)        ! Scattering macroscopic XSEC
    end type

    type, extends(xs_changes) :: xs_changes_rect  ! if xtab file is used
        real(dp), allocatable :: dc (:,:)
    end type

    type, public :: xs_change_type
        type(xs_changes), allocatable  :: ftem
        type(xs_changes), allocatable  :: mtem
        type(xs_changes), allocatable  :: cden
        type(xs_changes), allocatable  :: bcon
        type(xs_changes), allocatable  :: crod
    end type

    type :: xs_mat  ! Material wise XS [dimension(nmat, ng)]
        real(dp), pointer :: sigtr(:,:)          ! Transport macroscopic XSEC
        real(dp), pointer :: siga (:,:)          ! Absorption macroscopic XSEC
        real(dp), pointer :: nuf  (:,:)          ! nu*fission macroscopic XSEC
        real(dp), pointer :: sigf (:,:)          ! fission macroscopic XSEC
        real(dp), pointer :: sigs (:,:,:)        ! Scattering macroscopic XSEC
        real(dp), pointer :: chi  (:,:)          ! neutron fission spectrum
        real(dp), pointer :: dc   (:,:,:)        ! ADF
    end type

    type :: xs_auxilary    ! XS auxilary data used in this module
        integer, pointer      :: ind_mat(:)             ! node-wise material map
        integer               :: nnod, ng               ! number of nodes and group
        ! Below used for Control rod 
        integer, pointer      :: fbmap(:,:)             ! Radial control rod bank map (node wise)
        integer, pointer      :: xyz(:,:,:)             ! node index for given index in x, y, and z directions
        real(dp), pointer     :: zdel(:)                ! delta z
        integer               :: nxx, nyy, nzz
        real(dp)              :: coreh                  ! Core Height
        real(dp)              :: zero_pos, step_size    ! Zero step position and step size
    
    end type
    
    type(xs_mat)      :: xsm   ! Material wise XS
    type(xs_auxilary) :: xsa   ! XS auxilary data used in this module


    public :: set_xs_data, xsec_setup, set_xs_change, xsec_update, set_xs_crod

    contains

    !===============================================================================================!
    ! set xs changes due to a parameter change
    !===============================================================================================!

    subroutine set_xs_change(s, ng, nmat, ref, sigtr, siga, nuf, sigf, sigs)

        class(xs_changes)     :: s
        integer, intent(in)   :: ng ! number of energy group
        integer, intent(in)   :: nmat   ! number of material
        ! Below are material-wise XSEC changes
        real(dp), intent(in)  :: ref                 ! Reference parameter
        real(dp), intent(in)  :: sigtr(:,:)          ! Transport macroscopic XSEC
        real(dp), intent(in)  :: siga (:,:)          ! Absorption macroscopic XSEC
        real(dp), intent(in)  :: nuf  (:,:)          ! nu*fission macroscopic XSEC
        real(dp), intent(in)  :: sigf (:,:)          ! fission macroscopic XSEC
        real(dp), intent(in)  :: sigs (:,:,:)        ! Scattering macroscopic XSEC

        call xs_time % on

        allocate(s % sigtr(nmat, ng))
        allocate(s % siga (nmat, ng))
        allocate(s % nuf  (nmat, ng))
        allocate(s % sigf (nmat, ng))
        allocate(s % sigs (nmat, ng, ng))

        s % ref   = ref
        s % sigtr = sigtr
        s % siga  = siga
        s % nuf   = nuf
        s % sigf  = sigf
        s % sigs  = sigs

        call xs_time % off

    end subroutine

    !===============================================================================================!
    ! set XS related data
    !===============================================================================================!

    subroutine set_xs_data(nnod, ng, ind_mat)

        integer, intent(in)                :: ng, nnod
        integer, target, intent(in)        :: ind_mat(:)

        xsa % nnod    = nnod
        xsa % ng      = ng
        xsa % ind_mat => ind_mat

    end subroutine

    !===============================================================================================!
    ! set control rod related data
    !===============================================================================================!

    subroutine set_xs_crod(fbmap, xyz, zdel, nxx, nyy, nzz, coreh, zero_pos, step_size)

        integer, target, intent(in)   :: fbmap(:,:)             ! Radial control rod bank map (node wise)
        integer, target, intent(in)   :: xyz(:,:,:)             ! node index for given index in x, y, and z directions
        real(dp), target, intent(in)  :: zdel(:)                ! delta z
        integer, intent(in)           :: nxx, nyy, nzz
        real(dp), intent(in)          :: coreh                  ! Core Height
        real(dp), intent(in)          :: zero_pos, step_size    ! Zero step position and step size

        xsa % fbmap => fbmap
        xsa % xyz   => xyz
        xsa % zdel  => zdel
        xsa % nxx       = nxx
        xsa % nyy       = nyy
        xsa % nzz       = nzz
        xsa % coreh     = coreh
        xsa % zero_pos  = zero_pos
        xsa % step_size = step_size

    end subroutine

    !===============================================================================================!
    ! Set material wise xsec (except for chi and ADF)
    !===============================================================================================!

    subroutine xsec_setup(xs, sigtr, siga, nuf, sigf, sigs, chi, dc)

        ! Below are material-wise XSEC
        class(xs_base)                :: xs
        real(dp), target, intent(in)  :: sigtr(:,:)          ! Transport macroscopic XSEC
        real(dp), target, intent(in)  :: siga (:,:)          ! Absorption macroscopic XSEC
        real(dp), target, intent(in)  :: nuf  (:,:)          ! nu*fission macroscopic XSEC
        real(dp), target, intent(in)  :: sigf (:,:)          ! fission macroscopic XSEC
        real(dp), target, intent(in)  :: sigs (:,:,:)        ! Scattering macroscopic XSEC
        real(dp),         intent(in)  :: chi  (:,:)          ! fission macroscopic XSEC
        real(dp), target, intent(in)  :: dc   (:,:,:)        ! Scattering macroscopic XSEC

        integer :: n, nn
        integer :: ng, nnod

        call xs_time % on

        xsm % sigtr => sigtr
        xsm % siga  => siga 
        xsm % nuf   => nuf  
        xsm % sigf  => sigf 
        xsm % sigs  => sigs

        ng   = xsa % ng
        nnod = xsa % nnod

        allocate(xs % sigtr(nnod, ng), xs % siga(nnod, ng), xs % nuf(nnod, ng))
        allocate(xs % sigf(nnod, ng), xs % sigs(nnod, ng, ng), xs % chi(nnod, ng))
        allocate(xs % D(nnod, ng), xs % sigr(nnod, ng))

        if (.not. associated(xsa % ind_mat)) stop 'ind_mat not associated in xsec_setup routine'
        do n = 1, nnod
            nn = xsa % ind_mat(n)
            xs % chi (n,:) = chi (nn,:)
        end do

        select type(xs)
        type is (xs_rect)
            xs % dc => dc
        end select

        call xs_time % off

    end subroutine

    !===============================================================================================!
    !  To update current XS to base XS
    !===============================================================================================!

    subroutine xsec_update(xs, xsc, bcon, ftem, mtem, cden, bpos)

        class(xs_base)        :: xs
        class(xs_change_type) :: xsc
        real(dp), intent(in)  :: bcon
        real(dp), intent(in)  :: ftem(:)
        real(dp), intent(in)  :: mtem(:)
        real(dp), intent(in)  :: cden(:)
        real(dp), intent(in)  :: bpos(:)
        
        call xs_time % on
        
        call set_xs_base(xs)
        if (allocated(xsc % bcon)) call bcon_upd(xs, xsc % bcon, bcon)
        if (allocated(xsc % ftem)) call ftem_upd(xs, xsc % ftem, ftem)
        if (allocated(xsc % mtem)) call mod_upd (xs, xsc % mtem, mtem)
        if (allocated(xsc % cden)) call mod_upd (xs, xsc % cden, cden)
        if (allocated(xsc % crod)) call crod_upd(xs, xsc % crod, bpos)

        call upd_diff_sigr(xs)

        call xs_time % off

    end subroutine

    !===============================================================================================!
    !  To update current XS to base XS
    !===============================================================================================!

    subroutine set_xs_base(xs)

        class(xs_base)   :: xs
        
        integer   :: n, nn

        ! Set node wise xs
        do n = 1, xsa % nnod
            nn = xsa % ind_mat(n)
            xs % sigtr(n,:)   = xsm % sigtr(nn,:)
            xs % siga (n,:)   = xsm % siga (nn,:)
            xs % nuf  (n,:)   = xsm % nuf  (nn,:)
            xs % sigf (n,:)   = xsm % sigf (nn,:)
            xs % sigs (n,:,:) = xsm % sigs (nn,:,:)
        end do
        
    end subroutine

    !===============================================================================================!
    ! ! Calculate diff coeff and removal XS                                                         !
    !===============================================================================================!

    subroutine upd_diff_sigr(xs)

        type(xs_base) :: xs

        integer :: g, h, n

        xs % sigr = xs % siga

        do h = 1, xsa % ng
            do g = 1, xsa % ng
                do n = 1, xsa % nnod
                    if (g /= h) xs % sigr(n, g) = xs % sigr(n, g) + xs % sigs(n, g, h)
                end do
            end do
        end do

        do g = 1, xsa % ng
            do n = 1, xsa % nnod
                xs % D(n, g) = 1 / (3. * xs % sigtr(n, g))
            end do
        end do
        
    end subroutine

    !===============================================================================================!
    !  To update XS for given boron concentration
    !===============================================================================================!

    subroutine bcon_upd(xs, s, par)
      
        class(xs_base)        :: xs
        class(xs_changes)     :: s
        real(dp), intent(in)  :: par
        
        integer :: n, g, h

        do g = 1, xsa % ng
            do n = 1, xsa % nnod
                xs % sigtr(n,g) = xs % sigtr(n,g) + s % sigtr(xsa % ind_mat(n),g) * (par - s % ref)
                xs % siga(n,g)  = xs % siga(n,g)  + s % siga (xsa % ind_mat(n),g) * (par - s % ref)
                xs % nuf(n,g)   = xs % nuf(n,g)   + s % nuf  (xsa % ind_mat(n),g) * (par - s % ref)
                xs % sigf(n,g)  = xs % sigf(n,g)  + s % sigf (xsa % ind_mat(n),g) * (par - s % ref)
            end do
            do h = 1, xsa % ng
                do n = 1, xsa % nnod
                    xs % sigs(n,h,g) = xs % sigs(n,h,g) + s % sigs(xsa % ind_mat(n),h,g) * (par - s % ref)
                end do
            end do
        end do
      
    end subroutine

    !===============================================================================================!
    !  To update XS for given fuel temperature
    !===============================================================================================!

    subroutine ftem_upd(xs, s, par)
      
        class(xs_base)        :: xs
        class(xs_changes)     :: s
        real(dp), intent(in)  :: par(:)
        
        integer :: n, g, h

        do g = 1, xsa % ng
            do n = 1, xsa % nnod
                xs % sigtr(n,g) = xs % sigtr(n,g) + s % sigtr(xsa % ind_mat(n),g) * (sqrt(par(n)) - sqrt(s % ref))
                xs % siga(n,g)  = xs % siga(n,g)  + s % siga (xsa % ind_mat(n),g) * (sqrt(par(n)) - sqrt(s % ref))
                xs % nuf(n,g)   = xs % nuf(n,g)   + s % nuf  (xsa % ind_mat(n),g) * (sqrt(par(n)) - sqrt(s % ref))
                xs % sigf(n,g)  = xs % sigf(n,g)  + s % sigf (xsa % ind_mat(n),g) * (sqrt(par(n)) - sqrt(s % ref))
            end do
            do h = 1, xsa % ng
                do n = 1, xsa % nnod
                    xs % sigs(n,h,g) = xs % sigs(n,h,g) + s % sigs(xsa % ind_mat(n),h,g) * (sqrt(par(n)) - sqrt(s % ref))
                end do
            end do
        end do
      
    end subroutine

    !===============================================================================================!
    !  To update XS for given moderator temperature or density
    !===============================================================================================!

    subroutine mod_upd(xs, s, par)
      
        class(xs_base)        :: xs
        class(xs_changes)     :: s
        real(dp), intent(in)  :: par(:)
        
        integer :: n, g, h

        do g = 1, xsa % ng
            do n = 1, xsa % nnod
                xs % sigtr(n,g) = xs % sigtr(n,g) + s % sigtr(xsa % ind_mat(n),g) * (par(n) - s % ref)
                xs % siga(n,g)  = xs % siga(n,g)  + s % siga (xsa % ind_mat(n),g) * (par(n) - s % ref)
                xs % nuf(n,g)   = xs % nuf(n,g)   + s % nuf  (xsa % ind_mat(n),g) * (par(n) - s % ref)
                xs % sigf(n,g)  = xs % sigf(n,g)  + s % sigf (xsa % ind_mat(n),g) * (par(n) - s % ref)
            end do
            do h = 1, xsa % ng
                do n = 1, xsa % nnod
                    xs % sigs(n,h,g) = xs % sigs(n,h,g) + s % sigs(xsa % ind_mat(n),h,g) * (par(n) - s % ref)
                end do
            end do
        end do
      
    end subroutine

    !===============================================================================================!
    !  TO UPDATE AND CALCUALTE VOLUME WEIGHTED HOMOGENIZED CX FOR RODDED NODE
    !===============================================================================================!

    SUBROUTINE crod_upd(xs, s, bpos)
      
        class(xs_base)        :: xs
        class(xs_changes)     :: s
        REAL(DP), INTENT(IN)  :: bpos(:)
      
        INTEGER ::i, j, k, g, h
        REAL(DP) :: rodh, vfrac
        REAL(DP) :: dum

        integer,  pointer :: fbmap(:,:)
        integer,  pointer :: xyz(:,:,:)
        real(dp), pointer :: zdel(:)
        integer,  pointer :: ind_mat(:)

        fbmap    => xsa % fbmap
        xyz      => xsa % xyz
        zdel     => xsa % zdel
        ind_mat  => xsa % ind_mat
      
        ! For each node
        DO j = 1, xsa % nyy
            DO i = 1, xsa % nxx
                IF (fbmap(i,j) > 0) THEN
                   !!!(rodh -> posistion the tip of the control rod the top of core)
                    rodh = xsa % coreh - xsa % zero_pos  - bpos(fbmap(i,j))*xsa % step_size
                    dum = 0._dp
                    DO k = xsa % nzz, 1, -1
                        ! For partially rodded node, get volume weighted homogenized CX (0 < vfrac < 1._DP)
                        IF (rodh >= dum .AND. rodh <= dum+zdel(k)) THEN   ! If this node partially rodded
                            vfrac = (rodh - dum) / zdel(k)
                            xs % sigtr(xyz(i,j,k),:) = xs % sigtr(xyz(i,j,k),:) + &
                            vfrac * s % sigtr(ind_mat(xyz(i,j,k)),:)
                            xs % siga(xyz(i,j,k),:)  = xs % siga(xyz(i,j,k),:) + &
                            vfrac * s % siga(ind_mat(xyz(i,j,k)),:)
                            xs % nuf(xyz(i,j,k),:)   = xs % nuf(xyz(i,j,k),:) + &
                            vfrac * s % nuf(ind_mat(xyz(i,j,k)),:)
                            xs % sigf(xyz(i,j,k),:)  = xs % sigf(xyz(i,j,k),:) + &
                            vfrac * s % sigf(ind_mat(xyz(i,j,k)),:)
                            xs % sigs(xyz(i,j,k),:,:)  = xs % sigs(xyz(i,j,k),:,:) + &
                            vfrac * s % sigs(ind_mat(xyz(i,j,k)),:,:)
                          exit
                        END IF
                        ! For fully rodded node, vfrac = 1.
                        xs %sigtr(xyz(i,j,k),:)  = xs % sigtr(xyz(i,j,k),:)  + s % sigtr(ind_mat(xyz(i,j,k)),:)
                        xs %siga(xyz(i,j,k),:)   = xs % siga(xyz(i,j,k),:)   + s % siga(ind_mat(xyz(i,j,k)),:)
                        xs %nuf(xyz(i,j,k),:)    = xs % nuf(xyz(i,j,k),:)    + s % nuf(ind_mat(xyz(i,j,k)),:)
                        xs %sigf(xyz(i,j,k),:)   = xs % sigf(xyz(i,j,k),:)   + s % sigf(ind_mat(xyz(i,j,k)),:)
                        xs %sigs(xyz(i,j,k),:,:) = xs % sigs(xyz(i,j,k),:,:) + s % sigs(ind_mat(xyz(i,j,k)),:,:)
           
                        dum = dum + zdel(k)
                    END DO
                    ! if negative CX found, Surpress CX to zero and calculate D and sigr
                    DO k = xsa % nzz, 1, -1
                        DO g = 1, xsa % ng
                            IF (xs % siga(xyz(i,j,k),g) < 0.) xs % siga(xyz(i,j,k),g) = 0.
                            IF (xs % nuf(xyz(i,j,k),g) < 0.)  xs % nuf(xyz(i,j,k),g) = 0.
                            IF (xs % sigf(xyz(i,j,k),g) < 0.) xs % sigf(xyz(i,j,k),g) = 0.
                            DO h = 1, xsa % ng
                                IF (xs % sigs(xyz(i,j,k),g,h) < 0.) xs % sigs(xyz(i,j,k),g,h) = 0.
                            END DO
                        END DO
                    END DO
                END IF
            END DO
        END DO
      
      
        END SUBROUTINE
    
end module


    