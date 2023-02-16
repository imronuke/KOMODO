module xsec

    use iso_fortran_env, only: real64, error_unit
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
        logical               :: is_tabular          ! tabular xs is used
        logical               :: is_rod              ! control rod is used
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
        integer               :: nnod, ng               ! number of nodes, energy group, and delayed precursor family
        ! Below used for Control rod 
        integer, pointer      :: fbmap(:,:)             ! Radial control rod bank map (node wise)
        integer, pointer      :: xyz(:,:,:)             ! node index for given index in x, y, and z directions
        real(dp), pointer     :: zdel(:)                ! delta z
        integer               :: nxx, nyy, nzz
        real(dp)              :: coreh                  ! Core Height
        real(dp)              :: zero_pos, step_size    ! Zero step position and step size
    end type

    ! Data type for branch xsec data used if XTAB file present
    type, public :: branch_xs
        real(dp), dimension(:), allocatable :: sigtr, siga, nuf, sigf   !XSEC
        real(dp), allocatable               :: sigs(:,:)
        real(dp), allocatable               :: dc(:,:)        !ASSEMBLY DISCONTINUITY FACTOR
    end type
    ! Data Type to store data in XTAB file
    type :: branch
        integer                             :: n_den, n_boron, n_fuel, n_mod  ! Branch parameter dimension (Coolant dens., boron conc., fuel and moderator temp.)
        real(dp), dimension(:), allocatable :: pd, pb, pf, pm  ! Branch paramaters (Coolant dens., boron conc., fuel and moderator temp.)
        type(branch_xs), allocatable        :: xsec(:,:,:,:)   ! Unrodded XSEC
        type(branch_xs), allocatable        :: rxsec(:,:,:,:)  ! Rodded XSEC
        real(dp), allocatable               :: velo(:)         ! Neutron velocity
        real(dp), dimension(:), allocatable :: beta, lambda    ! beta and decay constant
        integer                             :: tadf            ! Control input: adf
        integer                             :: trod            ! Control input: control rod
    end type
    
    type(xs_mat)                            :: xsm   ! Material wise XS
    type(xs_auxilary)                       :: xsa   ! XS auxilary data used in this module
    type(branch), allocatable, public       :: m(:)  ! branch xs

    !Define + and - operators for XBRANCH type addition and substitution respectively
    interface operator (.add.)
        module procedure branch_add
    end interface
    interface operator (.subs.)
        module procedure branch_substract
    end interface
    interface operator (.mult.)
        module procedure branch_multiply
    end interface

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

    subroutine xsec_setup(xs, sigtr, siga, nuf, sigf, sigs, chi, dc, bxtab, bcrod)

        ! Below are material-wise XSEC
        class(xs_base)                :: xs
        real(dp), target, intent(in)  :: sigtr(:,:)          ! Transport macroscopic XSEC
        real(dp), target, intent(in)  :: siga (:,:)          ! Absorption macroscopic XSEC
        real(dp), target, intent(in)  :: nuf  (:,:)          ! nu*fission macroscopic XSEC
        real(dp), target, intent(in)  :: sigf (:,:)          ! fission macroscopic XSEC
        real(dp), target, intent(in)  :: sigs (:,:,:)        ! Scattering macroscopic XSEC
        real(dp),         intent(in)  :: chi  (:,:)          ! fission macroscopic XSEC
        real(dp), target, intent(in)  :: dc   (:,:,:)        ! Scattering macroscopic XSEC
        integer, intent(in)           :: bxtab               ! if tabular xs is used
        integer, intent(in)           :: bcrod               ! if control rod is used

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

        if (bxtab == 0) then
            xs % is_tabular = .false.
        else
            xs % is_tabular = .true.
        end if

        if (bcrod == 0) then
            xs % is_rod = .false.
        else
            xs % is_rod = .true.
        end if

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

        integer :: n
        
        call xs_time % on
        
        if (xs % is_tabular) then
            select type(xs)
            type is(xs_rect)
                do n = 1, xsa % nnod
                    call branch_interp(.false., xsa % ind_mat(n), cden(n),  bcon, ftem(n), mtem(n), xs % sigtr(n,:), &
                    xs % siga(n,:), xs % nuf(n,:), xs % sigf(n,:), xs % sigs(n,:,:), xs % dc(n,:,:))
                end do
            class default
                do n = 1, xsa % nnod
                    call branch_interp(.false., xsa % ind_mat(n), cden(n),  bcon, ftem(n), mtem(n), xs % sigtr(n,:), &
                    xs % siga(n,:), xs % nuf(n,:), xs % sigf(n,:), xs % sigs(n,:,:))
                end do
            end select
            if (allocated(xsc % crod)) call crod_updt_tabular(xs, bpos, cden, bcon, ftem, mtem)          
        else
            call set_xs_base(xs)
            if (allocated(xsc % bcon)) call bcon_upd(xs, xsc % bcon, bcon)
            if (allocated(xsc % ftem)) call ftem_upd(xs, xsc % ftem, ftem)
            if (allocated(xsc % mtem)) call mod_upd (xs, xsc % mtem, mtem)
            if (allocated(xsc % cden)) call mod_upd (xs, xsc % cden, cden)
            if (allocated(xsc % crod)) call crod_upd(xs, xsc % crod, bpos)
        end if

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
                    if (g /= h) xs % sigr(n,g) = xs % sigr(n,g) + xs % sigs(n,g,h)
                end do
            end do
        end do

        do g = 1, xsa % ng
            do n = 1, xsa % nnod
                xs % D(n,g) = 1 / (3. * xs % sigtr(n,g))
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

    subroutine crod_upd(xs, s, bpos)
      
        class(xs_base)        :: xs
        class(xs_changes)     :: s
        real(dp), intent(in)  :: bpos(:)
      
        integer ::i, j, k, g, h
        real(dp) :: rodh, vfrac
        real(dp) :: dum

        integer,  pointer :: fbmap(:,:)
        integer,  pointer :: xyz(:,:,:)
        real(dp), pointer :: zdel(:)
        integer,  pointer :: ind_mat(:)

        fbmap    => xsa % fbmap
        xyz      => xsa % xyz
        zdel     => xsa % zdel
        ind_mat  => xsa % ind_mat
      
        ! For each node
        do j = 1, xsa % nyy
            do i = 1, xsa % nxx
                if (fbmap(i,j) > 0) then
                   !!!(rodh -> posistion the tip of the control rod the top of core)
                    rodh = xsa % coreh - xsa % zero_pos  - bpos(fbmap(i,j))*xsa % step_size
                    dum = 0._dp
                    do k = xsa % nzz, 1, -1
                        ! For partially rodded node, get volume weighted homogenized CX (0 < vfrac < 1._DP)
                        if (rodh >= dum .AND. rodh <= dum+zdel(k)) then   ! If this node partially rodded
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
                        end if
                        ! For fully rodded node, vfrac = 1.
                        xs %sigtr(xyz(i,j,k),:)  = xs % sigtr(xyz(i,j,k),:)  + s % sigtr(ind_mat(xyz(i,j,k)),:)
                        xs %siga(xyz(i,j,k),:)   = xs % siga(xyz(i,j,k),:)   + s % siga(ind_mat(xyz(i,j,k)),:)
                        xs %nuf(xyz(i,j,k),:)    = xs % nuf(xyz(i,j,k),:)    + s % nuf(ind_mat(xyz(i,j,k)),:)
                        xs %sigf(xyz(i,j,k),:)   = xs % sigf(xyz(i,j,k),:)   + s % sigf(ind_mat(xyz(i,j,k)),:)
                        xs %sigs(xyz(i,j,k),:,:) = xs % sigs(xyz(i,j,k),:,:) + s % sigs(ind_mat(xyz(i,j,k)),:,:)
           
                        dum = dum + zdel(k)
                    end do
                    ! if negative CX found, Surpress CX to zero and calculate D and sigr
                    do k = xsa % nzz, 1, -1
                        do g = 1, xsa % ng
                            if (xs % siga(xyz(i,j,k),g) < 0.) xs % siga(xyz(i,j,k),g) = 0.
                            if (xs % nuf(xyz(i,j,k),g) < 0.)  xs % nuf(xyz(i,j,k),g) = 0.
                            if (xs % sigf(xyz(i,j,k),g) < 0.) xs % sigf(xyz(i,j,k),g) = 0.
                            do h = 1, xsa % ng
                                if (xs % sigs(xyz(i,j,k),g,h) < 0.) xs % sigs(xyz(i,j,k),g,h) = 0.
                            end do
                        end do
                    end do
                end if
            end do
        end do
      
    end subroutine

    !===============================================================================================!
    ! To interpolate the xsec data from xtab file for given bcon, ftem, mtem and cden
    !===============================================================================================!

    subroutine branch_interp (is_rod, mn, xcden, xbcon, xftem, xmtem, sigtr, siga, nuf, &
        sigf, sigs, dc)
    
        logical, intent(in)                   :: is_rod                      !
        integer, intent(in)                   :: mn                          ! material number
        real(dp), intent(in)                  :: xbcon, xftem, xmtem, xcden  ! TH Parameters
        real(dp), dimension(:), intent(out)   :: sigtr, siga, nuf, sigf
        real(dp), intent(out)                 :: sigs(:,:)
        real(dp), optional, intent(out)       :: dc(:,:)
      
        integer, parameter  :: nx = 8
        type(branch_xs)     :: xs(nx)   !Temporary xsec for interpolation
        integer             :: s, t, u, v, mx
        integer             :: s1, s2, t1, t2, u1, u2, v1, v2  ! Dimension Position of the given parameters
        integer             :: i
        real(dp)            :: radx
      
        ! Set to default
        s1=1; s2=1; t1=1; t2=1; u1=1; u2=1; v1=1; v2=1
      
        ! Get 2 closest points for interpolation
        ! FOR COOLANT DENSITY
        if (m(mn) % n_den > 1) then
            mx = m(mn) % n_den
            if (xcden >= m(mn) % pd(1) .AND. xcden <= m(mn) % pd(mx)) then
                do s = 2, mx
                    if (xcden >= m(mn) % pd(s-1) .AND. xcden <= m(mn) % pd(s)) then
                        s1 = s-1
                        s2 = s
                        exit
                    end if
                end do
            else if (xcden < m(mn) % pd(1) .AND. (m(mn) % pd(1) - xcden) / m(mn) % pd(1) < 0.2) then
                s1 = 1
                s2 = 2
            else if (xcden > m(mn) % pd(mx) .AND. (xcden - m(mn) % pd(mx)) / m(mn) % pd(mx) < 0.2) then
                s1 = mx - 1
                s2 = mx
            else
                write(error_unit,1567) xcden
                stop
                1567 format(2X, '  ERROR: COOLANT DENSITY ', F7.3 ,' g/cm3 &
                & IS OUT OF THE RANGE OF THE BRANCH PARAMETER')
            end if
        end if
      
        ! FOR BORON CONCENTRATION
        if (m(mn) % n_boron > 1) then
            mx = m(mn) % n_boron
            if (xbcon >= m(mn) % pb(1) .AND. xbcon <= m(mn) % pb(mx)) then
                do t = 2, mx
                    if (xbcon >= m(mn) % pb(t-1) .AND. xbcon <= m(mn) % pb(t)) then
                        t1 = t-1
                        t2 = t
                        EXIT
                    end if
                end do
            else if (xbcon < m(mn) % pb(1) .AND. (m(mn) % pb(1) - xbcon) < 100.) then
                t1 = 1
                t2 = 2
            else if (xbcon > m(mn) % pb(mx) .AND. (xbcon - m(mn) % pb(mx)) < 100.) then
                t1 = mx - 1
                t2 = mx
            else
                write(error_unit,1568) xbcon
                stop
                1568 format(2X, '  ERROR: BORON CONCENTRATION ', F8.1 ,' PPM &
                & IS OUT OF THE RANGE OF THE BRANCH PARAMETER')
            end if
        end if
      
        ! FOR FUEL TEMPERATURE
        if (m(mn) % n_fuel > 1) then
            mx = m(mn) % n_fuel
            if (xftem >= m(mn) % pf(1) .AND. xftem <= m(mn) % pf(mx)) then
                do u = 2, mx
                    if (xftem >= m(mn) % pf(u-1) .AND. xftem <= m(mn) % pf(u)) then
                        u1 = u-1
                        u2 = u
                        EXIT
                    end if
                end do
            else if (xftem < m(mn) % pf(1) .AND. (m(mn) % pf(1) - xftem) / m(mn) % pf(1) < 0.2) then
                u1 = 1
                u2 = 2
            else if (xftem > m(mn) % pf(mx) .AND. (xftem - m(mn) % pf(mx)) / m(mn) % pf(mx) < 0.2) then
                u1 = mx - 1
                u2 = mx
            else
                write(error_unit,1570) xftem
                stop
                1570 format(2X, '  ERROR: FUEL TEMPERATURE ', F7.1 ,' K &
                & IS OUT OF THE RANGE OF THE BRANCH PARAMETER')
            end if
        end if
      
        ! FOR MODERATOR TEMPERATURE
        if (m(mn) % n_mod > 1) then
            mx = m(mn) % n_mod
            if (xmtem >= m(mn) % pm(1) .AND. xmtem <= m(mn) % pm(mx)) then
                do v = 2, mx
                    if (xmtem >= m(mn) % pm(v-1) .AND. xmtem <= m(mn) % pm(v)) then
                        v1 = v-1
                        v2 = v
                        EXIT
                    end if
                end do
            else if (xmtem < m(mn) % pm(1) .AND. (m(mn) % pm(1) - xmtem) / m(mn) % pm(1) < 0.2) then
                v1 = 1
                v2 = 2
            else if (xmtem > m(mn) % pm(mx) .AND. (xmtem - m(mn) % pm(mx)) / m(mn) % pm(mx) < 0.2) then
                v1 = mx - 1
                v2 = mx
            else
                write(error_unit,1569) xmtem
                stop
                1569 format(2X, '  ERROR: MODERATOR TEMPERATURE ', F7.1 ,' K &
                & IS OUT OF THE RANGE OF THE BRANCH PARAMETER')
            end if
        end if
      
        !Start doing interpolation
        do i = 1, nx   !Allocate memory to xs
            allocate(xs(i) % sigtr(xsa % ng))
            allocate(xs(i) % siga (xsa % ng))
            allocate(xs(i) % nuf  (xsa % ng))
            allocate(xs(i) % sigf (xsa % ng))
            allocate(xs(i) % dc   (xsa % ng, 6))
            allocate(xs(i) % sigs (xsa % ng, xsa % ng))
        end do
      
        if (.not. is_rod) then   ! For Unrodded XSEC
            !interpolation on Moderator Temperature
            if (m(mn) % n_mod > 1) then
                radx = (xmtem - m(mn) % pm(v1)) / (m(mn) % pm(v2) - m(mn) % pm(v1))
                xs(1) = m(mn) % xsec(s1,t1,u1,v1) .add. &
                (radx .mult. (m(mn) % xsec(s1,t1,u1,v2) .subs. m(mn) % xsec(s1,t1,u1,v1)))
                xs(2) = m(mn) % xsec(s1,t1,u2,v1) .add. &
                (radx .mult. (m(mn) % xsec(s1,t1,u2,v2) .subs. m(mn) % xsec(s1,t1,u2,v1)))
                xs(3) = m(mn) % xsec(s1,t2,u1,v1) .add. &
                (radx .mult. (m(mn) % xsec(s1,t2,u1,v2) .subs. m(mn) % xsec(s1,t2,u1,v1)))
                xs(4) = m(mn) % xsec(s1,t2,u2,v1) .add. &
                (radx .mult. (m(mn) % xsec(s1,t2,u2,v2) .subs. m(mn) % xsec(s1,t2,u2,v1)))
                xs(5) = m(mn) % xsec(s2,t1,u1,v1) .add. &
                (radx .mult. (m(mn) % xsec(s2,t1,u1,v2) .subs. m(mn) % xsec(s2,t1,u1,v1)))
                xs(6) = m(mn) % xsec(s2,t1,u2,v1) .add. &
                (radx .mult. (m(mn) % xsec(s2,t1,u2,v2) .subs. m(mn) % xsec(s2,t1,u2,v1)))
                xs(7) = m(mn) % xsec(s2,t2,u1,v1) .add. &
                (radx .mult. (m(mn) % xsec(s2,t2,u1,v2) .subs. m(mn) % xsec(s2,t2,u1,v1)))
                xs(8) = m(mn) % xsec(s2,t2,u2,v1) .add. &
                (radx .mult. (m(mn) % xsec(s2,t2,u2,v2) .subs. m(mn) % xsec(s2,t2,u2,v1)))
            else
                xs(1) = m(mn) % xsec(s1,t1,u1,v1)
                xs(2) = m(mn) % xsec(s1,t1,u2,v1)
                xs(3) = m(mn) % xsec(s1,t2,u1,v1)
                xs(4) = m(mn) % xsec(s1,t2,u2,v1)
                xs(5) = m(mn) % xsec(s2,t1,u1,v1)
                xs(6) = m(mn) % xsec(s2,t1,u2,v1)
                xs(7) = m(mn) % xsec(s2,t2,u1,v1)
                xs(8) = m(mn) % xsec(s2,t2,u2,v1)
            end if
            !interpolation on Fuel Temperature
            if (m(mn) % n_fuel > 1) then
                radx = (xftem - m(mn) % pf(u1)) / (m(mn) % pf(u2) - m(mn) % pf(u1))
                xs(1) = xs(1) .add. (radx .mult. (xs(2) .subs. xs(1)))
                xs(3) = xs(3) .add. (radx .mult. (xs(4) .subs. xs(3)))
                xs(5) = xs(5) .add. (radx .mult. (xs(6) .subs. xs(5)))
                xs(7) = xs(7) .add. (radx .mult. (xs(8) .subs. xs(7)))
            end if
            !interpolation on Boron concentration
            if (m(mn) % n_boron > 1) then
                radx = (xbcon - m(mn) % pb(t1)) / (m(mn) % pb(t2) - m(mn) % pb(t1))
                xs(1) = xs(1) .add. (radx .mult. (xs(3) .subs. xs(1)))
                xs(5) = xs(5) .add. (radx .mult. (xs(7) .subs. xs(5)))
            end if
            !interpolation on coolant density
            if (m(mn) % n_den > 1) then
                xs(1) = xs(1) .add. ((xcden - m(mn) % pd(s1)) / (m(mn) % pd(s2) - m(mn) % pd(s1)) .mult. &
                (xs(5) .subs. xs(1)))
            end if
        else   ! For Rodded XSEC
            !interpolation on Moderator Temperature
            if (m(mn) % n_mod > 1) then
                radx = (xmtem - m(mn) % pm(v1)) / (m(mn) % pm(v2) - m(mn) % pm(v1))
                xs(1) = m(mn) % rxsec(s1,t1,u1,v1) .add. &
                (radx .mult. (m(mn) % rxsec(s1,t1,u1,v2) .subs. m(mn) % rxsec(s1,t1,u1,v1)))
                xs(2) = m(mn) % rxsec(s1,t1,u2,v1) .add. &
                (radx .mult. (m(mn) % rxsec(s1,t1,u2,v2) .subs. m(mn) % rxsec(s1,t1,u2,v1)))
                xs(3) = m(mn) % rxsec(s1,t2,u1,v1) .add. &
                (radx .mult. (m(mn) % rxsec(s1,t2,u1,v2) .subs. m(mn) % rxsec(s1,t2,u1,v1)))
                xs(4) = m(mn) % rxsec(s1,t2,u2,v1) .add. &
                (radx .mult. (m(mn) % rxsec(s1,t2,u2,v2) .subs. m(mn) % rxsec(s1,t2,u2,v1)))
                xs(5) = m(mn) % rxsec(s2,t1,u1,v1) .add. &
                (radx .mult. (m(mn) % rxsec(s2,t1,u1,v2) .subs. m(mn) % rxsec(s2,t1,u1,v1)))
                xs(6) = m(mn) % rxsec(s2,t1,u2,v1) .add. &
                (radx .mult. (m(mn) % rxsec(s2,t1,u2,v2) .subs. m(mn) % rxsec(s2,t1,u2,v1)))
                xs(7) = m(mn) % rxsec(s2,t2,u1,v1) .add. &
                (radx .mult. (m(mn) % rxsec(s2,t2,u1,v2) .subs. m(mn) % rxsec(s2,t2,u1,v1)))
                xs(8) = m(mn) % rxsec(s2,t2,u2,v1) .add. &
                (radx .mult. (m(mn) % rxsec(s2,t2,u2,v2) .subs. m(mn) % rxsec(s2,t2,u2,v1)))
            else
                xs(1) = m(mn) % rxsec(s1,t1,u1,v1)
                xs(2) = m(mn) % rxsec(s1,t1,u2,v1)
                xs(3) = m(mn) % rxsec(s1,t2,u1,v1)
                xs(4) = m(mn) % rxsec(s1,t2,u2,v1)
                xs(5) = m(mn) % rxsec(s2,t1,u1,v1)
                xs(6) = m(mn) % rxsec(s2,t1,u2,v1)
                xs(7) = m(mn) % rxsec(s2,t2,u1,v1)
                xs(8) = m(mn) % rxsec(s2,t2,u2,v1)
            end if
            !interpolation on Fuel Temperature
            if (m(mn) % n_fuel > 1) then
                radx = (xftem - m(mn) % pf(u1)) / (m(mn) % pf(u2) - m(mn) % pf(u1))
                xs(1) = xs(1) .add. (radx .mult. (xs(2) .subs. xs(1)))
                xs(3) = xs(3) .add. (radx .mult. (xs(4) .subs. xs(3)))
                xs(5) = xs(5) .add. (radx .mult. (xs(6) .subs. xs(5)))
                xs(7) = xs(7) .add. (radx .mult. (xs(8) .subs. xs(7)))
            end if
            !interpolation on Boron concentration
            if (m(mn) % n_boron > 1) then
                radx = (xbcon - m(mn) % pb(t1)) / (m(mn) % pb(t2) - m(mn) % pb(t1))
                xs(1) = xs(1) .add. (radx .mult. (xs(3) .subs. xs(1)))
                xs(5) = xs(5) .add. (radx .mult. (xs(7) .subs. xs(5)))
            end if
        
            !interpolation on coolant density
            if (m(mn) % n_den > 1) then
                xs(1) = xs(1) .add. ((xcden - m(mn) % pd(s1)) / (m(mn) % pd(s2) - m(mn) % pd(s1)) .mult. &
                (xs(5) .subs. xs(1)))
            end if
        end if
        sigtr = xs(1) % sigtr
        siga  = xs(1) % siga
        nuf   = xs(1) % nuf
        sigf  = xs(1) % sigf
        sigs  = xs(1) % sigs
        dc    = xs(1) % dc
      
        do i = 1, nx   !DeAllocate memory to xs
            deallocate(xs(i) % sigtr)
            deallocate(xs(i) % siga)
            deallocate(xs(i) % nuf)
            deallocate(xs(i) % sigf)
            deallocate(xs(i) % dc)
            deallocate(xs(i) % sigs)
        end do
    
    end subroutine

    !===============================================================================================!
    ! TO UPDATE AND CALCUALTE VOLUME WEIGHTED HOMOGENIZED CX FOR RODDED NODE (WITH TABULAR XS)
    !===============================================================================================!

    subroutine crod_updt_tabular (xs, bpos, cden, bcon, ftem, mtem)
      
        class(xs_base)        :: xs
        real(dp), intent(in) :: bpos(:)
        real(dp), intent(in) :: cden(:)
        real(dp), intent(in) :: bcon
        real(dp), intent(in) :: ftem(:)
        real(dp), intent(in) :: mtem(:)
      
        integer  :: i, j, k, g, h
        real(dp) :: rodh, vfrac
        real(dp) :: dum
      
        integer :: n
        integer,  pointer :: fbmap(:,:)
        integer,  pointer :: xyz(:,:,:)
        real(dp), pointer :: zdel(:)
        integer,  pointer :: ind_mat(:)
      
        real(dp), dimension(xsa % nnod, xsa % ng) :: rsigtr, rsiga, rnuf, rsigf
        real(dp)                                  :: rsigs(xsa % nnod, xsa % ng, xsa % ng)
        real(dp)                                  :: rdc(xsa % nnod, xsa % ng, 6)


        fbmap    => xsa % fbmap
        xyz      => xsa % xyz
        zdel     => xsa % zdel
        ind_mat  => xsa % ind_mat
      
        do j = 1, xsa % nyy
            do i = 1, xsa % nxx
                if (fbmap(i,j) > 0) then
                    !(rodh -> posistion the tip of the control rod the top of core)
                    rodh = xsa % coreh - xsa % zero_pos - bpos(fbmap(i,j))*xsa % step_size
                    dum = 0._DP
                    do k = xsa % nzz, 1, -1
                        n = xyz(i,j,k)                                  ! Node number
                        if (m(ind_mat(n)) % trod == 1) then
                            call branch_interp(.true., ind_mat(n), cden(n), bcon, ftem(n), &
                            mtem(n), rsigtr(n,:), rsiga(n,:), rnuf(n,:), &
                            rsigf(n,:), rsigs(n,:,:), rdc(n,:,:))
                        else
                            write(error_unit,1671) fbmap(i,j), ind_mat(n)
                            stop
                            1671 format (2X, 'CONTROL ROD BANK NUMBER ', I4, &
                            ' COINCIDES WITH MATERIAL NUMBER ', I4, &
                            ' THAT DOES NOT HAVE CONTROL ROD DATA IN XTAB FILE')
                        end if
                        ! For partially rodded node, get volume weighted homogenized CX (0 < vfrac < 1._DP)
                        if (rodh >= dum .AND. rodh <= dum+zdel(k)) then   ! If this node partially rodded
                            vfrac = (rodh - dum) / zdel(k)                    ! Rodded fraction
                            xs % sigtr(n,:)  = (1. - vfrac) * xs % sigtr(n,:)  + vfrac * rsigtr(n,:)
                            xs % siga(n,:)   = (1. - vfrac) * xs % siga(n,:)   + vfrac * rsiga(n,:)
                            xs % nuf(n,:)    = (1. - vfrac) * xs % nuf(n,:)    + vfrac * rnuf(n,:)
                            xs % sigf(n,:)   = (1. - vfrac) * xs % sigf(n,:)   + vfrac * rsigf(n,:)
                            xs % sigs(n,:,:) = (1. - vfrac) * xs % sigs(n,:,:) + vfrac * rsigs(n,:,:)
                            select type(xs)
                            type is (xs_rect)
                                xs % dc(n,:,:)   = (1. - vfrac) * xs % dc(n,:,:)   + vfrac * rdc(n,:,:)
                            end select
                            exit
                        end if
                        ! For fully rodded node, vfrac = 1.
                        xs % sigtr(n,:)  = rsigtr(n,:)
                        xs % siga(n,:)   = rsiga(n,:)
                        xs % nuf(n,:)    = rnuf(n,:)
                        xs % sigf(n,:)   = rsigf(n,:)
                        xs % sigs(n,:,:) = rsigs(n,:,:)
                        select type(xs)
                        type is (xs_rect)
                            xs % dc(n,:,:)   = rdc(n,:,:)
                        end select
      
                        dum = dum + zdel(k)
                    end do
                    ! if negative CX found, Surpress CX to zero  and calculate D and sigr
                    do k = xsa % nzz, 1, -1
                        n = xyz(i,j,k)
                        do g = 1, xsa % ng
                            if (xs % siga(n,g) < 0.) xs % siga(n,g)  = 0.
                            if (xs % nuf(n,g)  < 0.) xs % nuf(n,g)   = 0.
                            if (xs % sigf(n,g) < 0.) xs % sigf(n,g)  = 0.
                            do h = 1, xsa % ng
                                if (xs % sigs(n,g,h) < 0.) xs % sigs( n,g,h) = 0.
                            end do
                            select type(xs)
                            type is (xs_rect)
                                do h = 1, 6
                                    if (xs % dc(n,g,h) < 0.) xs % dc(n,g,h) = 0.
                                end do
                            end select
                        end do
                    end do
                end if
            end do
        end do
      
    end subroutine
    
    !===============================================================================================!
    ! To perform branch data type addition
    !===============================================================================================!
    
    function branch_add(A, B) result (C)
        
        type(branch_xs), intent(in) :: A, B
        type(branch_xs)             :: C

        C % sigtr = A % sigtr + B % sigtr
        C % siga  = A % siga  + B % siga
        C % nuf   = A % nuf   + B % nuf
        C % sigf  = A % sigf  + B % sigf
        C % sigs  = A % sigs  + B % sigs
        C % dc    = A % dc    + B % dc
    
    end function
    
    !===============================================================================================!
    ! To perform branch data type substraction
    !===============================================================================================!
    
    function branch_substract(A, B) result (C)
    
        type(branch_xs), intent(in) :: A, B
        type(branch_xs)             :: C

        C % sigtr = A % sigtr - B % sigtr
        C % siga  = A % siga  - B % siga
        C % nuf   = A % nuf   - B % nuf
        C % sigf  = A % sigf  - B % sigf
        C % sigs  = A % sigs  - B % sigs
        C % dc    = A % dc    - B % dc
    
    end function
    
    !===============================================================================================!
    ! ! To perform branch data type substraction
    !===============================================================================================!
    
    function branch_multiply(x, A) result (B)
    
        real(dp), intent(in)     :: x
        type(branch_xs), intent(in) :: A
        type(branch_xs)             :: B

        B % sigtr = x * A % sigtr
        B % siga  = x * A % siga
        B % nuf   = x * A % nuf
        B % sigf  = x * A % sigf
        B % sigs  = x * A % sigs
        B % dc    = x * A % dc
    
    end function
    
end module


    