module fdm

    use omp_lib
    use iso_fortran_env, only: real64, output_unit
    use time
    use xsec
    use stagger
    use linear_system
    use utilities
    use nodal

    implicit none

    private

    save

    integer, parameter :: dp = real64

    integer, parameter  :: n_rect_surf = 6           ! order is => east, west, north, top, bottom

    integer, parameter  :: zero_flux = 0
    integer, parameter  :: zero_incoming = 1
    integer, parameter  :: reflective = 2

    integer :: fid                          ! file unit number

    integer :: xwest, xeast                 ! West and East Boundary conditions
    integer :: ysouth, ynorth               ! South and North Boundary conditions
    integer :: zbott, ztop                  ! Bottom and Top Boundary conditions

    ! Finite Difference Method (fdm) data type
    type, public :: fdm_base
        character(4)        :: kernel        ! Nodal kernel used
        integer             :: nnod          ! Number of nodes
        real(dp), pointer   :: vdel(:)       ! Delta nodes' volume in cm3
    
        integer             :: ng            ! number of neutron energy groups

        integer             :: nmat          ! number of materials
    
        real(dp)            :: Keff = 1.0    ! Multiplication factor

        real(dp), pointer     :: fsrc(:) => null()     ! fission source
        real(dp), pointer     :: flux(:,:) => null()   ! flux
        real(dp), pointer     :: exsrc(:,:) => null()  ! external source
        real(dp), allocatable :: scat(:,:)             ! scattering source
    
        type(xs_rect)       :: xs            ! node-wise xs
        type(matrix)        :: m             ! FDM matrix for linear system solver

        real(dp)            :: fsrc_diff     ! Fission source difference between outer iteration
        real(dp)            :: flux_diff     ! Flux difference between outer iteration

        ! used only for transient
        real(dp), pointer   :: total_beta(:) => null()         ! total beta for each material
        real(dp), pointer   :: dfis(:) => null()  

        integer             :: nod_intv      ! nodal coupling coefficients update interval
        real(dp)            :: max_flux_err  ! max flux error
        real(dp)            :: max_fiss_err  ! max fission source error
        integer             :: max_outer     ! max number of outer iteration
        integer             :: max_inner     ! max number of outer iteration
        integer             :: extrp         ! interval for fission source extrapolation
    end type

    ! Rectangular FDM data type
    type, extends(fdm_base), public :: fdm_rect
        real(dp), pointer :: xdel(:) => null()        ! Delta x in cm
        real(dp), pointer :: ydel(:) => null()        ! Delta y in cm
        real(dp), pointer :: zdel(:) => null()        ! Delta z in cm
    
        integer  :: nxx                               ! Number of nodes in x direction
        integer  :: nyy                               ! Number of nodes in y direction
        integer  :: nzz                               ! Number of nodes z direction

        integer, pointer :: ix(:) => null()           ! index in x direction for given node index
        integer, pointer :: iy(:) => null()           ! index in y direction for given node index
        integer, pointer :: iz(:) => null()           ! index in z direction for given node index
        integer, pointer :: xyz(:,:,:) => null()      ! node index for given index in x, y, and z directions

        integer, pointer :: ind_mat(:) => null()      ! material map / node-wise material index
        
        ! coupling coeff (order => east, west, north, south, top, bottom)
        real(dp), allocatable      :: df(:,:,:)   ! D tilde : FDM nodal coupling coefficients
        real(dp), allocatable      :: dn(:,:,:)   ! D hat   : Corrected (higher order) nodal coupling coefficients

        type(staggered), pointer   :: xstag(:) => null()    ! Staggered mesh data
        type(staggered), pointer   :: ystag(:) => null()    ! Staggered mesh data

        type(nodal_type), allocatable :: d        ! nodal update module

    end type

    public :: set_fdm_data, set_fdm_pointer, set_fdm_transient, deallocate_fdm
    public :: outer_iter, outer_adjoint, outer_fixed_src, th_outer, transient_outer

    contains

    !===============================================================================================!
    ! Set FDM data for rectangular geometry
    !===============================================================================================!

    subroutine set_fdm_data(c, output_fid, kernel, ng, nnod, nxx, nyy, nzz, &
        nmat, east, west, north, south, top, bottom, nod_intv, max_flux_err, &
        max_fiss_err, max_outer, max_inner, extrp)

        class(fdm_rect)              :: c
        integer, intent(in)           :: output_fid                               ! output file unit number
        character(*), intent(in)      :: kernel
        integer, intent(in)           :: nnod                                     ! total number of nodes
        integer, intent(in)           :: ng                                       ! number of group
        integer, intent(in)           :: nxx, nyy, nzz                            ! number of nodes in x, y, and z directions
        integer, intent(in)           :: nmat                                     ! number of material
        integer, intent(in)           :: east, west, north, south, bottom, top    ! BCs
        integer, intent(in)           :: nod_intv      ! nodal coupling coefficients update interval
        real(dp), intent(in)          :: max_flux_err  ! max flux error
        real(dp), intent(in)          :: max_fiss_err  ! max fission source error
        integer, intent(in)           :: max_outer     ! max number of outer iteration
        integer, intent(in)           :: max_inner     ! max number of outer iteration
        integer, intent(in)           :: extrp         ! interval for fission source extrapolation

        !set output file unit number
        fid = output_fid

        !set nodal kernel
        c % kernel = kernel
        
        ! set number of energy group
        c % ng = ng

        ! set number of nodes in x, y, and z directions
        c % nnod = nnod
        c % nxx  = nxx
        c % nyy  = nyy
        c % nzz  = nzz
        c % nmat = nmat

        ! set BCs
        xeast  = east
        xwest  = west
        ynorth = north
        ysouth = south
        zbott  = bottom
        ztop   = top

        ! set iterations
        c % nod_intv     = nod_intv
        c % max_flux_err = max_flux_err
        c % max_fiss_err = max_fiss_err
        c % max_outer    = max_outer   
        c % max_inner    = max_inner   
        c % extrp        = extrp
        
        ! allocate scattering source
        allocate(c % scat(nnod, ng))
        
    end subroutine

    !===============================================================================================!
    ! Set FDM pointers for rectangular geometry
    !===============================================================================================!

    subroutine set_fdm_pointer(c, flux, fsrc, exsrc, ind_mat, xdel, ydel, zdel, vdel, &
        xstag, ystag, ix, iy, iz, xyz)

        class(fdm_rect)       :: c
        real(dp), intent(in), target           :: flux(:,:)
        real(dp), intent(in), target           :: fsrc(:), exsrc(:,:)
        integer, intent(in), target            :: ind_mat(:)      ! 3D Material map
        real(dp), intent(in), target           :: xdel(:), ydel(:), zdel(:), vdel(:)
        integer, intent(in), target            :: ix(:), iy(:), iz(:), xyz(:,:,:)
        type(staggered), intent(in), target    :: xstag(:), ystag(:)

        c % flux => flux
        c % flux = 1.0  ! Initialize flux
        if (size(c % flux, dim=1) .ne. c % nnod) stop 'Flux first-dimension size does not match'
        if (size(c % flux, dim=2) .ne. c % ng)   stop 'Flux second-dimension size does not match'

        c % fsrc => fsrc
        if (size(c % fsrc) .ne. c % nnod) stop 'Fission source size does not match'

        c % exsrc => exsrc
        if (size(c % exsrc, dim=1) .ne. c % nnod) stop 'exsrc first-dimension size does not match'
        if (size(c % exsrc, dim=2) .ne. c % ng)   stop 'exsrc second-dimension size does not match'

        c % ind_mat => ind_mat
        if (size(c % ind_mat) .ne. c % nnod) stop 'material index ind_mat size does not match'

        c % xdel => xdel
        c % ydel => ydel
        c % zdel => zdel
        c % vdel => vdel
        c % ix => ix
        c % iy => iy
        c % iz => iz
        
        c % xyz => xyz
        if (size(c % xyz, dim=1) .ne. c % nxx) stop 'xyz first-dimension size does not match'
        if (size(c % xyz, dim=2) .ne. c % nyy) stop 'xyz second-dimension size does not match'
        if (size(c % xyz, dim=3) .ne. c % nzz) stop 'xyz third-dimension size does not match'
        
        c % xstag => xstag
        c % ystag => ystag

        ! final FDM setup
        call setup_rectangular(c)

    end subroutine

    !===============================================================================================!
    ! Set FDM pointers for rectangular geometry
    !===============================================================================================!

    subroutine set_fdm_transient(c, total_beta, dfis)

        class(fdm_rect)       :: c
        real(dp), intent(in), target        :: dfis(:)
        real(dp), intent(in), target        :: total_beta(:)

        c % dfis => dfis
        c % total_beta => total_beta

        if (c % kernel .ne. ' FDM') then
            call set_nodal_transient(c % d, c % ind_mat, c % total_beta, c % dfis)
        end if

    end subroutine

    !===============================================================================================!
    ! Set FDM pointers for rectangular geometry
    !===============================================================================================!

    subroutine deallocate_fdm(c)

        class(fdm_rect) :: c

        deallocate(c % scat)
        deallocate(c % df, c % dn)
        if (allocated(c % d)) deallocate(c % d)

        c % flux => null()
        c % fsrc => null()
        c % exsrc => null()
        c % ind_mat => null()
        c % xdel => null()
        c % ydel => null()
        c % zdel => null()
        c % vdel => null()
        c % ix => null()
        c % iy => null()
        c % iz => null()
        c % xyz => null()
        c % xstag => null()
        c % ystag => null()

        c % dfis => null()
        c % total_beta => null()


    end subroutine

    !===============================================================================================!
    ! to perfom forward outer iteration
    !===============================================================================================!

    subroutine outer_iter(c, print_iter, is_converge)

        class(fdm_rect)               :: c
        logical, intent(in)            :: print_iter
        logical, intent(out)           :: is_converge

        integer  :: i, g
        real(dp) :: fc, fco                                  ! new and old integrated fission SOURCES
        real(dp) :: Keo
        real(dp) :: flux_old(c % nnod, c % ng)
        real(dp) :: fsrc_old(c % nnod)
        real(dp) :: bs(c % nnod)
        real(dp) :: e1, e2
        real(dp) :: errn(c % nnod), erro(c % nnod)            ! current and past error vectors

        call fdm_time % on

        ! Build FDM matrix
        call build_matrix(c, upd_fdm_coupling=.true.)

        ! Initialize fission source
        call get_fission_src(c)
        fc = integrate(c, c % fsrc)

        ! Init error vector
        errn = 1._dp
        e1 = integrate(c, errn)

        call fdm_time % off
    
        !Start outer iteration
        is_converge    = .false.
        do i = 1, c % max_outer

            call fdm_time % on
            fco       = fc               ! Save old integrated fission source
            Keo       = c % Keff         ! Save old keff
            flux_old  = c % flux         ! Save old flux
            fsrc_old  = c % fsrc         ! Save old fission source
            erro      = errn             ! Save old fission source error/difference
            
            ! Inner iteration
            do g = 1, c % ng
                call get_total_src(c, g, bs)
                call linear_solver(c % m % mtx(g), c % m % ind, c % max_inner, bs, c % flux(:,g))
            end do

            call get_fission_src(c)                      !Update fission source
            errn = c % fsrc - fsrc_old
            e2 = norm2(errn)
            if (mod(i, c % extrp) == 0) then
                call fsrc_extrp(e1, e2, erro, errn, c % fsrc, print_iter) ! Fission source extrapolation
            end if
            e1 = e2
            fc = integrate(c, c % fsrc)
            c % Keff = Keo * fc / fco             ! Update Keff
            call get_difference(c % flux, flux_old, c % fsrc, fsrc_old, c % flux_diff, c % fsrc_diff)
            if (print_iter) then
                write(fid,'(I5,F13.6,2ES15.5)') i, c % Keff, c % fsrc_diff, c % flux_diff
                write(output_unit,'(I5,F13.6,2ES15.5)')   i, c % Keff, c % fsrc_diff, c % flux_diff
            end if
            if ((c % max_flux_err > c % flux_diff) .and. (c % max_fiss_err > c % fsrc_diff)) then
                is_converge = .true.
                exit
            end if
            call fdm_time % off

            if (mod(i, c % nod_intv) == 0) then
                call nodal_update(c, 'static', print_iter)  ! Nodal coefficients update
            end if
        end do

    end subroutine

    !===============================================================================================!
    ! To perform transient fixed source outer iteration
    !===============================================================================================!

    subroutine transient_outer(c, mod_sigr, is_converge)

        class(fdm_rect)               :: c
        real(dp), intent(in)           :: mod_sigr(:,:)  ! modified sigma removal for transient
        logical, intent(out)           :: is_converge

        integer  :: i, g
        real(dp) :: flux_old(c % nnod, c % ng)
        real(dp) :: fsrc_old(c % nnod)
        real(dp) :: bs(c % nnod)
        real(dp) :: e1, e2
        real(dp) :: errn(c % nnod), erro(c % nnod)            ! current and past error vectors

        call fdm_time % on

        ! Build FDM matrix
        call build_matrix(c, upd_fdm_coupling=.true., sigr=mod_sigr)

        ! Init error vector
        errn = 1._dp
        e1 = integrate(c, errn)

        call fdm_time % off
    
        !Start outer iteration
        is_converge = .false.
        do i = 1, c % max_outer

            call fdm_time % on
            flux_old  = c % flux         ! Save old flux
            fsrc_old  = c % fsrc         ! Save old fission source
            erro      = errn             ! Save old fission source error/difference
            
            ! Inner iteration
            do g = 1, c % ng
                call get_total_transient(c, g, bs)
                call linear_solver(c % m % mtx(g), c % m % ind, c % max_inner, bs, c % flux(:,g))
            end do

            call get_fission_src(c)                      !Update fission source
            errn = c % fsrc - fsrc_old
            e2 = norm2(errn)
            if (mod(i, c % extrp) == 0) then
                call fsrc_extrp(e1, e2, erro, errn, c % fsrc, .false.) ! Fission source extrapolation
            end if
            e1 = e2
            call get_difference(c % flux, flux_old, c % fsrc, fsrc_old, c % flux_diff, c % fsrc_diff)
            if ((c % max_flux_err > c % flux_diff) .AND. (c % max_fiss_err > c % fsrc_diff)) then
                is_converge = .true.
                exit
            end if
            call fdm_time % off

            if (mod(i, c % nod_intv) == 0) then
                call nodal_update(c, 'transient', print_iter=.false., sigr=mod_sigr)  ! Nodal coefficients update
            end if
        end do

    end subroutine

    !===============================================================================================!
    ! to perfom outer iteration (for TH feedback)
    ! Has fixed number of outer iterations
    !===============================================================================================!

    subroutine th_outer(c)

        class(fdm_rect)               :: c

        integer  :: i, g
        real(dp) :: fc, fco                                  ! new and old integrated fission SOURCES
        real(dp) :: Keo
        real(dp) :: flux_old(c % nnod, c % ng)
        real(dp) :: fsrc_old(c % nnod)
        real(dp) :: bs(c % nnod)
        real(dp) :: e1, e2
        real(dp) :: errn(c % nnod), erro(c % nnod)            ! current and past error vectors

        call fdm_time % on

        ! Build FDM matrix
        call build_matrix(c, upd_fdm_coupling=.true.)

        c % nod_intv = ceiling(0.45 * real(c % max_outer))

        ! Initialize fission source
        call get_fission_src(c)
        fc = integrate(c, c % fsrc)

        ! Init error vector
        errn = 1._dp
        e1 = integrate(c, errn)

        call fdm_time % off
    
        !Start outer iteration
        do i = 1, c % max_outer

            call fdm_time % on
            fco       = fc               ! Save old integrated fission source
            Keo       = c % Keff         ! Save old keff
            flux_old  = c % flux         ! Save old flux
            fsrc_old  = c % fsrc         ! Save old fission source
            erro      = errn             ! Save old fission source error/difference
            
            ! Inner iteration
            do g = 1, c % ng
                call get_total_src(c, g, bs)
                call linear_solver(c % m % mtx(g), c % m % ind, c % max_inner, bs, c % flux(:,g))
            end do

            ! if (any(c % flux < 0.)) then
            !     write(*,*) i
            !     stop
            ! end if

            call get_fission_src(c)      !Update fission source
            errn = c % fsrc - fsrc_old
            e2 = norm2(errn)
            if (mod(i, c % extrp) == 0) then
                call fsrc_extrp(e1, e2, erro, errn, c % fsrc, .false.) ! Fission source extrapolation
            end if
            e1 = e2
            fc = integrate(c, c % fsrc)
            c % Keff = Keo * fc / fco    ! Update Keff
            write(*,*) c % Keff
            call get_difference(c % flux, flux_old, c % fsrc, fsrc_old, c % flux_diff, c % fsrc_diff)
            call fdm_time % off

            if (mod(i, c % nod_intv) == 0) then
                call nodal_update(c, 'static', print_iter=.false.)  ! Nodal coefficients update
            end if
        end do

    end subroutine

    !===============================================================================================!
    ! to perfom outer iteration  (fixed source)
    !===============================================================================================!

    subroutine outer_fixed_src(c, print_iter, is_converge)
    
        class(fdm_rect)               :: c
        logical, intent(in)            :: print_iter
        logical, intent(out)           :: is_converge

        integer  :: i, g
        real(dp) :: flux_old(c % nnod, c % ng)
        real(dp) :: fsrc_old(c % nnod)
        real(dp) :: bs(c % nnod)
        real(dp) :: e1, e2
        real(dp) :: errn(c % nnod), erro(c % nnod)            ! current and past error vectors

        call fdm_time % on

        ! Build FDM matrix
        call build_matrix(c, upd_fdm_coupling=.true.)
    
        ! Initialize fission source
        call get_fission_src(c)

        ! Init error vector
        errn = 1._dp
        e1 = integrate(c, errn)

        call fdm_time % off
    
        !Start outer iteration
        is_converge = .false.
        do i = 1, c % max_outer

            call fdm_time % on
            flux_old  = c % flux          ! Save old flux
            fsrc_old  = c % fsrc         ! Save old fission source
            erro       = errn             ! Save old fission source error/difference
            
            ! Inner iteration
            do g = 1, c % ng
                call get_total_src(c, g, bs)
                call linear_solver(c % m % mtx(g), c % m % ind, c % max_inner, bs, c % flux(:,g))
            end do

            call get_fission_src(c)                      !Update fission source
            errn = c % fsrc - fsrc_old
            e2 = norm2(errn)
            if (mod(i, c % extrp) == 0) then
                call fsrc_extrp(e1, e2, erro, errn, c % fsrc, print_iter) ! Fission source extrapolation
            end if
            e1 = e2
            call get_difference(c % flux, flux_old, c % fsrc, fsrc_old, c % flux_diff, c % fsrc_diff)
            if (print_iter) then
                write(fid,'(I5,2ES15.5)') i, c % fsrc_diff, c % flux_diff
                write(output_unit,'(I5,2ES15.5)')   i, c % fsrc_diff, c % flux_diff
            end if
            if ((c % max_flux_err > c % flux_diff) .AND. (c % max_fiss_err > c % fsrc_diff)) then
                is_converge = .true.
                exit
            end if
            call fdm_time % off

            if (mod(i, c % nod_intv) == 0) then
                call nodal_update(c, 'static', print_iter)  ! Nodal coefficients update
            end if
        end do
    
    end subroutine

    !===============================================================================================!
    ! to perfom outer iteration (adjoint)
    !===============================================================================================!

    subroutine outer_adjoint(c, print_iter, is_converge)
    
        class(fdm_rect)               :: c
        logical, intent(in)            :: print_iter
        logical, intent(out)           :: is_converge

        integer  :: i, g
        real(dp) :: fc, fco                                  ! new and old integrated fission SOURCES
        real(dp) :: Keo
        real(dp) :: flux_old(c % nnod, c % ng)
        real(dp) :: fsrc_old(c % nnod)
        real(dp) :: bs(c % nnod)
        real(dp) :: e1, e2
        real(dp) :: errn(c % nnod), erro(c % nnod)            ! current and past error vectors

        call fdm_time % on

        ! Build FDM matrix
        call build_matrix(c, upd_fdm_coupling=.true.)
    
        ! Initialize fission source
        call get_fission_adjoint(c)
        fc = integrate(c, c % fsrc)

        ! Init error vector
        errn = 1._dp
        e1   = integrate(c, errn)

        call fdm_time % off
    
        !Start outer iteration
        is_converge = .false.
        do i = 1, c % max_outer

            call fdm_time % on
            fco        = fc               ! Save old integrated fission source
            Keo        = c % Keff         ! Save old keff
            flux_old  = c % flux          ! Save old flux
            fsrc_old  = c % fsrc          ! Save old fission source
            erro       = errn             ! Save old fission source error/difference
            
            ! Inner iteration
            do g = c % ng, 1, -1
                call get_total_adjoint(c, g, bs)
                call linear_solver(c % m % mtx(g), c % m % ind, c % max_inner, bs, c % flux(:,g))
            end do

            call get_fission_adjoint(c)                      !Update fission source
            errn = c % fsrc - fsrc_old
            e2 = norm2(errn)
            if (mod(i, c % extrp) == 0) then
                call fsrc_extrp(e1, e2, erro, errn, c % fsrc, print_iter) ! Fission source extrapolation
            end if
            e1 = e2
            fc = integrate(c, c % fsrc)
            c % Keff = Keo * fc / fco             ! Update Keff
            call get_difference(c % flux, flux_old, c % fsrc, fsrc_old, c % flux_diff, c % fsrc_diff)
            if (print_iter) then
                write(fid,'(I5,F13.6,2ES15.5)') i, c % Keff, c % fsrc_diff, c % flux_diff
                write(output_unit,'(I5,F13.6,2ES15.5)')   i, c % Keff, c % fsrc_diff, c % flux_diff
            end if
            if ((c % max_flux_err > c % flux_diff) .AND. (c % max_fiss_err > c % fsrc_diff)) then
                is_converge = .true.
                exit
            end if
            call fdm_time % off

            if (mod(i, c % nod_intv) == 0) then
                call nodal_update(c, 'adjoint', print_iter)  ! Nodal coefficients update
            end if
        end do
    
    end subroutine

    !===============================================================================================!
    ! To update nodal coupling coefficients (D hat)
    !===============================================================================================!

    subroutine nodal_update(c, calc_mode, print_iter, sigr)
      
        class(fdm_rect), intent(in)   :: c
        character(*), intent(in)      :: calc_mode
        logical, intent(in)           :: print_iter
        real(dp), optional            :: sigr(:,:)

        call nodal_time % on
      
        ! Update nodal coupling coefficients
        if (c % kernel == ' FDM') then
            return
        elseif (c % kernel == 'SANM') then
            call nodal_update_sanm(c % d, calc_mode)
        else
            call nodal_update_pnm(c % d, calc_mode)
        end if

        ! Update CMFD matrix
        if (present(sigr)) then
            call build_matrix(c, upd_fdm_coupling=.false., sigr=sigr)
        else
            call build_matrix(c, upd_fdm_coupling=.false.)
        end if

        if (print_iter) then
            write(fid,*) '    .....NODAL COUPLING UPDATED..... '
            write(output_unit,*) '    .....NODAL COUPLING UPDATED..... '
        end if

        call nodal_time % off
      
      end subroutine

    !===============================================================================================!
    ! Search maximum point wise fission source Relative difference, and
    ! search maximum point wise flux relative difference
    !===============================================================================================!

    subroutine get_difference(flux, flux_last, fsrc, fsrc_last, flux_diff, fsrc_diff)

        real(dp), intent(in)  :: flux(:,:), flux_last(:,:)
        real(dp), intent(in)  :: fsrc(:), fsrc_last(:)
        real(dp), intent(out) :: flux_diff, fsrc_diff

        integer  :: nnod, ng
        integer  :: n, g
        real(dp) :: diff(size(fsrc)), diffx

        fsrc_diff = 0.
        flux_diff = 0.

        nnod      = size(flux, dim=1)
        ng        = size(flux, dim=2)

        do n= 1, nnod
            diff(n) = 0.
            if (fsrc(n) > 0.) diff(n) = abs(fsrc(n) - fsrc_last(n)) / fsrc(n)
        end do

        fsrc_diff = maxval(diff)

        do g = 1, ng
            do n = 1, nnod
                if (flux(n,g) > 0.) diffx = abs(flux(n,g) - flux_last(n,g)) / flux(n,g)
                flux_diff = max(flux_diff, diffx)
            end do
        end do
    
    end subroutine

    !===============================================================================================!
    ! calculate total source
    !===============================================================================================!

    pure subroutine get_total_src(c, g, total_src)

        class(fdm_rect), intent(inout)  :: c
        integer, intent(in)         :: g
        real(dp), intent(out)       :: total_src(:)

        integer                  :: n, h
        real(dp)                 :: tmp_sum(c % nnod, c % ng)

        tmp_sum = 0.0
        do h = 1, g-1
            do concurrent (n = 1:c % nnod)
                tmp_sum(n,h) = c % xs % sigs(n, h, g) * c % flux(n, h)
            end do
        end do

        do h = g+1, c % ng
            do concurrent (n = 1:c % nnod)
                tmp_sum(n,h) = c % xs % sigs(n, h, g) * c % flux(n, h)
            end do
        end do

        c % scat(:,g) = sum(tmp_sum, dim=2) 

        do concurrent (n = 1:c % nnod)
            total_src(n) = c % xs % chi(n,g) * c % fsrc(n) / c % Keff &
            + c % scat(n,g) + c % exsrc(n,g)
        end do
    
    end subroutine

    !===============================================================================================!
    ! calculate total source
    !===============================================================================================!

    pure subroutine get_total_transient(c, g, total_src)

        class(fdm_rect), intent(inout)  :: c
        integer, intent(in)         :: g
        real(dp), intent(out)       :: total_src(:)

        integer                  :: n, h
        real(dp)                 :: tmp_sum(c % nnod, c % ng)

        tmp_sum = 0.0
        do h = 1, g-1
            do concurrent (n = 1:c % nnod)
                tmp_sum(n,h) = c % xs % sigs(n, h, g) * c % flux(n, h)
            end do
        end do

        do h = g+1, c % ng
            do concurrent (n = 1:c % nnod)
                tmp_sum(n,h) = c % xs % sigs(n, h, g) * c % flux(n, h)
            end do
        end do

        c % scat(:,g) = sum(tmp_sum, dim=2) 

        do concurrent (n = 1:c % nnod)
            total_src(n) = (1. - c % total_beta(c % ind_mat(n)) + c % dfis(n)) * c % xs % chi(n,g) * c % fsrc(n) &
            + c % scat(n,g) + c % exsrc(n,g)
        end do
    
    end subroutine

    !===============================================================================================!
    ! calculate total source (adjoint)
    !===============================================================================================!

    pure subroutine get_total_adjoint(c, g, total_src)

        class(fdm_rect), intent(inout)  :: c
        integer, intent(in)         :: g
        real(dp), intent(out)       :: total_src(:)

        integer                  :: n, h
        real(dp)                 :: tmp_sum(c % nnod, c % ng)

        tmp_sum = 0.0
        do h = 1, g-1
            do concurrent (n = 1:c % nnod)
                tmp_sum(n,h) = c % xs % sigs(n,g,h) * c % flux(n, h)
            end do
        end do

        do h = g+1, c % ng
            do concurrent (n = 1:c % nnod)
                tmp_sum(n,h) = c % xs % sigs(n,g,h) * c % flux(n, h)
            end do
        end do

        c % scat(:,g) = sum(tmp_sum, dim=2)

        do concurrent (n = 1:c % nnod)
            total_src(n) = c % xs % nuf(n,g) * c % fsrc(n) / c % Keff &
            + c % scat(n,g) + c % exsrc(n,g)
        end do
    
    end subroutine

    !===============================================================================================!
    ! calculate fission source
    !===============================================================================================!

    pure subroutine get_fission_src(c)

        class(fdm_rect), intent(inout) :: c

        integer                 :: g, n
        real(dp)                :: tmp_sum(c % nnod, c % ng)

        tmp_sum = 0.0
        do g = 1, c % ng
            do concurrent (n = 1:c % nnod)
                tmp_sum(n,g) = c % flux(n,g) * c % xs % nuf(n,g)
            end do
        end do
        
        c % fsrc = sum(tmp_sum, dim=2)
    
    end subroutine

    !===============================================================================================!
    ! calculate fission source (adjoint)
    !===============================================================================================!

    pure subroutine get_fission_adjoint(c)

        class(fdm_rect), intent(inout) :: c

        integer                 :: g, n
        real(dp)                :: tmp_sum(c % nnod, c % ng)

        tmp_sum = 0.0
        do g = 1, c % ng
            do concurrent (n = 1:c % nnod)
                tmp_sum(n,g) = c % flux(n,g) * c % xs % chi(n,g)
            end do
        end do
        
        c % fsrc = sum(tmp_sum, dim=2)
    
    end subroutine

    !===============================================================================================!
    ! calculate fission source and its integrated value                                             !
    !===============================================================================================!

    subroutine fsrc_extrp(e1, e2, erro, errn, fsrc, print_iter)

        real(dp), intent(in)                  :: e1, e2
        real(dp), dimension(:), intent(in)    :: erro, errn
        real(dp), dimension(:), intent(inout) :: fsrc
        logical, intent(in)                   :: print_iter

        real(dp) :: domiR, mval

        domiR = e2 / e1            ! Dominance ratio
        mval = maxval(abs(erro))
        if (mval * mval < 0.0) domiR = -domiR
        fsrc = fsrc + domiR / (1._dp - domiR) * errn

        if (print_iter) then
            write(fid,*) '    ...FISSION SOURCE EXTRAPOLATED...'
            write(output_unit,*)   '    ...FISSION SOURCE EXTRAPOLATED...'
        end if
    
    end subroutine

    !===============================================================================================!
    ! do volume integration on variable a
    !===============================================================================================!

    pure function integrate(c, a) result(res)

        class(fdm_rect), intent(in) :: c
        real(dp), intent(in)    :: a(:)
        real(dp)                :: res

        integer                 :: n

        res = 0.
        do n = 1, c % nnod
            res = res + a(n) * c % vdel(n)
        end do
    
    end function

    !===============================================================================================!
    ! Set node data and setup the matrix index for linear solver                                          !
    !===============================================================================================!

    subroutine setup_rectangular(c)

        class(fdm_rect)  :: c

        integer :: nnod
        integer :: nxx, nyy, nzz
        integer :: ng

        nnod = c % nnod
        nxx  = c % nxx
        nyy  = c % nyy
        nzz  = c % nzz
        ng   = c % ng

        ! check number material
        if (maxval(c % ind_mat) > c % nmat) call fatal_error(fid, 'Material number ' &
        // n2c(maxval(c % ind_mat)) // ' is greater than number of material')

        allocate(c % df(nnod, ng, 6), c % dn(nnod, ng, 6))
        c % dn = 0.0      ! Init higher-order coupling coef

        ! Setup matrix
        call set_fdm_matrix(c)

        ! set nodal update data
        if (c % kernel .ne. ' FDM') then
            allocate(c % d)
            
            call set_nodal_data(c % d, fid, c % ng, c % nnod, c % nxx, c % nyy, c % nzz, &
            xeast, xwest, ynorth, ysouth, ztop, zbott, .false.)

            call set_nodal_pointer(c % d, c % Keff, c % xs, c % flux, c % exsrc, &
            c % xdel, c % ydel, c % zdel, c % xstag, c % ystag, c % ix, c % iy, c % iz, &
            c % xyz, c % df, c % dn)
            
            call alloc_data(c % d)
        end if

    end subroutine

    !===============================================================================================!
    ! Set FDM Matrix                                                                           !
    !===============================================================================================!

    subroutine set_fdm_matrix(c)

        class(fdm_rect)          :: c

        integer  :: n, i, j, k
        integer  :: n_non_zero

        ! Count non zero elements
        n_non_zero = 0
        do n = 1, c % nnod

            i = c % ix(n); j = c % iy(n); k = c % iz(n)         ! Set i, j, k

            if (k /= 1)                   n_non_zero = n_non_zero + 1
            if (j /= c % xstag(i) % smin) n_non_zero = n_non_zero + 1
            if (i /= c % ystag(j) % smin) n_non_zero = n_non_zero + 1
            n_non_zero = n_non_zero + 1
            if (i /= c % ystag(j) % smax) n_non_zero = n_non_zero + 1
            if (j /= c % xstag(i) % smax) n_non_zero = n_non_zero + 1
            if (k /= c % nzz)             n_non_zero = n_non_zero + 1
        end do

        call allocate_matrix(c % m, c % ng, c % nnod, n_non_zero)
        call set_fdm_index(c)

    end subroutine

    !===============================================================================================!
    ! Set FDM Matrix Index (indexed in CSR)                                                         !
    !===============================================================================================!

    subroutine set_fdm_index(c)

        class(fdm_rect)         :: c

        integer  :: nodp(c % nxx, c % nyy)  !radial node position
        integer  :: n, np, idx
        integer  :: i, j, k

        integer :: nnod
        integer :: nxx, nyy, nzz

        nnod = c % nnod
        nxx  = c % nxx
        nyy  = c % nyy
        nzz  = c % nzz
    
        ! setup radial node position nodp
        nodp = 0
        idx = 0
        do j = 1, nyy
          do i = c % ystag(j) % smin, c % ystag(j) % smax
            idx = idx + 1
            nodp(i,j) = idx
          end do
        end do
    
        ! Calculate number of nodes for one planar
        np = idx
    
        idx = 1
        do n = 1, nnod

            ! Set i, j, k
            i = c % ix(n); j = c % iy(n); k = c % iz(n)
    
            c % m % ind % row(n) = idx
      
            ! Lower diagonal matrix element for z-direction
            if (k /= 1) then
                c % m % ind % col(idx) = n - np
                idx                    = idx + 1
            end if
      
            ! Lower diagonal matrix element for y-direction
            if (j /= c % xstag(i) % smin) then
                c % m % ind % col(idx) = n - (nodp(i,j) - nodp(i,j-1))
                idx                    = idx + 1
            end if
      
            ! Lower diagonal matrix element for x-direction
            if (i /= c % ystag(j) % smin) then
                c % m % ind % col(idx) = n - 1
                idx                    = idx + 1
            end if
      
            ! Diagonal matrix elementss
            c % m % ind % col(idx) = n
            c % m % ind % diag(n)  = idx
            idx                    = idx + 1
      
             ! Upper diagonal matrix element for x-direction
            if (i /= c % ystag(j) % smax) then
                c % m % ind % col(idx) = n + 1
                idx                    = idx + 1
            end if
      
            ! Upper diagonal matrix element for y-direction
            if (j /= c % xstag(i) % smax) then
                c % m % ind % col(idx) = n + (nodp(i,j+1) - nodp(i,j))
                idx                    = idx + 1
            end if
      
            ! Upper diagonal matrix element for z-direction
            if (k /= nzz) then
                c % m % ind % col(idx) = n + np
                idx                    = idx + 1
            end if
    
        end do

        c % m % ind % row(nnod+1) = idx

    end subroutine

    !===============================================================================================!
    ! Setup sparse penta-diagonal matrix indexed in CSR                                             !
    !===============================================================================================!

    subroutine build_matrix(c, upd_fdm_coupling, sigr)

        class(fdm_rect), target       :: c
        logical, intent(in)            :: upd_fdm_coupling
        real(dp), intent(in), optional :: sigr(:,:)

        integer                     :: n, g, idx

        real(dp) :: sigma_removal(c % nnod, c % ng)
        integer  :: nnod, ng
        integer  :: i, j, k

        nnod = c % nnod
        ng   = c % ng
    
        ! If need to calculate FDM coupling coefficients
        if (upd_fdm_coupling) call fdm_coupling_coef(c)

        if (present(sigr)) then
            sigma_removal = sigr
        else
            sigma_removal = c % xs % sigr
        end if
    
        ! Setup CMFD linear system
        do g = 1, ng
            idx = 0
            do n = 1, nnod

                ! Set i, j, k
                i = c % ix(n); j = c % iy(n); k = c % iz(n)
      
                ! Lower diagonal matrix element for z-direction
                if (k /= 1) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(c % df(n,g,6) - c % dn(n,g,6)) / c % zdel(k)
                end if
        
                ! Lower diagonal matrix element for y-direction
                if (j /= c % xstag(i) % smin) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(c % df(n,g,4) - c % dn(n,g,4)) / c % ydel(j)
                end if
        
                ! Lower diagonal matrix element for x-direction
                if (i /= c % ystag(j) % smin) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(c % df(n,g,2) - c % dn(n,g,2)) / c % xdel(i)
                end if
        
                ! Diagonal matrix elementss
                idx = idx + 1
                c % m % mtx(g) % elmn(idx) = &
                (c % df(n,g,1) + c % df(n,g,2) - c % dn(n,g,1) + c % dn(n,g,2)) / c % xdel(i) + &
                (c % df(n,g,3) + c % df(n,g,4) - c % dn(n,g,3) + c % dn(n,g,4)) / c % ydel(j) + &
                (c % df(n,g,5) + c % df(n,g,6) - c % dn(n,g,5) + c % dn(n,g,6)) / c % zdel(k) + &
                sigma_removal(n,g)
        
                 ! Upper diagonal matrix element for x-direction
                if (i /= c % ystag(j) % smax) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(c % df(n,g,1) + c % dn(n,g,1)) / c % xdel(i)
                end if
        
                ! Upper diagonal matrix element for y-direction
                if (j /= c % xstag(i) % smax) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(c % df(n,g,3) + c % dn(n,g,3)) / c % ydel(j)
                end if
        
                ! Upper diagonal matrix element for z-direction
                if (k /= c % nzz) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(c % df(n,g,5) + c % dn(n,g,5)) / c % zdel(k)
                end if
      
            end do
        end do

    end subroutine

    !===============================================================================================!
    ! calculate FDM nodal coupling coefficients                                                     !
    !===============================================================================================!

    pure subroutine fdm_coupling_coef(c)

        class(fdm_rect), intent(inout) :: c

        real(dp), parameter :: big = 1.e30
        integer  :: g
        real(dp) :: DN, D
        real(dp) :: xd, yd, zd
        real(dp) :: xdn, ydn, zdn
        integer  :: i, j, k
        integer  :: n

        do g = 1, c % ng
            do n = 1, c % nnod

                ! Set i, j, k
                i = c % ix(n); j = c % iy(n); k = c % iz(n)
    
                D  = c % xs % D(n,g)
                xd = c % xdel(i)
                yd = c % ydel(j)
                zd = c % zdel(k)
    
                ! Set FDM coupling coefficients in x direction
                if (i == c % ystag(j) % smax) then
                    if (xeast == zero_flux) then
                        c % df(n,g,1) = 2. * big * D / (2. * D + big * xd)
                    elseif (xeast == zero_incoming) then
                        c % df(n,g,1) =  D / (2. * D + 0.5 * xd)
                    else
                        c % df(n,g,1) =  0.
                    end if
                else
                    DN  = c % xs % D(c % xyz(i+1, j, k), g)
                    xdn = c % xdel(i+1)
                    c % df(n,g,1) =  2. * D * DN / (D * xdn +  DN * xd)
                endif
    
                if (i == c % ystag(j) % smin) then
                    if (xwest == zero_flux) then
                        c % df(n,g,2) = 2. * big * D / (2. * D + big * xd)
                    elseif (xwest == zero_incoming) then
                        c % df(n,g,2) =  D / (2. * D + 0.5 * xd)
                    else
                        c % df(n,g,2) =  0.
                    end if
                else
                    DN  = c % xs % D(c % xyz(i-1, j, k), g)
                    xdn = c % xdel(i-1)
                    c % df(n,g,2) =  2. * D * DN / (D * xdn +  DN * xd)
                endif
    
                ! Set FDM coupling coefficients in y direction
                if (j == c % xstag(i) % smax) then
                    if (ynorth == zero_flux) then
                        c % df(n,g,3) = 2. * big * D / (2. * D + big * yd)
                    elseif (ynorth == zero_incoming) then
                        c % df(n,g,3) =  D / (2. * D + 0.5 * yd)
                    else
                        c % df(n,g,3) =  0.
                    end if
                else
                    DN  = c % xs % D(c % xyz(i, j+1, k), g)
                    ydn = c % ydel(j+1)
                    c % df(n,g,3) =  2. * D * DN / (D * ydn +  DN * yd)
                endif
    
                if (j == c % xstag(i) % smin) then
                    if (ysouth == 0) then
                        c % df(n,g,4) = 2. * big * D / (2. * D + big * yd)
                    elseif (ysouth == 1) then
                        c % df(n,g,4) =  D / (2. * D + 0.5 * yd)
                    else
                        c % df(n,g,4) =  0.
                    end if
                else
                    DN  = c % xs % D(c % xyz(i, j-1, k), g)
                    ydn = c % ydel(j-1)
                    c % df(n,g,4) =  2. * D * DN / (D * ydn +  DN * yd)
                endif
    
                ! Set FDM coupling coefficients in z direction
                if (k == c % nzz) then
                    if (ztop == 0) then
                        c % df(n,g,5) = 2. * big * D / (2. * D + big * zd)
                    elseif (ztop == 1) then
                        c % df(n,g,5) =  D / (2. * D + 0.5 * zd)
                    else
                        c % df(n,g,5) =  0.
                    end if
                else
                    DN  = c % xs % D(c % xyz(i, j, k+1), g)
                    zdn = c % zdel(k+1)
                    c % df(n,g,5) =  2. * D * DN / (D * zdn +  DN * zd)
                endif
        
                if (k == 1) then
                    if (zbott == 0) then
                        c % df(n,g,6) = 2. * big * D / (2. * D + big * zd)
                    elseif (zbott == 1) then
                        c % df(n,g,6) =  D / (2. * D + 0.5 * zd)
                    else
                        c % df(n,g,6) =  0.
                    end if
                else
                    DN  = c % xs % D(c % xyz(i, j, k-1), g)
                    zdn = c % zdel(k-1)
                    c % df(n,g,6) =  2. * D * DN / (D * zdn +  DN * zd)
                endif
            end do
        end do

    end subroutine
    
end module
    