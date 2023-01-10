module fdm

    use iso_fortran_env, only: real64
    use time
    use xsec
    use stagger
    use linear_system
    use util
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

    ! core data type
    type :: core_base
        character(4)        :: kernel        ! Nodal kernel used
        integer             :: nnod          ! Number of nodes
        real(dp), pointer   :: vdel(:)       ! Delta nodes' volume in cm3
    
        integer             :: ng            ! number of neutron energy groups

        integer             :: nmat          ! number of materials
    
        real(dp)            :: Keff = 1.0    ! Multiplication factor

        real(dp), pointer   :: fsrc(:)       ! fission source
        real(dp), pointer   :: flux(:,:)     ! flux
        real(dp), pointer   :: exsrc(:,:)  ! external source
    
        type(xs_rect)       :: xs ! node-wise xs
        type(matrix)        :: m  ! FDM matrix for linear system solver

        logical             :: is_first = .true.
    end type

    ! Rectangilar core data type
    type, extends(core_base), public :: core_rect
        real(dp), pointer :: xdel(:)        ! Delta x in cm
        real(dp), pointer :: ydel(:)        ! Delta y in cm
        real(dp), pointer :: zdel(:)        ! Delta z in cm
    
        integer  :: nxx                     ! Number of nodes in x direction
        integer  :: nyy                     ! Number of nodes in y direction
        integer  :: nzz                     ! Number of nodes z direction

        integer, pointer :: ix(:)           ! index in x direction for given node index
        integer, pointer :: iy(:)           ! index in y direction for given node index
        integer, pointer :: iz(:)           ! index in z direction for given node index
        integer, pointer :: xyz(:,:,:)      ! node index for given index in x, y, and z directions

        integer, pointer :: mat_map(:,:,:)  ! material map
        
        ! coupling coeff (order => east, west, north, south, top, bottom)
        real(dp), allocatable      :: df(:,:,:)   ! D tilde : FDM nodal coupling coefficients
        real(dp), allocatable      :: dn(:,:,:)   ! D hat   : Corrected (higher order) nodal coupling coefficients

        type(staggered), pointer   :: xstag(:)    ! Staggered mesh data
        type(staggered), pointer   :: ystag(:)    ! Staggered mesh data

        type(nodal_type), allocatable :: d        ! nodal update module

    end type

    public :: set_rect_data, set_rect_pointer, outer_iter

    contains

    !===============================================================================================!
    ! Set FDM data for rectangular geometry
    !===============================================================================================!

    subroutine set_rect_data(c, kernel, ng, nnod, nxx, nyy, nzz, &
        nmat, east, west, north, south, bottom, top, output_fid)

        class(core_rect)              :: c
        character(*), intent(in)      :: kernel
        integer, intent(in)           :: nnod                                     ! total number of nodes
        integer, intent(in)           :: ng                                       ! number of group
        integer, intent(in)           :: nxx, nyy, nzz                            ! number of nodes in x, y, and z directions
        integer, intent(in)           :: nmat                                     ! number of material
        integer, intent(in), optional :: east, west, north, south, bottom, top    ! BCs
        integer, intent(in), optional :: output_fid                               ! output file unit number

        !set output file unit number
        if(present(output_fid)) fid = output_fid

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
        if (present(east))   xeast  = east
        if (present(west))   xwest  = west
        if (present(north))  ynorth = north
        if (present(south))  ysouth = south
        if (present(bottom)) zbott  = bottom
        if (present(top))    ztop   = top
        
    end subroutine

    !===============================================================================================!
    ! Set FDM pointers for rectangular geometry
    !===============================================================================================!

    subroutine set_rect_pointer(c, flux, fsrc, exsrc, mat_map, xdel, ydel, zdel, vdel, &
        xstag, ystag, ix, iy, iz, xyz)

        class(core_rect)       :: c
        real(dp), intent(inout), target        :: flux(:,:)
        real(dp), intent(in), target           :: fsrc(:), exsrc(:,:)
        integer, intent(in), target            :: mat_map(:,:,:)      ! 3D Material map
        real(dp), intent(in), target           :: xdel(:), ydel(:), zdel(:), vdel(:)
        integer, intent(in), target            :: ix(:), iy(:), iz(:), xyz(:,:,:)
        type(staggered), intent(in), target    :: xstag(:), ystag(:)

        flux = 1.0  ! Init flux
        c % flux => flux
        if (size(c % flux, dim=1) .ne. c % nnod) stop 'Flux first-dimension size does not match'
        if (size(c % flux, dim=2) .ne. c % ng)   stop 'Flux second-dimension size does not match'

        c % fsrc => fsrc
        if (size(c % fsrc) .ne. c % nnod) stop 'Fission source size does not match'

        c % exsrc => exsrc
        if (size(c % exsrc, dim=1) .ne. c % nnod) stop 'exsrc first-dimension size does not match'
        if (size(c % exsrc, dim=2) .ne. c % ng)   stop 'exsrc second-dimension size does not match'

        c % mat_map => mat_map
        if (size(c % mat_map, dim=1) .ne. c % nxx) stop 'mat_map first-dimension size does not match'
        if (size(c % mat_map, dim=2) .ne. c % nyy) stop 'mat_map second-dimension size does not match'
        if (size(c % mat_map, dim=3) .ne. c % nzz) stop 'mat_map third-dimension size does not match'

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

    end subroutine

    !===============================================================================================!
    ! to perfom outer iteration                                                                     !
    !===============================================================================================!

    subroutine outer_iter(c, max_flux_err, max_fiss_err, &
        max_outer, max_inner, extrp, nod_intv_inp, print_iter)

        class(core_rect)                   :: c
        real(dp), intent(in)               :: max_flux_err, max_fiss_err
        integer, intent(in)                :: max_outer, max_inner
        logical, intent(in)                :: print_iter
        integer, intent(in)                :: extrp          ! extrapolation interval
        integer, intent(in), optional      :: nod_intv_inp    ! nodal coupling coefficients update interval

        integer  :: i, g
        real(dp) :: fc, fco                                  ! new and old integrated fission SOURCES
        real(dp) :: Keo
        real(dp) :: flux_prev(c % nnod, c % ng)
        real(dp) :: fsrc_prev(c % nnod)
        real(dp) :: bs(c % nnod)
        real(dp) :: flux_diff, fsrc_diff
        real(dp) :: e1, e2
        real(dp) :: errn(c % nnod), erro(c % nnod)            ! current and past error vectors
        integer  :: nod_intv
        logical  :: converge

        call fdm_time % on

        if (c % is_first) then
            call setup_rectangular(c)
            c % is_first = .false.
        end if

        if (present(nod_intv_inp)) then
            nod_intv = nod_intv_inp
        else
            nod_intv = ceiling((c % nxx + c % nyy + c % nzz) / 2.5)
        end if
    
        ! Initialize fission source
        call get_fission_src(c)
        fc = integrate(c, c % fsrc)

        ! Init error vector
        errn = 1._DP
        e1 = integrate(c, errn)

        call fdm_time % off
    
        !Start outer iteration
        converge = .false.
        do i = 1, max_outer

            call fdm_time % on
            fco        = fc               ! Save old integrated fission source
            Keo        = c % Keff         ! Save old keff
            flux_prev  = c % flux         ! Save old flux
            fsrc_prev  = c % fsrc         ! Save old fission source
            erro       = errn             ! Save old fission source error/difference
            
            ! Inner iteration
            do g = 1, c % ng
                call get_total_src(c, g, bs)
                call linear_solver(c % m % mtx(g), c % m % ind, max_inner, bs, c % flux(:,g))
            end do

            call get_fission_src(c)                      !Update fission source
            errn = c % fsrc - fsrc_prev
            e2 = norm2(errn)
            if (MOD(i, extrp) == 0) call fsrc_extrp(e1, e2, erro, errn, c % fsrc, print_iter) ! Fission source extrapolation
            e1 = e2
            fc = integrate(c, c % fsrc)
            c % Keff = Keo * fc / fco             ! Update Keff
            call get_difference(c % flux, flux_prev, c % fsrc, fsrc_prev, flux_diff, fsrc_diff)
            if (print_iter) then
                write(fid,'(I5,F13.6,2ES15.5)') i, c % Keff, fsrc_diff, flux_diff
                write(*,'(I5,F13.6,2ES15.5)')   i, c % Keff, fsrc_diff, flux_diff
            end if
            if ((max_flux_err > flux_diff) .AND. (max_fiss_err > fsrc_diff)) then
                converge = .true.
                exit
            end if
            call fdm_time % off

            if (mod(i, nod_intv) == 0) call nodal_update(c, .true.)  ! Nodal coefficients update
        end do
    
        if (.not. converge) THEN
          write(*,*)
          write(*,*) '  MAXIMUM NUMBER OF OUTER ITERATION IS REACHED IN FORWARD CALCULATION.'
          write(*,*) '  CHECK PROBLEM SPECIFICATION OR CHANGE ITERATION CONTROL (%ITER).'
          write(*,*) '  PERHAPS BY MAKING FISSION SOURCE INTERPOLATION MORE FREQUENT'
          write(*,*) '  KOMODO IS STOPING...'
          stop
        end if
    
    end subroutine

    !===============================================================================================!
    ! To update nodal coupling coefficients (D hat)
    !===============================================================================================!

    subroutine nodal_update(c, print_iter)
      
        class(core_rect), intent(in)  :: c
        logical, intent(in)           :: print_iter

        call nodal_time % on
      
        ! Update nodal coupling coefficients
        if (c % kernel == ' FDM') then
            return
        elseif (c % kernel == 'SANM') then
            call nodal_update_sanm(c % d)
        else
            call nodal_update_pnm(c % d)
        end if

        ! Update CMFD matrix
        call setup_matrix(c, upd_fdm_coupling=.false.)

        if (print_iter) then
          write(fid,*) '    .....NODAL COUPLING UPDATED..... '
          write(*,*) '    .....NODAL COUPLING UPDATED..... '
        end if

        call nodal_time % off
      
      end subroutine

    !===============================================================================================!
    ! Search maximum point wise fission source Relative difference, and
    ! search maximum point wise flux relative difference
    !===============================================================================================!

    pure subroutine get_difference(flux, flux_last, fsrc, fsrc_last, flux_diff, fsrc_diff)

        real(dp), intent(in)  :: flux(:,:), flux_last(:,:)
        real(dp), intent(in)  :: fsrc(:), fsrc_last(:)
        real(dp), intent(out) :: flux_diff, fsrc_diff

        integer  :: nnod, ng
        integer  :: n, g
        real(dp) :: diff

        fsrc_diff = 0.
        flux_diff = 0.

        nnod      = size(flux, dim=1)
        ng        = size(flux, dim=2)

        do n= 1, nnod
            diff = abs(fsrc(n) - fsrc_last(n)) / max(fsrc(n), 1.e-30_dp)
            fsrc_diff = max(fsrc_diff, diff)
        end do

        do g = 1, ng
            do n = 1, nnod
                diff = abs(flux(n, g) - flux_last(n, g)) / max(flux(n, g), 1.e-30_dp)
                flux_diff = max(flux_diff, diff)
            end do
        end do
    
    end subroutine

    !===============================================================================================!
    ! calculate total source                                                     !
    !===============================================================================================!

    pure subroutine get_total_src(c, g, total_src)

        class(core_rect), intent(inout)  :: c
        integer, intent(in)         :: g
        real(dp), intent(out)       :: total_src(:)

        integer                  :: n, h
        real(dp)                 :: scat_src(c % nnod, c % ng)
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

        scat_src(:,g) = sum(tmp_sum, dim=2) 

        do concurrent (n = 1:c % nnod)
            total_src(n) = c % xs % chi(n, g) * c % fsrc(n) / c % Keff &
            + scat_src(n, g) + c % exsrc(n, g)
        end do
    
    end subroutine

    !===============================================================================================!
    ! calculate fission source
    !===============================================================================================!

    pure subroutine get_fission_src(c)

        class(core_rect), intent(inout) :: c

        integer                 :: g, n
        real(dp)                :: tmp_sum(c % nnod, c % ng)

        tmp_sum = 0.0
        do g = 1, c % ng
            do concurrent (n = 1:c % nnod)
                tmp_sum(n,g) = c % flux(n, g) * c % xs % nuf(n, g)
            end do
        end do
        
        c % fsrc = sum(tmp_sum, dim=2)
    
    end subroutine

    !===============================================================================================!
    ! calculate fission source (adjoint)
    !===============================================================================================!

    pure subroutine get_fission_adjoint(c)

        class(core_rect), intent(inout) :: c

        integer                 :: g, n
        real(dp)                :: tmp_sum(c % nnod, c % ng)

        tmp_sum = 0.0
        do g = 1, c % ng
            do concurrent (n = 1:c % nnod)
                tmp_sum(n,g) = c % flux(n, g) * c % xs % chi(n, g)
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
        fsrc = fsrc + domiR / (1._DP - domiR) * errn

        if (print_iter) then
            write(fid,*) '    ...FISSION SOURCE EXTRAPOLATED...'
            write(*,*)   '    ...FISSION SOURCE EXTRAPOLATED...'
        end if
    
    end subroutine

    !===============================================================================================!
    ! do volume integration on variable a
    !===============================================================================================!

    pure function integrate(c, a) result(res)

        class(core_rect), intent(in) :: c
        real(dp), intent(in)    :: a(:)
        real(dp)                :: res

        integer                 :: n

        res = 0.
        do n = 1, c % nnod
            res = res + a(n) * c % vdel(n)
        end do
    
    end function

    !===============================================================================================!
    ! Set node data and setup the matrix for linear solver                                          !
    !===============================================================================================!

    subroutine setup_rectangular(c)

        class(core_rect)  :: c

        integer :: nnod
        integer :: nxx, nyy, nzz
        integer :: ng

        nnod = c % nnod
        nxx  = c % nxx
        nyy  = c % nyy
        nzz  = c % nzz
        ng   = c % ng

        ! check number material
        if (maxval(c % mat_map) > c % nmat) call fatal_error(fid, 'Material number ' &
        // n2c(maxval(c % mat_map)) // ' is greater than number of material')

        allocate(c % df(nnod, ng, 6), c % dn(nnod, ng, 6))
        c % dn = 0.0      ! Init higher-order coupling coef

        ! Setup matrix
        call set_fdm_matrix(c)

        if (c % kernel .ne. ' FDM') then
            allocate(c % d)
            
            call set_nodal_data(c % d, c % ng, c % nnod, c % nxx, c % nyy, c % nzz, &
            xeast, xwest, ynorth, ysouth, ztop, zbott, 'static', fid, .false.)

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

        class(core_rect)          :: c

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
        call setup_matrix(c, upd_fdm_coupling=.true.)

    end subroutine

    !===============================================================================================!
    ! Set FDM Matrix Index (indexed in CSR)                                                         !
    !===============================================================================================!

    subroutine set_fdm_index(c)

        class(core_rect)         :: c

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

    subroutine setup_matrix(c, upd_fdm_coupling)

        class(core_rect), target   :: c
        logical, intent(in)   :: upd_fdm_coupling

        integer                     :: n, g, idx

        integer :: nnod, ng
        integer :: i, j, k

        nnod = c % nnod
        ng   = c % ng
    
        ! If need to calculate FDM coupling coefficients
        if (upd_fdm_coupling) call fdm_coupling_coef(c)
    
        ! Setup CMFD linear system
        do g = 1, ng
            idx = 0
            do n = 1, nnod

                ! Set i, j, k
                i = c % ix(n); j = c % iy(n); k = c % iz(n)
      
                ! Lower diagonal matrix element for z-direction
                if (k /= 1) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(c % df(n, g, 6) - c % dn(n, g, 6)) / c % zdel(k)
                end if
        
                ! Lower diagonal matrix element for y-direction
                if (j /= c % xstag(i) % smin) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(c % df(n, g, 4) - c % dn(n, g, 4)) / c % ydel(j)
                end if
        
                ! Lower diagonal matrix element for x-direction
                if (i /= c % ystag(j) % smin) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(c % df(n, g, 2) - c % dn(n, g, 2)) / c % xdel(i)
                end if
        
                ! Diagonal matrix elementss
                idx = idx + 1
                c % m % mtx(g) % elmn(idx) = &
                (c % df(n, g, 1) + c % df(n, g, 2) - c % dn(n, g, 1) + c % dn(n, g, 2)) / c % xdel(i) + &
                (c % df(n, g, 3) + c % df(n, g, 4) - c % dn(n, g, 3) + c % dn(n, g, 4)) / c % ydel(j) + &
                (c % df(n, g, 5) + c % df(n, g, 6) - c % dn(n, g, 5) + c % dn(n, g, 6)) / c % zdel(k) + &
                c % xs % sigr(n, g)
        
                 ! Upper diagonal matrix element for x-direction
                if (i /= c % ystag(j) % smax) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(c % df(n, g, 1) + c % dn(n, g, 1)) / c % xdel(i)
                end if
        
                ! Upper diagonal matrix element for y-direction
                if (j /= c % xstag(i) % smax) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(c % df(n, g, 3) + c % dn(n, g, 3)) / c % ydel(j)
                end if
        
                ! Upper diagonal matrix element for z-direction
                if (k /= c % nzz) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(c % df(n, g, 5) + c % dn(n, g, 5)) / c % zdel(k)
                end if
      
            end do
        end do

    end subroutine

    !===============================================================================================!
    ! calculate FDM nodal coupling coefficients                                                     !
    !===============================================================================================!

    pure subroutine fdm_coupling_coef(c)

        class(core_rect), intent(inout) :: c

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
    
                D  = c % xs % D(n, g)
                xd = c % xdel(i)
                yd = c % ydel(j)
                zd = c % zdel(k)
    
                ! Set FDM coupling coefficients in x direction
                if (i == c % ystag(j) % smax) then
                    if (xeast == zero_flux) then
                        c % df(n, g, 1) = 2. * big * D / (2. * D + big * xd)
                    elseif (xeast == zero_incoming) then
                        c % df(n, g, 1) =  D / (2. * D + 0.5 * xd)
                    else
                        c % df(n, g, 1) =  0.
                    end if
                else
                    DN  = c % xs % D(c % xyz(i+1, j, k), g)
                    xdn = c % xdel(i+1)
                    c % df(n, g, 1) =  2. * D * DN / (D * xdn +  DN * xd)
                endif
    
                if (i == c % ystag(j) % smin) then
                    if (xwest == zero_flux) then
                        c % df(n, g, 2) = 2. * big * D / (2. * D + big * xd)
                    elseif (xwest == zero_incoming) then
                        c % df(n, g, 2) =  D / (2. * D + 0.5 * xd)
                    else
                        c % df(n, g, 2) =  0.
                    end if
                else
                    DN  = c % xs % D(c % xyz(i-1, j, k), g)
                    xdn = c % xdel(i-1)
                    c % df(n, g, 2) =  2. * D * DN / (D * xdn +  DN * xd)
                endif
    
                ! Set FDM coupling coefficients in y direction
                if (j == c % xstag(i) % smax) then
                    if (ynorth == zero_flux) then
                        c % df(n, g, 3) = 2. * big * D / (2. * D + big * yd)
                    elseif (ynorth == zero_incoming) then
                        c % df(n, g, 3) =  D / (2. * D + 0.5 * yd)
                    else
                        c % df(n, g, 3) =  0.
                    end if
                else
                    DN  = c % xs % D(c % xyz(i, j+1, k), g)
                    ydn = c % ydel(j+1)
                    c % df(n, g, 3) =  2. * D * DN / (D * ydn +  DN * yd)
                endif
    
                if (j == c % xstag(i) % smin) then
                    if (ysouth == 0) then
                        c % df(n, g, 4) = 2. * big * D / (2. * D + big * yd)
                    elseif (ysouth == 1) then
                        c % df(n, g, 4) =  D / (2. * D + 0.5 * yd)
                    else
                        c % df(n, g, 4) =  0.
                    end if
                else
                    DN  = c % xs % D(c % xyz(i, j-1, k), g)
                    ydn = c % ydel(j-1)
                    c % df(n, g, 4) =  2. * D * DN / (D * ydn +  DN * yd)
                endif
    
                ! Set FDM coupling coefficients in z direction
                if (k == c % nzz) then
                    if (ztop == 0) then
                        c % df(n, g, 5) = 2. * big * D / (2. * D + big * zd)
                    elseif (ztop == 1) then
                        c % df(n, g, 5) =  D / (2. * D + 0.5 * zd)
                    else
                        c % df(n, g, 5) =  0.
                    end if
                else
                    DN  = c % xs % D(c % xyz(i, j, k+1), g)
                    zdn = c % zdel(k+1)
                    c % df(n, g, 5) =  2. * D * DN / (D * zdn +  DN * zd)
                endif
        
                if (k == 1) then
                    if (zbott == 0) then
                        c % df(n, g, 6) = 2. * big * D / (2. * D + big * zd)
                    elseif (zbott == 1) then
                        c % df(n, g, 6) =  D / (2. * D + 0.5 * zd)
                    else
                        c % df(n, g, 6) =  0.
                    end if
                else
                    DN  = c % xs % D(c % xyz(i, j, k-1), g)
                    zdn = c % zdel(k-1)
                    c % df(n, g, 6) =  2. * D * DN / (D * zdn +  DN * zd)
                endif
            end do
        end do

    end subroutine




        ! real(dp), allocatable :: fuel_temp         ! Fuel temperature
        ! real(dp), allocatable :: mod_temp          ! Moderator temperature
        ! real(dp), allocatable :: mod_dens          ! Moderator density
        ! real(dp), allocatable :: node_power        ! nodes power (watt)
        ! real(dp), allocatable :: tfm(:)            ! Fuel pin mesh temperature for each nodes
        ! real(dp), allocatable :: ent               ! Coolant Enthalpy (J/Kg)
        ! real(dp), allocatable :: heatf             ! Heat flux (W/m2
        ! real(dp), allocatable :: frate             ! coolant mass flow rate

        ! real(dp), allocatable :: prec_dens(:)      ! precusor_density
        ! real(dp), allocatable :: trans_src(:)      ! source for transient mode
        ! real(dp), allocatable :: exsrc(:)        ! external source
        ! real(dp), allocatable :: omega(:)          ! Exponential transformation constant
        ! real(dp), allocatable :: sigrp(:)          ! Initial removal cross sections before added by parameters required for transient
        ! real(dp), allocatable :: L    (:)          ! Total leakages for node n and group g
        ! real(dp), allocatable :: dfis               

        ! integer               :: i_mat                ! material index

    
end module
    