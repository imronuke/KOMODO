module fdm

    use iso_fortran_env, only: real64
    use xsec
    use node
    use linear_system
    use util
    use nodal

    implicit none

    private

    save

    integer, parameter :: dp = real64

    integer, parameter  :: n_rect_surf = 6           ! order is => east, west, north, top, bottom
    real(dp), parameter :: big = 1.e30

    integer, parameter  :: zero_flux = 0
    integer, parameter  :: zero_incoming = 1
    integer, parameter  :: reflective = 2

    integer :: fid                          ! file unit number

    integer :: xwest, xeast                 ! West and East Boundary conditions
    integer :: ysouth, ynorth               ! South and North Boundary conditions
    integer :: zbott, ztop                  ! Bottom and Top Boundary conditions

    ! core data type
    type :: core_base
        character(4)                 :: kernel        ! Nodal kernel used
        integer                      :: nnod          ! Number of nodes
        real(dp), allocatable        :: vdel(:)       ! Delta nodes' volume in cm3
    
        integer                      :: ng            ! number of neutron energy groups
    
        real(dp)                     :: Keff = 1.0    ! Multiplication factor
        real(dp)                     :: int_fsrc      ! intgerated fission source

        real(dp), allocatable         :: fsrc(:)             ! fission source
        real(dp), pointer             :: flux(:,:)           ! flux
        real(dp), allocatable         :: scat_src(:,:)       ! scattering source
        real(dp), allocatable         :: ext_src(:,:)        ! external source
    
        type(xs_rect)                  :: xs          ! node-wise xs
        type(node_rect), allocatable   :: node(:)     ! Node data
        type(matrix)                   :: m           ! FDM matrix for linear system solver
        type(nodal_type), allocatable  :: d

        logical                        :: first = .true.

    end type

    type, extends(core_base), public :: core_rect
        real(dp), allocatable :: xdel(:)        ! Delta x in cm
        real(dp), allocatable :: ydel(:)        ! Delta y in cm
        real(dp), allocatable :: zdel(:)        ! Delta z in cm
    
        integer  :: nxx                          ! Number of nodes in x direction
        integer  :: nyy                          ! Number of nodes in y direction
        integer  :: nzz                          ! Number of nodes z direction

        type(staggered), allocatable   :: xstag(:)    ! Staggered mesh data
        type(staggered), allocatable   :: ystag(:)    ! Staggered mesh data
    end type

    public :: set_fdm_parameters, set_nodes, outer_iter

    contains

    !===============================================================================================!
    ! to perfom outer iteration                                                                     !
    !===============================================================================================!

    subroutine outer_iter(c, max_flux_err, max_fiss_err, flux, &
        max_outer, max_inner, extrp, nod_intv_inp, print_iter)

        class(core_rect)                   :: c
        real(dp), intent(in)               :: max_flux_err, max_fiss_err
        real(dp), intent(inout), target    :: flux(:,:)
        integer, intent(in)                :: max_outer, max_inner
        logical, intent(in)                :: print_iter
        integer, intent(in)                :: extrp          ! extrapolation interval
        integer, intent(in), optional      :: nod_intv_inp    ! nodal coupling coefficients update interval

        integer  :: i, g
        real(dp) :: fc, Keo                                  ! old integrated fission sources and keff
        real(dp) :: flux_prev(c % nnod, c % ng)
        real(dp) :: fsrc_prev(c % nnod)
        real(dp) :: bs(c % nnod)
        real(dp) :: flux_diff, fsrc_diff
        real(dp) :: e1, e2
        real(dp) :: errn(c % nnod), erro(c % nnod)            ! current and past error vectors
        integer  :: nod_intv

        if (present(nod_intv_inp)) then
            nod_intv = nod_intv_inp
        else
            nod_intv = ceiling((c % nxx + c % nyy + c % nzz) / 2.5)
        end if

        if (c % first) then
            flux = 1.0
            c % flux => flux
            c % d % flux => flux
            c % first = .false.
        end if
    
        ! Initialize fission source
        call get_fission_src(c, flux)
        c % int_fsrc = integrate(c, c % fsrc)

        ! Init error vector
        errn = 1._DP
        e1 = integrate(c, errn)
    
        !Start outer iteration
        do i = 1, max_outer
            fc         = c % int_fsrc     ! Save old integrated fission source
            Keo        = c % Keff         ! Save old keff
            flux_prev  = flux             ! Save old flux
            fsrc_prev  = c % fsrc         ! Save old fission source
            erro       = errn             ! Save old fission source error/difference
            
            ! Inner iteration
            do g = 1, c % ng
                call get_total_src(c, g, flux, bs)
                call linear_solver(c % m % mtx(g), c % m % ind, max_inner, bs, flux(:,g))
            end do

            call get_fission_src(c, flux)                      !Update fission source
            errn = c % fsrc - fsrc_prev
            e2 = norm2(errn)
            if (MOD(i, extrp) == 0) call fsrc_extrp(e1, e2, erro, errn, c % fsrc, print_iter) ! Fission source extrapolation
            e1 = e2
            c % int_fsrc = integrate(c, c % fsrc)
            c % Keff = Keo * c % int_fsrc / fc             ! Update Keff
            call get_difference(flux, flux_prev, c % fsrc, fsrc_prev, flux_diff, fsrc_diff)
            if (print_iter) then
                write(fid,'(I5,F13.6,2ES15.5)') i, c % Keff, fsrc_diff, flux_diff
                write(*,'(I5,F13.6,2ES15.5)')   i, c % Keff, fsrc_diff, flux_diff
            end if
            if ((max_flux_err > flux_diff) .AND. (max_fiss_err > fsrc_diff)) exit
            if (mod(i, nod_intv) == 0) call nodal_update(c, .true.)  ! Nodal coefficients update
        end do
    
        if (i-1 == max_outer) THEN
          write(*,*)
          write(*,*) '  MAXIMUM NUMBER OF OUTER ITERATION IS REACHED IN FORWARD CALCULATION.'
          write(*,*) '  CHECK PROBLEM SPECIFICATION OR CHANGE ITERATION CONTROL (%ITER).'
          write(*,*) '  PERHAPS BY MAKING FISSION SOURCE INTERPOLATION MORE FREQUENT'
          write(*,*) '  KOMODO IS STOPING...'
          STOP
        end if
    
    end subroutine

    !===============================================================================================!
    ! To update nodal coupling coefficients (D hat)
    !===============================================================================================!

    subroutine nodal_update(c, print_iter)
      
        class(core_rect), intent(in)  :: c
        logical, intent(in)           :: print_iter
      
        !Update nodal coupling coefficients
        if (c % kernel == ' FDM') then
            return
        elseif (c % kernel == 'SANM') then
            call nodal_update_sanm(c % d)
        else
            call nodal_update_pnm(c % d)
        end if

        !Update CMFD matrix
        call setup_matrix(c, upd_fdm_coupling=.false.)

        if (print_iter) then
          write(fid,*) '    .....NODAL COUPLING UPDATED..... '
          write(*,*) '    .....NODAL COUPLING UPDATED..... '
        end if
      
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

    pure subroutine get_total_src(c, g, flux, total_src)

        class(core_rect), intent(inout)  :: c
        integer, intent(in)         :: g
        real(dp), intent(in)        :: flux(:,:)
        real(dp), intent(out)       :: total_src(:)

        integer                  :: n, h
        real(dp)                 :: tmp_sum(c % nnod, c % ng)

        tmp_sum = 0.0
        do h = 1, g-1
            do concurrent (n = 1:c % nnod)
                tmp_sum(n,h) = c % xs % sigs(n, h, g) * flux(n, h)
            end do
        end do

        do h = g+1, c % ng
            do concurrent (n = 1:c % nnod)
                tmp_sum(n,h) = c % xs % sigs(n, h, g) * flux(n, h)
            end do
        end do

        c % scat_src(:,g) = sum(tmp_sum, dim=2) 

        do concurrent (n = 1:c % nnod)
            total_src(n) = c % xs % chi(n, g) * c % fsrc(n) / c % Keff &
            + c % scat_src(n, g) + c % ext_src(n, g)
        end do
    
    end subroutine

    !===============================================================================================!
    ! calculate fission source
    !===============================================================================================!

    pure subroutine get_fission_src(c, flux)

        class(core_rect), intent(inout) :: c

        integer                 :: g, n
        real(dp), intent(in)    :: flux(:,:)
        real(dp)                :: tmp_sum(c % nnod, c % ng)

        tmp_sum = 0.0
        do g = 1, c % ng
            do concurrent (n = 1:c % nnod)
                tmp_sum(n,g) = flux(n, g) * c % xs % nuf(n, g)
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
    ! calculate FDM nodal coupling coefficients                                                     !
    !===============================================================================================!

    pure subroutine fdm_coupling_coef(c)

        class(core_rect), intent(inout) :: c

        integer  :: g
        real(dp) :: DN, D
        real(dp) :: xd, yd, zd
        real(dp) :: xdn, ydn, zdn
        integer  :: n, np

        do g = 1, c % ng
            do n = 1, c % nnod
    
                D  = c % xs % D(n, g)
                xd = c % xdel(n)
                yd = c % ydel(n)
                zd = c % zdel(n)
    
                ! Set FDM coupling coefficients in x direction
                if (.not. associated(c % node(n) % east)) then
                    if (xeast == zero_flux) then
                        c % node(n) % sf(1) % dt(g) = 2. * big * D / (2. * D + big * xd)
                    elseif (xeast == zero_incoming) then
                        c % node(n) % sf(1) % dt(g) =  D / (2. * D + 0.5 * xd)
                    else
                        c % node(n) % sf(1) % dt(g) =  0.
                    end if
                else
                    np  = c % node(n) % east % n
                    DN  = c % xs % D(np, g)
                    xdn = c % xdel(np)
                    c % node(n) % sf(1) % dt(g) =  2. * D * DN / (D * xdn +  DN * xd)
                endif
    
                if (.not. associated(c % node(n) % west)) then
                    if (xwest == zero_flux) then
                        c % node(n) % sf(2) % dt(g) = 2. * big * D / (2. * D + big * xd)
                    elseif (xwest == zero_incoming) then
                        c % node(n) % sf(2) % dt(g) =  D / (2. * D + 0.5 * xd)
                    else
                        c % node(n) % sf(2) % dt(g) =  0.
                    end if
                else
                    np  = c % node(n) % west % n
                    DN  = c % xs % D(np, g)
                    xdn = c % xdel(np)
                    c % node(n) % sf(2) % dt(g) =  2. * D * DN / (D * xdn +  DN * xd)
                endif
    
                ! Set FDM coupling coefficients in y direction
                if (.not. associated(c % node(n) % north)) then
                    if (ynorth == zero_flux) then
                        c % node(n) % sf(3) % dt(g) = 2. * big * D / (2. * D + big * yd)
                    elseif (ynorth == zero_incoming) then
                        c % node(n) % sf(3) % dt(g) =  D / (2. * D + 0.5 * yd)
                    else
                        c % node(n) % sf(3) % dt(g) =  0.
                    end if
                else
                    np  = c % node(n) % north % n
                    DN  = c % xs % D(np, g)
                    ydn = c % ydel(np)
                    c % node(n) % sf(3) % dt(g) =  2. * D * DN / (D * ydn +  DN * yd)
                endif
    
                if (.not. associated(c % node(n) % south)) then
                    if (ysouth == 0) then
                        c % node(n) % sf(4) % dt(g) = 2. * big * D / (2. * D + big * yd)
                    elseif (ysouth == 1) then
                        c % node(n) % sf(4) % dt(g) =  D / (2. * D + 0.5 * yd)
                    else
                        c % node(n) % sf(4) % dt(g) =  0.
                    end if
                else
                    np  = c % node(n) % south % n
                    DN  = c % xs % D(np, g)
                    ydn = c % ydel(np)
                    c % node(n) % sf(4) % dt(g) =  2. * D * DN / (D * ydn +  DN * yd)
                endif
    
                ! Set FDM coupling coefficients in z direction
                if (.not. associated(c % node(n) % top)) then
                    if (ztop == 0) then
                        c % node(n) % sf(5) % dt(g) = 2. * big * D / (2. * D + big * zd)
                    elseif (ztop == 1) then
                        c % node(n) % sf(5) % dt(g) =  D / (2. * D + 0.5 * zd)
                    else
                        c % node(n) % sf(5) % dt(g) =  0.
                    end if
                else
                    np  = c % node(n) % top % n
                    DN  = c % xs % D(np, g)
                    zdn = c % zdel(np)
                    c % node(n) % sf(5) % dt(g) =  2. * D * DN / (D * zdn +  DN * zd)
                endif
        
                if (.not. associated(c % node(n) % bottom)) then
                    if (zbott == 0) then
                        c % node(n) % sf(6) % dt(g) = 2. * big * D / (2. * D + big * zd)
                    elseif (zbott == 1) then
                        c % node(n) % sf(6) % dt(g) =  D / (2. * D + 0.5 * zd)
                    else
                        c % node(n) % sf(6) % dt(g) =  0.
                    end if
                else
                    np  = c % node(n) % bottom % n
                    DN  = c % xs % D(np, g)
                    zdn = c % zdel(np)
                    c % node(n) % sf(6) % dt(g) =  2. * D * DN / (D * zdn +  DN * zd)
                endif
            end do
        end do
    
    end subroutine

    !===============================================================================================!
    ! Set node data and setup the matrix for linear solver                                          !
    !===============================================================================================!

    subroutine set_nodes(c, kernel, mat_map, xdel, ydel, zdel, &
        sigtr, siga, nuf, sigf, sigs, chi)

        class(core_rect), target       :: c
        character(*), intent(in)       :: kernel
        integer, intent(in)            :: mat_map(:,:,:)
        real(dp), intent(in)           :: xdel(:), ydel(:), zdel(:)
        real(dp), intent(in)           :: sigtr(:,:)          ! Transport macroscopic XSEC
        real(dp), intent(in)           :: siga (:,:)          ! Absorption macroscopic XSEC
        real(dp), intent(in)           :: nuf  (:,:)          ! nu*fission macroscopic XSEC
        real(dp), intent(in)           :: sigf (:,:)          ! fission macroscopic XSEC
        real(dp), intent(in)           :: sigs (:,:,:)        ! Scattering macroscopic XSEC
        real(dp), intent(in)           :: chi (:,:)           ! fission spectrum

        integer :: i, j, k, n, ss
        integer :: ix(c % nnod), iy(c % nnod), iz(c % nnod)
        integer :: xyz(c % nxx, c % nyy, c % nzz)

        integer :: nnod
        integer :: nxx, nyy, nzz
        integer :: ng

        nnod = c % nnod
        nxx  = c % nxx
        nyy  = c % nyy
        nzz  = c % nzz
        ng   = c % ng

        c % kernel = kernel

        ! Set ix, iy, iz and xyz
        n = 0
        xyz = 0
        do k = 1, nzz
            do j = 1, nyy
                do i = c % ystag(j) % smin, c % ystag(j) % smax
                     n = n + 1
                     ix(n) = i
                     iy(n) = j
                     iz(n) = k
                     xyz(i,j,k) = n
                end do
            end do
        end do

        allocate(c % node(nnod))
        allocate(c % fsrc(c % nnod))

        ! set delta mesh
        allocate(c % vdel(nnod))
        allocate(c % xdel(nnod))
        allocate(c % ydel(nnod))
        allocate(c % zdel(nnod))
        n = 0
        do k = 1, nzz
            do j = 1, nyy
                do i = c % ystag(j) % smin, c % ystag(j) % smax
                    n = n + 1
                    c % xdel(n)     = xdel(i)
                    c % ydel(n)     = ydel(j)
                    c % zdel(n)     = zdel(k)
                    c % vdel(n)     = xdel(i) * ydel(j) * zdel(k)
                    c % node(n) % n = n
                end do
            end do
        end do

        ! allocate data
        allocate(c % scat_src(nnod, ng))
        allocate(c % ext_src(nnod, ng))
        c % ext_src = 0.0
        do ss = 1, n_rect_surf
            do n = 1, nnod
                allocate(c % node(n) % sf(ss) % dt(ng))
            end do
        end do
        
        ! to do: for FDM not necessary to allocate this
        do n = 1, nnod
            do ss = 1, n_rect_surf
                allocate(c % node(n) % sf(ss) % dh(ng))
                c % node(n) % sf(ss) % dt = 0.0
            end do
        end do

        ! set linked list
        do n = 1, nnod
            i = ix(n); j = iy(n); k = iz(n)

            if (i .ne. c % ystag(j) % smax) c % node(n) % east   => c % node(xyz(i+1, j, k))
            if (i .ne. c % ystag(j) % smin) c % node(n) % west   => c % node(xyz(i-1, j, k))
            if (j .ne. c % xstag(i) % smax) c % node(n) % north  => c % node(xyz(i, j+1, k))
            if (j .ne. c % xstag(i) % smin) c % node(n) % south  => c % node(xyz(i, j-1, k))
            if (k .ne. nzz)                 c % node(n) % top    => c % node(xyz(i, j, k+1))
            if (k .ne. 1)                   c % node(n) % bottom => c % node(xyz(i, j, k-1))
        end do

        ! Allocate xs
        call alloc_xsec(c % xs, nnod, ng, n_rect_surf)
        call set_xsec(c % xs, nnod, ix, iy, iz, mat_map, sigtr, siga, nuf, sigf, sigs, chi)

        ! Setup matrix
        call set_fdm_matrix(c, ix, iy, iz)

        if (c % kernel .ne. ' FDM') then
            allocate(c % d)
            call setup_nodal(c % d, c % ng, c % nnod, c % nxx, c % nyy, c % nzz, c % node, c % Keff, &
            c % xs, c % flux, c % ext_src, c % xdel, c % ydel, c % zdel, fid, xeast, xwest, &
            ynorth, ysouth, ztop, zbott, c % xstag, c % ystag, xyz, 'static', .false.)
            
            call alloc_data(c % d)
        end if

    end subroutine

    !===============================================================================================!
    ! Set FDM Matrix                                                                           !
    !===============================================================================================!

    subroutine set_fdm_matrix(c, ix, iy, iz)

        class(core_rect)          :: c
        integer, intent(in)  :: ix(:), iy(:), iz(:)

        integer  :: n
        integer  :: n_non_zero

        ! Count non zero elements
        n_non_zero = 0
        do n = 1, c % nnod
            if (associated(c % node(n) % bottom)) n_non_zero = n_non_zero + 1
            if (associated(c % node(n) % south)) n_non_zero = n_non_zero + 1
            if (associated(c % node(n) % west)) n_non_zero = n_non_zero + 1
            n_non_zero = n_non_zero + 1
            if (associated(c % node(n) % east)) n_non_zero = n_non_zero + 1
            if (associated(c % node(n) % north)) n_non_zero = n_non_zero + 1
            if (associated(c % node(n) % top)) n_non_zero = n_non_zero + 1
        end do

        call allocate_matrix(c % m, c % ng, c % nnod, n_non_zero)
        call set_fdm_index(c, ix, iy, iz)
        call setup_matrix(c, upd_fdm_coupling=.true.)

    end subroutine

    !===============================================================================================!
    ! Set FDM Matrix Index (indexed in CSR)                                                         !
    !===============================================================================================!

    subroutine set_fdm_index(c, ix, iy, iz)

        class(core_rect)         :: c
        integer, intent(in) :: ix(:), iy(:), iz(:)

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
            i = ix(n); j = iy(n); k = iz(n)
    
            c % m % ind % row(n) = idx
      
            ! Lower diagonal matrix element for z-direction
            if (associated(c % node(n) % bottom)) then
                c % m % ind % col(idx) = n - np
                idx                 = idx + 1
            end if
      
            ! Lower diagonal matrix element for y-direction
            if (associated(c % node(n) % south)) then
                c % m % ind % col(idx) = n - (nodp(i,j) - nodp(i,j-1))
                idx                 = idx + 1
            end if
      
            ! Lower diagonal matrix element for x-direction
            if (associated(c % node(n) % west)) then
                c % m % ind % col(idx) = n - 1
                idx                 = idx + 1
            end if
      
            ! Diagonal matrix elementss
            c % m % ind % col(idx) = n
            idx                 = idx + 1
      
             ! Upper diagonal matrix element for x-direction
            if (associated(c % node(n) % east)) then
                c % m % ind % col(idx) = n + 1
                idx                 = idx + 1
            end if
      
            ! Upper diagonal matrix element for y-direction
            if (associated(c % node(n) % north)) then
                c % m % ind % col(idx) = n + (nodp(i,j+1) - nodp(i,j))
                idx                    = idx + 1
            end if
      
            ! Upper diagonal matrix element for z-direction
            if (associated(c % node(n) % top)) then
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
        real(dp)                    :: xd, yd, zd
        type(surface_type), pointer :: su(:)

        integer :: nnod, ng

        nnod = c % nnod
        ng   = c % ng
    
        ! If need to calculate FDM coupling coefficients
        if (upd_fdm_coupling) call fdm_coupling_coef(c)
    
        ! Setup CMFD linear system
        do g = 1, ng
            idx = 0
            do n = 1, nnod

                xd = c % xdel(n)
                yd = c % ydel(n)
                zd = c % zdel(n)

                su => c % node(n) % sf
      
                ! Lower diagonal matrix element for z-direction
                if (associated(c % node(n) % bottom)) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(su(6) % dt(g) - su(6) % dh(g)) / zd
                end if
        
                ! Lower diagonal matrix element for y-direction
                if (associated(c % node(n) % south)) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(su(4) % dt(g) - su(4) % dh(g)) / yd
                end if
        
                ! Lower diagonal matrix element for x-direction
                if (associated(c % node(n) % west)) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(su(2) % dt(g) - su(2) % dh(g)) / xd
                end if
        
                ! Diagonal matrix elementss
                idx = idx + 1
                c % m % mtx(g) % elmn(idx) = &
                (su(1) % dt(g) + su(2) % dt(g) - su(1) % dh(g) + su(2) % dh(g)) / xd + &
                (su(3) % dt(g) + su(4) % dt(g) - su(3) % dh(g) + su(4) % dh(g)) / yd + &
                (su(5) % dt(g) + su(6) % dt(g) - su(5) % dh(g) + su(6) % dh(g)) / zd + &
                c % xs % sigr(n, g)
        
                 ! Upper diagonal matrix element for x-direction
                if (associated(c % node(n) % east)) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(su(1) % dt(g) + su(1) % dh(g)) / xd
                end if
        
                ! Upper diagonal matrix element for y-direction
                if (associated(c % node(n) % north)) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(su(3) % dt(g) + su(3) % dh(g)) / yd
                end if
        
                ! Upper diagonal matrix element for z-direction
                if (associated(c % node(n) % top)) then
                    idx = idx + 1
                    c % m % mtx(g) % elmn(idx) = -(su(5) % dt(g) + su(5) % dh(g)) / zd
                end if
      
            end do
        end do

    end subroutine

    !===============================================================================================!
    ! Set FDM parameters                                                                           !
    !===============================================================================================!

    subroutine set_fdm_parameters(c, inp_ng, inp_nxx, inp_nyy, inp_nzz, &
        nmat, mat_map, east, west, north, south, bottom, top, output_fid)

        class(core_rect)                   :: c
        integer, intent(in)           :: inp_ng                                   ! number of group
        integer, intent(in)           :: inp_nxx, inp_nyy, inp_nzz                ! number of nodes in x, y, and z directions
        integer, intent(in)           :: nmat                                     ! number of material
        integer, intent(in)           :: mat_map(:,:,:)                           ! 3D Material map
        integer, intent(in), optional :: east, west, north, south, bottom, top    ! BCs
        integer, intent(in), optional :: output_fid                               ! output file unit number

        integer :: i, j, k, n

        !set output file unit number
        if(present(output_fid)) fid = output_fid

        ! set number of energy group
        c % ng = inp_ng

        ! set number of nodes in x, y, and z directions
        c % nxx  = inp_nxx
        c % nyy  = inp_nyy
        c % nzz  = inp_nzz

        ! set BCs
        if (present(east))   xeast  = east
        if (present(west))   xwest  = west
        if (present(north))  ynorth = north
        if (present(south))  ysouth = south
        if (present(bottom)) zbott  = bottom
        if (present(top))    ztop   = top

        ! check number material
        if (maxval(mat_map) > nmat ) call fatal_error(fid, 'Material number ' &
        // n2c(maxval(mat_map)) // ' is greater than number of material')

        ! -Indexing non zero material for staggered mesh-
        allocate(c % ystag(c % nyy), c % xstag(c % nxx))

        !Indexing non zero material for staggered mesh along y direction
        do j= 1, c % nyy
            c % ystag(j) % smin = c % nxx
            do i = 1, c % nxx
                if (mat_map(i,j,1) /= 0) then
                    c % ystag(j) % smin = i
                    exit
                end if
            end do
        end do
        
        do j= 1, c % nyy
            c % ystag(j) % smax = 0
            do i = c % nxx, 1, -1
                if (mat_map(i,j,1) /= 0) then
                    c % ystag(j) % smax = i
                    exit
                end if
            end do
        end do
        
        !Indexing non zero material for staggered mesh along x direction
        do i= 1, c % nxx
            c % xstag(i) % smin = c % nyy
            do j = 1, c % nyy
                if (mat_map(i,j,1) /= 0) then
                    c % xstag(i) % smin = j
                    exit
                end if
            end do
        end do
        
        do i= 1, c % nxx
            c % xstag(i) % smax = 0
            do j = c % nyy, 1, -1
                if (mat_map(i,j,1) /= 0) then
                    c % xstag(i)%smax = j
                    exit
                end if
            end do
        end do
        
        ! Checking zero material between non-zero material
        do k = 1, c % nzz
            do j = 1, c % nyy
                do i = c % ystag(j) % smin, c % ystag(j) % smax
                    if (mat_map(i,j,k) == 0) call fatal_error(fid, &
                    'Zero material found inside core_rect. Check material assignment')
                end do
            end do
        end do
        do k = 1, c % nzz
            do i = 1, c % nxx
                do j = c % xstag(i) % smin, c % xstag(i) % smax
                    if (mat_map(i,j,k) == 0) call fatal_error(fid, &
                    'Zero material found inside core_rect. Check material assignment')
                end do
            end do
        end do

        ! set total number of nodes
        n = 0
        do k = 1, c % nzz
            do j = 1, c % nyy
                do i = c % ystag(j) % smin, c % ystag(j) % smax
                     n = n + 1
                end do
            end do
        end do

        c % nnod = n
        
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
        ! real(dp), allocatable :: ext_src(:)        ! external source
        ! real(dp), allocatable :: omega(:)          ! Exponential transformation constant
        ! real(dp), allocatable :: sigrp(:)          ! Initial removal cross sections before added by parameters required for transient
        ! real(dp), allocatable :: L    (:)          ! Total leakages for node n and group g
        ! real(dp), allocatable :: dfis               

        ! integer               :: i_mat                ! material index

    
end module
    