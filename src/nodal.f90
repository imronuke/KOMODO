module nodal

    !-------- References -----------
    ! 1. Zimin V.G. and Ninokata, H., (1997) 
    ! Nonlinear Iteration Procedure Based on Legendre Polynomials, 
    ! Trans. Am. Nucl. Soc., 76, 162.
    
    ! 2. Zimin, V. G., & Ninokata, H. (1998). 
    ! Nodal neutron kinetics model based on nonlinear iteration procedure for LWR analysis. 
    ! Annals of Nuclear Energy, 25(8), 507–528. https://doi.org/10.1016/S0306-4549(97)00078-9
    ! ------------------------------

    use iso_fortran_env, only: real64
    use node
    use xsec
    use util

    implicit none

    private

    save

    integer, parameter :: dp = real64

    integer, parameter  :: zero_flux = 0
    integer, parameter  :: zero_incoming = 1
    integer, parameter  :: reflective = 2

    integer :: fid                          ! file unit number

    integer :: xwest, xeast                 ! West and East Boundary conditions
    integer :: ysouth, ynorth               ! South and North Boundary conditions
    integer :: zbott, ztop                  ! Bottom and Top Boundary conditions
  
  ! Nodal data for D hat update
    type, public :: nodal_type
        character(9)         :: calc_mode
        logical              :: only_axial
        integer              :: nnod, ng
        integer              :: nxx, nyy, nzz
        integer, allocatable :: xyz(:,:,:)

        real(dp), allocatable   :: a1n(:), a2n(:), a3n(:), a4n(:)  ! Flux expansion coefficients for current node
        real(dp), allocatable   :: a1p(:), a2p(:), a3p(:), a4p(:)  ! Flux expansion coefficients for following node

        real(dp), allocatable :: Bcn(:,:), Bcp(:,:)                       ! Buckling for left and right nodes
        real(dp), allocatable :: An(:), Bn(:), En(:), Fn(:), Gn(:), Hn(:)
        real(dp), allocatable :: Ap(:), Bp(:), Ep(:), Fp(:), Gp(:), Hp(:)
        real(dp), allocatable :: Lm2(:)                                   ! Second order source
        real(dp), allocatable :: Sx(:,:), Sy(:,:), Sz(:,:)                ! zeroth order source in x, y, and z directions
        real(dp), allocatable :: Ln1(:), Lp1(:)                           ! First Transverse leakages moments
        
        type(node_rect), pointer :: node(:)
        real(dp), pointer        :: exsrc(:,:)
        real(dp), pointer        :: flux(:,:)
        real(dp), pointer        :: Ke
        type(xs_rect), pointer   :: xs            ! node-wise xs
        real(dp), pointer        :: xdel(:)       ! Delta x in cm
        real(dp), pointer        :: ydel(:)       ! Delta y in cm
        real(dp), pointer        :: zdel(:)       ! Delta z in cm
        type(staggered), pointer :: xstag(:)    ! Staggered mesh data
        type(staggered), pointer :: ystag(:)    ! Staggered mesh data
    end type

    public :: setup_nodal, alloc_data, dealloc_data, nodal_update_sanm, nodal_update_pnm

    contains

    !===============================================================================================!
    ! Allocate expansion coefficients and zeroth-moment source
    ! Aliasing external data
    !===============================================================================================!

    subroutine setup_nodal(d, ng, nnod, nxx, nyy, nzz, nod, Keff, xs, flux, exsrc, &
        xdel, ydel, zdel, fid_inp, east, west, north, south, top, bottom, xstag, ystag, xyz, &
        calc_mode, only_axial)
        
        class(nodal_type), intent(inout)       :: d
        integer, intent(in)                    :: ng, nnod
        integer, intent(in)                    :: nxx, nyy, nzz
        type(node_rect), intent(inout), target :: nod(:)
        real(dp), intent(in), target           :: Keff
        type(xs_rect), intent(inout), target   :: xs
        real(dp), intent(in), target           :: flux(:,:)
        real(dp), intent(in), target           :: exsrc(:,:)
        real(dp), intent(in), target           :: xdel(:), ydel(:), zdel(:)
        integer, intent(in)                    :: fid_inp
        integer, intent(in)                    :: east, west
        integer, intent(in)                    :: north, south
        integer, intent(in)                    :: top, bottom
        type(staggered), intent(in), target    :: xstag(:), ystag(:)
        integer, intent(in)                    :: xyz(:,:,:)
        character(*), intent(in)               :: calc_mode    
        logical, intent(in)                    :: only_axial

        d % nnod       = nnod
        d % ng         = ng
        d % nxx        = nxx
        d % nyy        = nyy
        d % nzz        = nzz
        d % calc_mode  = trim(adjustl(calc_mode))
        d % only_axial = only_axial

        allocate(d % xyz(nxx, nyy, nzz))
        d % xyz = xyz

        d % node  => nod
        d % exsrc => exsrc
        d % flux  => flux
        d % Ke    => Keff
        d % xs    => xs
        d % xdel  => xdel
        d % ydel  => ydel
        d % zdel  => zdel
        d % xstag => xstag
        d % ystag => ystag

        allocate(d % a1n(ng), d % a2n(ng), d % a3n(ng), d % a4n(ng))
        allocate(d % a1p(ng), d % a2p(ng), d % a3p(ng), d % a4p(ng))
        allocate(d % Sx(nnod,ng), d % Sy(nnod,ng), d % Sz(nnod,ng))

       ! set BCs
        xeast  = east
        xwest  = west
        ynorth = north
        ysouth = south
        zbott  = bottom
        ztop   = top

        fid = fid_inp

    end subroutine

    !===============================================================================================!
    ! Allocate necessary data for nodal coef update
    !===============================================================================================!

    subroutine alloc_data(d)
        
        class(nodal_type), intent(inout) :: d

        integer :: ng

        ng = d % ng

        allocate(d % Ln1(ng), d % Lp1(ng))
        allocate(d % Bcn(ng, ng), d % Bcp(ng, ng))
        allocate(d % An(ng), d % Bn(ng), d % En(ng), d % Fn(ng), d % Gn(ng), d % Hn(ng))
        allocate(d % Ap(ng), d % Bp(ng), d % Ep(ng), d % Fp(ng), d % Gp(ng), d % Hp(ng))
        allocate(d % Lm2(ng))

    end subroutine

    !===============================================================================================!
    ! Dellocate necessary data for nodal coef update
    !===============================================================================================!

    subroutine dealloc_data(d)
        
        class(nodal_type), intent(inout) :: d
        
        deallocate(d % Ln1, d % Lp1)
        deallocate(d % Bcn, d % Bcp)
        deallocate(d % An, d % Bn, d % En, d % Fn, d % Gn, d % Hn)
        deallocate(d % Ap, d % Bp, d % Ep, d % Fp, d % Gp, d % Hp)
        deallocate(d % Lm2)

    end subroutine

    !===============================================================================================!
    ! to update nodal coupling coefficients (D hat) for entire core using
    ! Polynomial Nodal Method (PNM)
    !===============================================================================================!

    subroutine nodal_update_pnm(d)

        class(nodal_type) :: d
    
        integer  :: n, p
        integer  :: i, j, k

        d % An = 1. / 15.; d % Bn = 1. / 35.; d % En = 2. / 7.
        d % Fn = 2. / 5.; d % Gn = 10.; d % Hn = 6.
        d % An = d % Ap; d % Bn = d % Bp; d % En = d % Ep
        d % Fn = d % Fp; d % Gn = d % Gp; d % Hn = d % Hp
    
        call get_source(d)
    
        !Node sweeps in x-direction
        do k = 1, d % nzz
            do j = 1, d % nyy
                p = d % xyz(d % ystag(j) % smin, j, k)
                d % Bcp = get_buckling(d, 1, p)
                call get_coefs_first(d, xwest, 1, p)
                call nodal_coup_upd(d, 1, d % a1p, d % a2p, d % a3p, d % a4p, p = p)
                do i = d % ystag(j) % smin, d % ystag(j) % smax-1
                    n = d % xyz(i,j,k); p = d % xyz(i+1,j,k)
                    d % Bcn = d % Bcp
                    d % An = d % Ap; d % Bn = d % Bp; d % En = d % Ep
                    d % Fn = d % Fp; d % Gn = d % Gp; d % Hn = d % Hp
                    d % Bcp = get_buckling(d, 1, p)
                    call get_coefs(d, 1, n, p)
                    call nodal_coup_upd(d, 1, d % a1n, d % a2n, d % a3n, d % a4n, n, p)
                end do
                n = d % xyz(d % ystag(j) % smax, j, k)
                d % Bcn = d % Bcp
                d % An = d % Ap; d % Bn = d % Bp; d % En = d % Ep
                d % Fn = d % Fp; d % Gn = d % Gp; d % Hn = d % Hp
                call get_coefs_last(d, xeast, 1, n)
                call nodal_coup_upd(d, 1, d % a1n, d % a2n, d % a3n, d % a4n, n = n)
            end do
        end do
    
        !Node sweeps in y-direction
        do k = 1, d % nzz
            do i = 1, d % nxx
                p = d % xyz(i, d % xstag(i) % smin, k)
                d % Bcp = get_buckling(d, 2, p)
                call get_coefs_first(d, ysouth, 2, p)
                call nodal_coup_upd(d, 2, d % a1p, d % a2p, d % a3p, d % a4p, p = p)
                do j = d % xstag(i) % smin, d % xstag(i) % smax-1
                    n = d % xyz(i,j,k); p = d % xyz(i,j+1,k)
                    d % Bcn = d % Bcp
                    d % An = d % Ap; d % Bn = d % Bp; d % En = d % Ep
                    d % Fn = d % Fp; d % Gn = d % Gp; d % Hn = d % Hp
                    d % Bcp = get_buckling(d, 2, p)
                    call get_coefs(d, 2, n, p)
                    call nodal_coup_upd(d, 2, d % a1n, d % a2n, d % a3n, d % a4n, n, p)
                end do
                n = d % xyz(i, d % xstag(i) % smax, k)
                d % Bcn = d % Bcp
                d % An = d % Ap; d % Bn = d % Bp; d % En = d % Ep
                d % Fn = d % Fp; d % Gn = d % Gp; d % Hn = d % Hp
                call get_coefs_last(d, ynorth, 2, n)
                call nodal_coup_upd(d, 2, d % a1n, d % a2n, d % a3n, d % a4n, n = n)
            end do
        end do
    
        !Node sweeps in z-direction
        do j = 1, d % nyy
            do i = d % ystag(j) % smin, d % ystag(j) % smax
                p = d % xyz(i, j, 1)
                d % Bcp = get_buckling(d, 3, p)
                call get_coefs_first(d, zbott, 3, p)
                call nodal_coup_upd(d, 3, d % a1p, d % a2p, d % a3p, d % a4p, p = p)
                do k = 1, d % nzz-1
                    n = d % xyz(i, j, k); p = d % xyz(i, j, k+1)
                    d % Bcn = d % Bcp
                    d % An = d % Ap; d % Bn = d % Bp; d % En = d % Ep
                    d % Fn = d % Fp; d % Gn = d % Gp; d % Hn = d % Hp
                    d % Bcp = get_buckling(d, 3, p)
                    call get_coefs(d, 3, n, p)
                    call nodal_coup_upd(d, 3, d % a1n, d % a2n, d % a3n, d % a4n, n, p)
                end do
                n = d % xyz(i, j, d % nzz)
                d % Bcn = d % Bcp
                d % An = d % Ap; d % Bn = d % Bp; d % En = d % Ep
                d % Fn = d % Fp; d % Gn = d % Gp; d % Hn = d % Hp
                call get_coefs_last(d, ztop, 3, n)
                call nodal_coup_upd(d, 3, d % a1n, d % a2n, d % a3n, d % a4n, n = n)
            end do
        end do
     
    end subroutine

    !===============================================================================================!
    ! to update nodal coupling coefficients (D hat) for entire core using
    ! Semi-Analytic Nodal Method (SANM)
    !===============================================================================================!

    subroutine nodal_update_sanm(d)

        class(nodal_type) :: d
    
        integer  :: n, p
        integer  :: i, j, k
    
        call get_source(d)
    
        !Node sweeps in x-direction
        do k = 1, d % nzz
            do j = 1, d % nyy
                p = d % xyz(d % ystag(j) % smin, j, k)
                d % Bcp = get_buckling(d, 1, p)
                call get_ABEFGH (d, p, 1, d % Ap, d % Bp, d % Ep, d % Fp, d % Gp, d % Hp)
                call get_coefs_first(d, xwest, 1, p)
                call nodal_coup_upd(d, 1, d % a1p, d % a2p, d % a3p, d % a4p, p = p)
                do i = d % ystag(j) % smin, d % ystag(j) % smax-1
                    n = d % xyz(i,j,k); p = d % xyz(i+1,j,k)
                    d % Bcn = d % Bcp
                    d % An = d % Ap; d % Bn = d % Bp; d % En = d % Ep
                    d % Fn = d % Fp; d % Gn = d % Gp; d % Hn = d % Hp
                    d % Bcp = get_buckling(d, 1, p)
                    call get_ABEFGH (d, p, 1, d % Ap, d % Bp, d % Ep, d % Fp, d % Gp, d % Hp)
                    call get_coefs(d, 1, n, p)
                    call nodal_coup_upd(d, 1, d % a1n, d % a2n, d % a3n, d % a4n, n, p)
                end do
                n = d % xyz(d % ystag(j) % smax, j, k)
                d % Bcn = d % Bcp
                d % An = d % Ap; d % Bn = d % Bp; d % En = d % Ep
                d % Fn = d % Fp; d % Gn = d % Gp; d % Hn = d % Hp
                call get_coefs_last(d, xeast, 1, n)
                call nodal_coup_upd(d, 1, d % a1n, d % a2n, d % a3n, d % a4n, n = n)
            end do
        end do
    
        !Node sweeps in y-direction
        do k = 1, d % nzz
            do i = 1, d % nxx
                p = d % xyz(i, d % xstag(i) % smin, k)
                d % Bcp = get_buckling(d, 2, p)
                call get_ABEFGH (d, p, 2, d % Ap, d % Bp, d % Ep, d % Fp, d % Gp, d % Hp)
                call get_coefs_first(d, ysouth, 2, p)
                call nodal_coup_upd(d, 2, d % a1p, d % a2p, d % a3p, d % a4p, p = p)
                do j = d % xstag(i) % smin, d % xstag(i) % smax-1
                    n = d % xyz(i,j,k); p = d % xyz(i,j+1,k)
                    d % Bcn = d % Bcp
                    d % An = d % Ap; d % Bn = d % Bp; d % En = d % Ep
                    d % Fn = d % Fp; d % Gn = d % Gp; d % Hn = d % Hp
                    d % Bcp = get_buckling(d, 2, p)
                    call get_ABEFGH (d, p, 2, d % Ap, d % Bp, d % Ep, d % Fp, d % Gp, d % Hp)
                    call get_coefs(d, 2, n, p)
                    call nodal_coup_upd(d, 2, d % a1n, d % a2n, d % a3n, d % a4n, n, p)
                end do
                n = d % xyz(i, d % xstag(i) % smax, k)
                d % Bcn = d % Bcp
                d % An = d % Ap; d % Bn = d % Bp; d % En = d % Ep
                d % Fn = d % Fp; d % Gn = d % Gp; d % Hn = d % Hp
                call get_coefs_last(d, ynorth, 2, n)
                call nodal_coup_upd(d, 2, d % a1n, d % a2n, d % a3n, d % a4n, n = n)
            end do
        end do
    
        !Node sweeps in z-direction
        do j = 1, d % nyy
            do i = d % ystag(j) % smin, d % ystag(j) % smax
                p = d % xyz(i, j, 1)
                d % Bcp = get_buckling(d, 3, p)
                call get_ABEFGH (d, p, 3, d % Ap, d % Bp, d % Ep, d % Fp, d % Gp, d % Hp)
                call get_coefs_first(d, zbott, 3, p)
                call nodal_coup_upd(d, 3, d % a1p, d % a2p, d % a3p, d % a4p, p = p)
                do k = 1, d % nzz-1
                    n = d % xyz(i, j, k); p = d % xyz(i, j, k+1)
                    d % Bcn = d % Bcp
                    d % An = d % Ap; d % Bn = d % Bp; d % En = d % Ep
                    d % Fn = d % Fp; d % Gn = d % Gp; d % Hn = d % Hp
                    d % Bcp = get_buckling(d, 3, p)
                    call get_ABEFGH (d, p, 3, d % Ap, d % Bp, d % Ep, d % Fp, d % Gp, d % Hp)
                    call get_coefs(d, 3, n, p)
                    call nodal_coup_upd(d, 3, d % a1n, d % a2n, d % a3n, d % a4n, n, p)
                end do
                n = d % xyz(i, j, d % nzz)
                d % Bcn = d % Bcp
                d % An = d % Ap; d % Bn = d % Bp; d % En = d % Ep
                d % Fn = d % Fp; d % Gn = d % Gp; d % Hn = d % Hp
                call get_coefs_last(d, ztop, 3, n)
                call nodal_coup_upd(d, 3, d % a1n, d % a2n, d % a3n, d % a4n, n = n)
            end do
        end do
     
    end subroutine

    !===============================================================================================!
    ! to update nodal coupling coefficients (D hat)
    !===============================================================================================!

    subroutine nodal_coup_upd(d, u, a1, a2, a3, a4, n, p)
    
        class(nodal_type), intent(inout), target :: d
        integer, intent(in)                      :: u
        integer, intent(in), optional            :: n, p
        real(dp), dimension(:), intent(in)       :: a1, a2, a3, a4
    
        integer    :: g, s
        real(dp)   :: dh, jp
        type(surface_type), pointer :: sf, sfp
    
        if (u == 1) then
            if (present(n)) then
                dh = d % xdel(n)
            else
                dh = d % xdel(p)
            end if
            s = 1
        else if (u == 2) then
            if (present(n)) then
                dh = d % ydel(n)
            else
                dh = d % ydel(p)
            end if
            s = 3
        else
            if (present(n)) then
              dh = d % zdel(n)
            else
              dh = d % zdel(p)
            end if
            s = 5
        end if
    
        if (present(n) .and. present(p)) then

            sf => d % node(n) % sf(s)     ! n means current node
            sfp => d % node(p) % sf(s+1)  ! p means next node

            do g = 1, d % ng

                jp = -2. * d % xs % D(n,g) / dh &
                * (a1(g) + 3. * a2(g) + d % Hn(g) * a3(g) + d % Gn(g) * a4(g))
        
                ! Update nodal coupling
                sf % dh(g) = (sf % dt(g) * (d % flux(n,g) - d % flux(p,g)) - jp) &
                / (d % flux(n,g) + d % flux(p,g))
                
                sfp % dh(g) = sf % dh(g)

            end do

        else if (present(p)) then

            sfp => d % node(p) % sf(s+1)  ! p means next node

            do g = 1, d % ng
                jp = -2. * d % xs % D(p,g) / dh &
                * (a1(g) - 3. * a2(g) + d % Hp(g) * a3(g) - d % Gp(g) * a4(g))
        
                ! Update nodal coupling
                sfp % dh(g) = -(jp / d % flux(p,g) + sfp % dt(g))

            end do

        else

            sf => d % node(n) % sf(s)

            do g = 1, d % ng
                jp = -2. * d % xs % D(n,g) / dh &
                * (a1(g) + 3. * a2(g) + d % Hn(g) * a3(g) + d % Gn(g) * a4(g))
        
                ! Update nodal coupling
                sf % dh(g) = -(jp / d % flux(n,g) - sf % dt(g))
        
            end do

        end if
    
      end subroutine nodal_coup_upd

    !===============================================================================================!
    ! calculate flux expansion coefficients for the most left node (first one) (see ref 1)
    !===============================================================================================!

    subroutine get_coefs_first(d, bc, u, p)

        class(nodal_type), intent(inout) :: d
        integer, intent(in)              :: bc
        integer, intent(in)              :: u, p
    
        real(dp) :: A(d % ng, d % ng)   ! GxG Matrix
        real(dp) :: b(d % ng)           ! G vector
    
        !Setup GxG matrix and G vector to obtain a2(g) for node p
        call get_a2matvec(d, u, p, A, b)

        !calculate a2 expansion coefficients4
        d % a2p = LU_solve(p, d % ng, A, b)

        !calculate a4 expansion coefficients
        d % a4p = get_a4(d, d % a2p)
    
        !Setup GxG matrix and G vector to obtain a1(g) for left node
        call get_a1matvec_first(d, bc, u, p, d % a2p, d % a4p, A, b)
    
        !calculate a1 expansion coefficients
        d % a1p = LU_solve(p, d % ng, A, b)
    
        !calculate a3 expansion coefficients
        d % a3p = get_a3(d, 2, d % a1p, d % Lp1)

    end subroutine

    !===============================================================================================!
    ! calculate flux expansion coefficients for the most right node (last one) (see ref 1)
    !===============================================================================================!

    subroutine get_coefs_last(d, bc, u, n)

        class(nodal_type), intent(inout) :: d
        integer, intent(in)              :: bc
        integer, intent(in)              :: u, n
    
        real(dp) :: A(d % ng, d % ng)   ! GxG Matrix
        real(dp) :: b(d % ng)           ! G vector
    
        !get a2n and a4n expansion coefficients
        d % a2n = d % a2p
        d % a4n = d % a4p
    
        !Setup GxG matrix and G vector to obtain a1(g) for most right node
        call get_a1matvec_last(d, bc, u, n, d % a2n, d % a4n, A, b)
    
        !calculate a1 expansion coefficients
        d % a1n = LU_solve(n, d % ng, A, b)
    
        !calculate a3 expansion coefficients
        d % a3n = get_a3(d, 1, d % a1n, d % Ln1)

    end subroutine

    !===============================================================================================!
    ! calculate flux expansion coefficients (see ref 1)
    !===============================================================================================!

    subroutine get_coefs(d, u, n, p)

        class(nodal_type), intent(inout) :: d
        integer, intent(in)              :: u, n, p
    
        real(dp)  :: A(d % ng, d % ng)          ! GxG Matrix
        real(dp)  :: b(d % ng)                  ! G vector
        real(dp)  :: R(2 * d % ng, 2 * d % ng)  ! 2Gx2G Matrix
        real(dp)  :: s(2 * d % ng)              ! 2G vector
    
        !get a2n and a4n expansion coefficients
        d % a2n = d % a2p
        d % a4n = d % a4p
    
        !Setup GxG matrix and G vector to obtain a2(g) for node p
        call get_a2matvec(d, u, p, A, b)

        !calculate a2 expansion coefficients
        d % a2p = LU_solve(p, d % ng, A, b)

        !calculate a4 expansion coefficients
        d % a4p = get_a4(d, d % a2p)
    
        !Setup 2Gx2G matrix and 2G vector to obtain a1(g) for node n
        call get_a1matvec(d, u, n, p, d % a2n, d % a4n, d % a2p, d % a4p, R, s)
    
        !calculate a1 expansion coefficients
        s = LU_solve(n, 2 * d % ng, R, s)
        d % a1n = s(1 : d % ng)
    
        !calculate a3 expansion coefficients
        d % a3n = get_a3(d, 1, d % a1n, d % Ln1)

    end subroutine

    !===============================================================================================!
    ! To get a3 expansion coefficients
    !===============================================================================================!

    function get_a3(d, cp, a1, Lmn1) result (a3)

        class(nodal_type), intent(in) :: d
        integer, intent(in)           :: cp       ! to indicate whether it is first node (=2 means first node)
        real(dp), intent(in)          :: a1(:)       ! a1 expansion coefficients
        real(dp), intent(in)          :: Lmn1(:)     ! First transverse leakage moments
        real(dp)                      :: a3(d % ng)       ! a3 expansion coefficients
    
        real(dp) :: Bf
        real(dp), dimension(d % ng, d % ng)  :: Bc
        real(dp), dimension(d % ng)          :: Ac
        integer :: g, h
    
        if (cp == 1) then
            Bc = d % Bcn
            Ac = d % An
        else
            Bc = d % Bcp    ! Take buckling values from the left node
            Ac = d % Ap
        end if
    
        do g = 1, d % ng
            Bf = 0.
            do h = 1, d % ng
                Bf = Bf + Bc(g,h)*a1(h)
            end do
            a3(g) = Ac(g) * (Bf + Lmn1(g))
        end do
    
      end function

    !===============================================================================================!
    ! To get a4 expansion coefficients
    !===============================================================================================!

    function get_a4(d, a2) result (a4)

        class(nodal_type), intent(in)    :: d
        real(dp), intent(in)             :: a2(:)       ! a2 expansion coefficients
        real(dp)                         :: a4(d % ng)  ! a4 expansion coefficients
    
        real(dp) :: Bf
        integer :: g, h
    
        do g = 1, d % ng
            Bf = 0.
            do h = 1, d % ng
                Bf = Bf + d % Bcp(g,h) * a2(h)
            end do
      
            a4(g) = d % Bp(g) * (Bf + d % Lm2(g))
        end do
    
    end function

    !===============================================================================================!
    ! construct LHS and RHS of tp obtain a1 expansion coefficients for most left node (see ref 1)
    !===============================================================================================!

    subroutine get_a1matvec_first(d, bc, u, p, a2p, a4p, M, rhs)
    
        class(nodal_type), intent(inout) :: d
        
        integer, intent(in)              :: bc
        integer, intent(in)              :: u, p              ! Direction and node n umber
        real(dp), intent(in)             :: a2p(:), a4p(:)    ! a2 and a4 expansion coefficients
        real(dp), intent(out)            :: M(:,:)            ! 2Gx2G Matrix (LHS)
        real(dp), intent(out)            :: rhs(:)              ! 2G vector (RHS)
       
        real(dp)  :: dn
        real(dp)  :: dc
        real(dp)  :: Ap, Fp, Gp, Hp, a2, a4, L1
        real(dp)  :: Pp                       !dif. coef/dx
        integer   :: g, h, sf
    
        ! define node size in direction u
        if (u == 1) then
            dn = d % xdel(p)
            sf = 2
        else if (u == 2) then
            dn = d % ydel(p)
            sf = 4
        else
            dn = d % zdel(p)
            sf = 6
        end if
    
        do g = 1, d % ng
            !Calculate transverse leakage first moment
            call TLUpd1(d, u, p, g, d % Lp1(g))
            Pp = 2. * d % xs % D(p,g) / dn

            dc = d % xs % dc(p,g,sf)
            Ap = d % Ap(g)
            Fp = d % Fp(g)
            Gp = d % Gp(g)
            Hp = d % Hp(g)
            a2 = a2p(g)
            a4 = a4p(g)
            L1 = d % Lp1(g)
      
            !Setup GxG matrix A and G vector B to obtain a2(g)
            if (bc == reflective) then
                do h = 1, d % ng
                    if (h == g) then
                        M(g,g) = Pp * (d % Bcp(g,h) * Fp + 1.)
                    else
                        M(g,h) = Pp * d % Bcp(g,h) * Fp
                    end if
                end do
                rhs(g) = Pp * (3. * a2 + Gp * a4 - Fp * L1)
            else if (bc == zero_incoming) then
                do h = 1, d % ng
                    if (h == g) then
                        M(g,g) = -dc * (1. + Ap * d % Bcp(g,h)) - 2. * Pp * (Ap * d % Bcp(g,h) * Hp + 1.)
                    else
                        M(g,h) = -dc * Ap * d % Bcp(g,h) - 2. * Pp * Ap * d % Bcp(g,h) * Hp
                    end if
                end do
                rhs(g) = 2. * Pp * (Ap * Hp * L1 - 3. * a2 - Gp * a4) - dc * (a2 + a4 + d % flux(p,g) - Ap * L1)
            else
                do h = 1, d % ng
                    if (h == g) then
                        M(g,g) = dc * (1. + Ap * d % Bcp(g,h))
                    else
                        M(g,h) = dc * Ap * d % Bcp(g,h)
                    end if
                end do
                rhs(g) =  dc * (a2 + a4 + d % flux(p,g) - Ap * L1)
            end if
        end do
    
    end subroutine

    !===============================================================================================!
    ! construct LHS and RHS of tp obtain a1 expansion coefficients for most right node (see ref 1)
    !===============================================================================================!

    subroutine get_a1matvec_last(d, bc, u, n, a2n, a4n, M, rhs)
    
        class(nodal_type), intent(inout) :: d
        
        integer, intent(in)              :: bc
        integer, intent(in)              :: u, n              ! Direction and node n umber
        real(dp), intent(in)             :: a2n(:), a4n(:)    ! a2 and a4 expansion coefficients
        real(dp), intent(out)            :: M(:,:)            ! 2Gx2G Materix (LHS)
        real(dp), intent(out)            :: rhs(:)              ! 2G vector (RHS)
       
        real(dp)   :: dn
        real(dp)   :: dc
        real(dp)   :: An, Fn, Hn, Gn, L1, a2, a4
        real(dp)   :: Pn                       !dif. coef/dx
        integer    :: g, h, sf
    
        ! define node size in direction u
        if (u == 1) then
            dn = d % xdel(n)
            sf = 1
        else if (u == 2) then
            dn = d % ydel(n)
            sf = 3
        else
            dn = d % zdel(n)
            sf = 5
        end if
  
        do g = 1, d % ng
            d % Ln1(g) = d % Lp1(g)
            Pn  = 2. * d % xs % D(n,g) / dn

            dc = d % xs % dc(n,g,sf)
            An = d % An(g)
            Fn = d % Fn(g)
            Hn = d % Hn(g)
            Gn = d % Gn(g)
            a2 = a2n(g)
            a4 = a4n(g)
            L1 = d % Ln1(g)
      
            !Setup 2Gx2G matrix A and 2G vector B to obtain a2(g)
            if (bc == reflective) then
                do h = 1, d % ng
                    if (h == g) then
                        M(g,g) = -Pn * (d % Bcn(g,h) * Fn + 1.)
                    else
                        M(g,h) = -Pn * d % Bcn(g,h) * Fn
                    end if
                end do
                rhs(g) = Pn * (3. * a2 + Gn * a4 + Fn * L1)
            else if (bc == zero_incoming) then
                do h = 1, d % ng
                    if (h == g) then
                        M(g,g) = dc * (1. + An * d % Bcn(g,h)) + 2. * Pn * (An * d % Bcn(g,h) * Hn + 1.)
                    else
                        M(g,h) = dc * An * d % Bcn(g,h) + 2. * Pn * An * d % Bcn(g,h) * Hn
                    end if
                end do
                rhs(g) = -2. * Pn * (An * Hn * L1 + 3. * a2 + Gn * a4) - dc * (a2 + a4 + d % flux(n,g) + An * L1)
            else
                do h = 1, d % ng
                    if (h == g) then
                        M(g,g) = dc * (1. + An * d % Bcn(g,h))
                    else
                        M(g,h) = dc * An * d % Bcn(g,h)
                    end if
                end do
                rhs(g) = -dc * (a2 + a4 + d % flux(n,g) + An * L1)
            end if
        end do
    
    end subroutine

    !===============================================================================================!
    ! construct LHS and RHS of tp obtain a1 expansion coefficients (see ref 1)
    !===============================================================================================!

    subroutine get_a1matvec(d, u, n, p, a2n, a4n, a2p, a4p, M, rhs)
    
        class(nodal_type), intent(inout) :: d
        integer, intent(in)               :: u, n, p                  ! Direction and node number
        real(dp), intent(in)              :: a2n(:), a4n(:)           ! a2 and a4 expansion coefficients
        real(dp), intent(in)              :: a2p(:), a4p(:)           ! a2 and a4 expansion coefficients
        real(dp), intent(out)             :: M(:,:)           ! 2Gx2G Materix
        real(dp), intent(out)             :: rhs(:)             ! 2G vector
    
        real(dp)    :: hn, hp
        real(dp)    :: Pn(d % ng), Pp(d % ng)              !dif. coef/dx
        integer     :: ng, g, h, sf

        ng = d % ng
    
        ! define node size in direction u
        if (u == 1) then
            hn = d % xdel(n); hp = d % xdel(p)
            sf = 1
        else if (u == 2) then
            hn = d % ydel(n); hp = d % ydel(p)
            sf = 3
        else
            hn = d % zdel(n); hp = d % zdel(p)
            sf = 5
        end if
  
        do g = 1, ng
            ! Calculate transverse leakage moments
            d % Ln1(g) = d % Lp1(g)
            call TLUpd1 (d, u, p, g, d % Lp1(g))
      
            !Setup 2Gx2G matrix A and 2G vector B to obtain a2(g)
            Pn(g)  = 2. * d % xs % D(n,g) / hn; Pp(g) = 2. * d % xs % D(p,g) / hp

            do h = 1, ng
                if (h == g) then
                    M(g,g)    = -Pn(g) * (d % Bcn(g,h) * d % Fn(g) + 1.)
                    M(g,g+ng) =  Pp(g) * (d % Bcp(g,h) * d % Fp(g) + 1.)
                else
                    M(g,h)    = -Pn(g) * d % Bcn(g,h) * d % Fn(g)
                    M(g,h+ng) =  Pp(g) * d % Bcp(g,h) * d % Fp(g)
                end if
            end do
            rhs(g) = Pn(g) * (3. * a2n(g) + d % Gn(g) * a4n(g) + d % Fn(g) * d % Ln1(g)) &
            + Pp(g) * (3. * a2p(g) + d % Gp(g) * a4p(g) - d % Fp(g) * d % Lp1(g))
        end do

        do g = 1, ng
            do h = 1, ng
                if (h == g) then
                    M(g+ng,g)    = d % xs % dc(n,g,sf)   * (d % Bcn(g,h) * d % An(g) + 1.)
                    M(g+ng,g+ng) = d % xs % dc(p,g,sf+1) * (d % Bcp(g,h) * d % Ap(g) + 1.)
                else
                    M(g+ng,h)    = d % xs % dc(n,g,sf)   * d % Bcn(g,h) * d % An(g)
                    M(g+ng,h+ng) = d % xs % dc(p,g,sf+1) * d % Bcp(g,h) * d % Ap(g)
                end if
            end do
            ! Create vector b
            rhs(g+ng) = d % xs % dc(p,g,sf+1) * (a2p(g) + a4p(g) + d % flux(p,g) - d % An(g) * d % Ln1(g)) &
            - d % xs % dc(n,g,sf) * (a2n(g) + a4n(g) + d % flux(n,g) + d % Ap(g) * d % Lp1(g))
        end do
    
    end subroutine

    !===============================================================================================!
    ! construct LHS and RHS of tp obtain a2 expansion coefficients (see ref 1)
    !===============================================================================================!

    subroutine get_a2matvec(d, u, n, A, b)

        class(nodal_type), intent(inout) :: d
        integer, intent(in)              :: u, n        ! Direction and node number
        real(dp), intent(out)            :: A(:,:)      ! GxG Materix (LHS)
        real(dp), intent(out)            :: b(:)        ! G vector (RHS)
    
        real(dp) :: S
        real(dp) :: Bf(d % ng)
        integer  :: g, h
    
        !Setup GxG matrix and G vector to obtain a2(g)
        Bf = 0.
        do g = 1, d % ng
            !update zeroth source
            if (d % calc_mode == 'transient') then
                if (u == 1) then
                    S = 0.25 * d % xdel(n)**2 / d % xs % D(n,g) * d % Sx(n,g)
                else if (u == 2) then
                    S = 0.25 * d % ydel(n)**2 / d % xs % D(n,g) * d % Sy(n,g)
                else
                    S = 0.25 * d % zdel(n)**2 / d % xs % D(n,g) * d % Sz(n,g)
                end if
            else
                if (u == 1) then
                    S = 0.25 * d % xdel(n)**2 / d % xs % D(n,g) * (d % Sx(n,g) - d % exsrc(n,g))
                else if (u == 2) then
                    S = 0.25 * d % ydel(n)**2 / d % xs % D(n,g) * (d % Sy(n,g) - d % exsrc(n,g))
                else
                    S = 0.25 * d % zdel(n)**2 / d % xs % D(n,g) * (d % Sz(n,g) - d % exsrc(n,g))
                end if
            end if
    
    
            ! Create matrix A (LHS)
            do h = 1, d % ng
                if (h == g) then
                    A(g,g) = d % Bcp(g,h) * d % Ep(g) + 3.
                else
                    A(g,h) = d % Bcp(g,h) * d % Ep(g)
                end if
                Bf(g) = Bf(g) + d % Bcp(g,h) * d % flux(n,h)
            end do
    
            ! Get second moment transverse leakage
            call TLUpd2(d, u, n, g, d % Lm2(g))
      
            ! RHS
            b(g) = Bf(g) - d % Ep(g) * d % Lm2(g) + S
        end do

    end subroutine

    !===============================================================================================!
    ! To calaculate A,B,E,F,G,H parameters used to calculate matrix elements for
    ! semi-analytic nodal update
    !===============================================================================================!

    subroutine get_ABEFGH (d, n, u, A, B, E, F, G, H)

        class(nodal_type), intent(inout) :: d
        integer, intent(in)              :: u, n
        real(dp), intent(out)            :: A(:), B(:)
        real(dp), intent(out)            :: E(:), F(:)
        real(dp), intent(out)            :: G(:), H(:)

        real(dp) :: dn, alpha, alpha2
        real(dp) :: m1s, m0c, m2c
        integer  :: gg

        if (u == 1) then
            dn = d % xdel(n)
        else if (u == 2) then
            dn = d % ydel(n)
        else
            dn = d % zdel(n)
        end if

        do gg = 1, d % ng
            !Calculate alpha and othe parameters to calculate A,B,E,F,G,H
            alpha  = 0.5 * sqrt(d % xs % sigr(n, gg) / d % xs % D(n, gg)) * dn
            alpha2 = alpha**2
            m0c    = sinh(alpha) / alpha
            m1s    = 3. * (cosh(alpha)/alpha - sinh(alpha)/alpha2)
            m2c    = 5. * (sinh(alpha)/alpha - 3.*cosh(alpha)/alpha2 &
            + 3.*sinh(alpha)/alpha**3)
        
            A(gg) = (sinh(alpha) - m1s) / (alpha2 * m1s)
            B(gg) = (cosh(alpha) - m0c - m2c) / (alpha2 * m2c)
            E(gg) = (m0c/m2c - 3./alpha2)
            F(gg) = (alpha*cosh(alpha) - m1s) / (alpha2 * m1s)
            G(gg) = (alpha*sinh(alpha) - 3.*m2c) / (cosh(alpha) - m0c - m2c)
            H(gg) = (alpha*cosh(alpha) - m1s) / (sinh(alpha) - m1s)
        end do

    end subroutine

    !===============================================================================================!
    ! To calculate Buckling for node n and direction u
    !===============================================================================================!

    function get_buckling(d, u, n) result (B)

        class(nodal_type), target, intent(inout) :: d
        integer, intent(in)                      :: u,n                 ! direction and node number
        real(dp)                                 :: B(d % ng, d % ng)   ! Buckling B2
    
        real(dp) :: dum, dn
        integer :: g, h
        real(dp) :: Ke
    
        Ke = d % Ke
        
        if (u == 1) then
            dn = d % xdel(n)
        else if (u == 2) then
            dn = d % ydel(n)
        else
            dn = d % zdel(n)
        end if
    
        if (d % calc_mode == 'static') then
            do g = 1, d % ng
                do h = 1, d % ng
                    if (g == h) then
                        dum = d % xs % sigr(n, g) - d % xs % chi(n, g) * d % xs % nuf(n, h) / Ke
                    else
                        dum = -d % xs % sigs(n,h,g) - d % xs % chi(n, g) * d % xs % nuf(n, h) / Ke
                    end if
                    B(g, h) = 0.25 * dn**2 / d % xs % D(n, g) * dum
                end do
            end do
        else if (d % calc_mode == 'transient') then
            stop "transient not yet implemented"
        else
            do g = 1, d % ng
                do h = 1, d % ng
                    if (g == h) then
                        dum = d % xs % sigr(n,g) - d % xs % chi(n,g) * d % xs % nuf(n,h) / Ke
                    else
                        dum = -d % xs % sigs(n,g,h) - d % xs % chi(n,h) * d % xs % nuf(n,g) / Ke
                    end if
                    B(g, h) = 0.25 * dn**2 / d % xs % D(n,g) * dum
                end do
            end do
        end if
    
    end function

    !===============================================================================================!
    ! To calculate source for the nodal update
    !===============================================================================================!

    subroutine get_source(d)

        class(nodal_type), target, intent(inout) :: d
        
        real(dp) :: L1(d % nnod, d % ng)
        real(dp) :: L2(d % nnod, d % ng)
        real(dp) :: L3(d % nnod, d % ng)
        integer  :: g, n
        
        call TLUpd0(d, L1, L2, L3)

        do g = 1, d % ng
            do n = 1, d % nnod
                if (d % calc_mode == 'transient') then
                    d % Sx(n,g) =  L2(n,g) + L3(n,g) - d % exsrc(n,g)
                    d % Sy(n,g) =  L1(n,g) + L3(n,g) - d % exsrc(n,g)
                    d % Sz(n,g) =  L1(n,g) + L2(n,g) - d % exsrc(n,g)
                else
                    d % Sx(n,g) =  L2(n,g) + L3(n,g)
                    d % Sy(n,g) =  L1(n,g) + L3(n,g)
                    d % Sz(n,g) =  L1(n,g) + L2(n,g)
                end if
            end do
        end do

    end subroutine

    !===============================================================================================!
    ! To update zeroth transverse leakages (J_plus - J_minus) for group g and nod n
    !===============================================================================================!

    subroutine TLUpd0(d, L1, L2, L3)
        
        class(nodal_type), target, intent(inout) :: d
        real(dp), intent(out)                    :: L1(d % nnod, d % ng)
        real(dp), intent(out)                    :: L2(d % nnod, d % ng)
        real(dp), intent(out)                    :: L3(d % nnod, d % ng)

        type(node_rect), pointer :: s
        real(dp) :: jp, jm
        integer  :: p, m
        integer  :: n, g

        do g = 1, d % ng
            do n = 1, d % nnod
                
                s => d % node(n)
        
                ! x-direction zeroth transverse leakage
                ! west
                if (.not. associated(s % east)) then
                    if (xeast == reflective) then
                        jp = 0.0
                    else
                        jp =    s % sf(1) % dt(g) * d % flux(n, g) - &
                                s % sf(1) % dh(g) * d % flux(n, g)
                    endif
                else
                    p = s % east % n
                    jp =   -s % sf(1) % dt(g) * (d % flux(p, g) - d % flux(n, g)) - &
                            s % sf(1) % dh(g) * (d % flux(p, g) + d % flux(n, g))
                end if
                ! east
                if (.not. associated(s % west)) then
                    if (xwest == reflective) then
                        jm = 0.0
                    else
                        jm =   -s % sf(2) % dt(g) * d % flux(n, g) - &
                                s % sf(2) % dh(g) * d % flux(n, g)
                    endif
                else
                    m = s % west % n
                    jm =   -s % sf(2) % dt(g) * (d % flux(n, g) - d % flux(m, g)) - &
                            s % sf(2) % dh(g) * (d % flux(n, g) + d % flux(m, g))
                end if
        
                L1(n, g) = (jp - jm)  / d % xdel(n)
        
                ! y-direction zeroth transverse leakage
                ! north
                if (.not. associated(s % north)) then
                    if (ynorth == reflective) then
                        jp = 0.0
                    else
                        jp =    s % sf(3) % dt(g) * d % flux(n, g) - &
                                s % sf(3) % dh(g) * d % flux(n, g)
                    endif
                else
                    p = s % north % n
                    jp =   -s % sf(3) % dt(g) * (d % flux(p, g) - d % flux(n, g)) - &
                            s % sf(3) % dh(g) * (d % flux(p, g) + d % flux(n, g))
                end if
                ! south
                if (.not. associated(s % south)) then
                    if (ysouth == reflective) then
                        jm = 0.0
                    else
                        jm =   -s % sf(4) % dt(g) * d % flux(n, g) - &
                                s % sf(4) % dh(g) * d % flux(n, g)
                    endif
                else
                    m = s % south % n
                    jm =   -s % sf(4) % dt(g) * (d % flux(n, g) - d % flux(m, g)) - &
                            s % sf(4) % dh(g) * (d % flux(n, g) + d % flux(m, g))
                end if
        
                L2(n, g) = (jp - jm)  / d % ydel(n)
        
                ! z-direction zeroth transverse leakage
                ! top
                if (.not. associated(s % top)) then
                    if (ztop == reflective) then
                        jp = 0.0
                    else
                        jp =    s % sf(5) % dt(g) * d % flux(n, g) - &
                                s % sf(5) % dh(g) * d % flux(n, g)
                    endif
                else
                    p = s % top % n
                    jp =   -s % sf(5) % dt(g) * (d % flux(p, g) - d % flux(n, g)) - &
                            s % sf(5) % dh(g) * (d % flux(p, g) + d % flux(n, g))
                end if
                ! bottom
                if (.not. associated(s % bottom)) then
                    if (zbott == reflective) then
                        jm = 0.0
                    else
                        jm =   -s % sf(6) % dt(g) * d % flux(n, g) - &
                                s % sf(6) % dh(g) * d % flux(n, g)
                    endif
                else
                    m = s % bottom % n
                    jm =   -s % sf(6) % dt(g) * (d % flux(n, g) - d % flux(m, g)) - &
                            s % sf(6) % dh(g) * (d % flux(n, g) + d % flux(m, g))
                end if
        
                L3(n, g) = (jp - jm) / d % zdel(n)
        
            end do
        end do

    end subroutine

    !===============================================================================================!
    ! To calaculate transverse leakage first moments
    !===============================================================================================!

    subroutine TLUpd1(d, u, n, g, Lmom1)

        class(nodal_type), intent(inout) :: d
        integer, intent(in) :: u, n, g
        real(dp), intent(out) :: Lmom1

        type(node_rect), pointer :: s
        real(dp)                 :: tm, tp
        real(dp)                 :: p1m, p2m, p1p, p2p, hp
        integer                  :: p, m

        s => d % node(n)

        if (u == 1) then
            ! x-direction first moment transverse leakage
            if (.not. associated(s % east)) then
                m = s % west % n
                if (xeast == reflective) then
                    tm = d % xdel(m) / d % xdel(n)
                    tp = 1.
                    p1m = tm + 1.
                    p1p = tp + 1.; p2p = 2. * tp + 1.
                    hp = 2. * p1m * p1p * (tm + tp + 1.)
                    Lmom1 = (p1p * p2p * (d % Sx(n,g) - d % Sx(m,g))) / hp
                else
                    tm = d % xdel(m) / d % xdel(n)
                    p1m = tm + 1.
                    Lmom1 = (d % Sx(n,g) - d % Sx(m,g)) / p1m
                end if
            elseif (.not. associated(s % west)) then
                p = s % east % n
                if (xwest == reflective) then
                    tm = 1.
                    tp = d % xdel(p) / d % xdel(n)
                    p1m = tm + 1.; p2m = 2. * tm + 1.
                    p1p = tp + 1.
                    hp = 2. * p1m * p1p * (tm + tp + 1.)
                    Lmom1 = (p1m * p2m * (d % Sx(p,g) - d % Sx(n,g))) / hp
                else
                    tp = d % xdel(p) / d % xdel(n)
                    p1p = tp + 1.
                    Lmom1 = (d % Sx(p,g) - d % Sx(n,g)) / p1p
                end if
            else
                m = s % west % n
                p = s % east % n
                tm = d % xdel(m) / d % xdel(n)
                tp = d % xdel(p) / d % xdel(n)
                p1m = tm + 1.; p2m = 2. * tm + 1.
                p1p = tp + 1.; p2p = 2. * tp + 1.
                hp = 2.* p1m * p1p * (tm + tp + 1.)
                Lmom1 = (p1m * p2m * (d % Sx(p,g) - d % Sx(n,g)) &
                      +  p1p * p2p * (d % Sx(n,g) - d % Sx(m,g))) / hp
            endif

            Lmom1 = 0.25 * d % xdel(n)**2 / d % xs % D(n,g) * Lmom1

        elseif (u == 2) then
            ! y-direction first moment transverse leakage
            if (.not. associated(s % north)) then
                m = s % south % n
                if (ynorth == reflective) then
                    tm = d % ydel(m) / d % ydel(n)
                    tp = 1.
                    p1m = tm + 1.
                    p1p = tp + 1.; p2p = 2. * tp + 1.
                    hp = 2. * p1m * p1p * (tm + tp + 1.)
                    Lmom1 = (p1p * p2p * (d % Sy(n,g) - d % Sy(m,g))) / hp
                else
                    tm = d % ydel(m) / d % ydel(n)
                    p1m = tm + 1.
                    Lmom1 = (d % Sy(n,g) - d % Sy(m,g)) / p1m
                end if
            elseif (.not. associated(s % south)) then
                p = s % north % n
                if (ysouth == reflective) then
                    tm = 1.
                    tp = d % ydel(p) / d % ydel(n)
                    p1m = tm + 1.; p2m = 2. * tm + 1.
                    p1p = tp + 1.
                    hp = 2. * p1m * p1p * (tm + tp + 1.)
                    Lmom1 = (p1m * p2m * (d % Sy(p,g) - d % Sy(n,g))) / hp
                else
                    tp = d % ydel(p) / d % ydel(n)
                    p1p = tp + 1.
                    Lmom1 = (d % Sy(p,g) - d % Sy(n,g)) / p1p
                end if
            else
                m = s % south % n
                p = s % north % n
                tm = d % ydel(m) / d % ydel(n)
                tp = d % ydel(p) / d % ydel(n)
                p1m = tm + 1.; p2m = 2. * tm + 1.
                p1p = tp + 1.; p2p = 2. * tp + 1.
                hp = 2.* p1m * p1p * (tm + tp + 1.)
                Lmom1 = (p1m * p2m * (d % Sy(p,g) - d % Sy(n,g)) &
                      +  p1p * p2p * (d % Sy(n,g) - d % Sy(m,g))) / hp
            endif

            Lmom1 = 0.25 * d % ydel(n)**2 / d % xs % D(n,g) * Lmom1

        else
            ! z-direction first moment transverse leakage
            if (.not. associated(s % top)) then
                m = s % bottom % n
                if (ztop == reflective) then
                    tm = d % zdel(m) / d % zdel(n)
                    tp = 1.
                    p1m = tm + 1.
                    p1p = tp + 1.; p2p = 2. * tp + 1.
                    hp = 2. * p1m * p1p * (tm + tp + 1.)
                    Lmom1 = (p1p * p2p * (d % Sz(n,g) - d % Sz(m,g))) / hp
                else
                    tm = d % zdel(m) / d % zdel(n)
                    p1m = tm + 1.
                    Lmom1 = (d % Sz(n,g) - d % Sz(m,g)) / p1m
                end if
            elseif (.not. associated(s % bottom)) then
                p = s % top % n
                if (zbott == reflective) then
                    tm = 1.
                    tp = d % zdel(p) / d % zdel(n)
                    p1m = tm + 1.; p2m = 2. * tm + 1.
                    p1p = tp + 1.
                    hp = 2. * p1m * p1p * (tm + tp + 1.)
                    Lmom1 = (p1m * p2m * (d % Sz(p,g) - d % Sz(n,g))) / hp
                else
                    tp = d % zdel(p) / d % zdel(n)
                    p1p = tp + 1.
                    Lmom1 = (d % Sz(p,g) - d % Sz(n,g)) / p1p
                end if
            else
                m = s % bottom % n
                p = s % top % n
                tm = d % zdel(m) / d % zdel(n)
                tp = d % zdel(p) / d % zdel(n)
                p1m = tm + 1.; p2m = 2. * tm + 1.
                p1p = tp + 1.; p2p = 2. * tp + 1.
                hp = 2.* p1m * p1p * (tm + tp + 1.)
                Lmom1 = (p1m * p2m * (d % Sz(p,g) - d % Sz(n,g)) &
                      +  p1p * p2p * (d % Sz(n,g) - d % Sz(m,g))) / hp
            endif

            Lmom1 = 0.25 * d % zdel(n)**2 / d % xs % D(n,g) * Lmom1

        end if

    end subroutine

    !===============================================================================================!
    ! To calaculate transverse leakage first moments
    !===============================================================================================!

    subroutine TLUpd2(d, u, n, g, Lmom2)

        class(nodal_type), intent(inout) :: d
        integer, intent(in) :: u, n, g
        real(dp), intent(out) :: Lmom2

        type(node_rect), pointer :: s
        real(dp)                 :: tm, tp
        real(dp)                 :: p1m, p1p, hp
        integer                  :: p, m

        s => d % node(n)

        if (u == 1) then
            ! x-direction second moment transverse leakage
            if (.not. associated(s % east)) then
                m = s % west % n
                if (xeast == reflective) then
                    tm = d % xdel(m) / d % xdel(n)
                    tp = 1.
                    p1m = tm + 1.
                    p1p = tp + 1.
                    hp = 2. * p1m * p1p* (tm + tp + 1.)
                    Lmom2 = (p1p * (d % Sx(m,g) - d % Sx(n,g))) / hp
                else
                    Lmom2 = 0.
                end if
            elseif (.not. associated(s % west)) then
                p = s % east % n
                if (xwest == reflective) then
                    tm = 1.
                    tp = d % xdel(p) / d % xdel(p)
                    p1m = tm + 1.
                    p1p = tp + 1.
                    hp = 2. * p1m * p1p* (tm + tp + 1.)
                    Lmom2 = (p1m * (d % Sx(p,g) - d % Sx(n,g))) / hp
                else
                    Lmom2 = 0.
                end if
            else
                m = s % west % n
                p = s % east % n
                tm = d % xdel(m) / d % xdel(n)
                tp = d % xdel(p) / d % xdel(n)
                p1m = tm + 1.
                p1p = tp + 1.
                hp = 2. * p1m * p1p * (tm + tp + 1.)
                Lmom2 = (p1m * (d % Sx(p,g) - d % Sx(n,g)) + &
                         p1p * (d % Sx(m,g) - d % Sx(n,g))) / hp
            endif

            Lmom2 = 0.25 * d % xdel(n)**2 / d % xs % D(n,g) * Lmom2

        elseif (u == 2) then
            ! y-direction second moment transverse leakage
            if (.not. associated(s % north)) then
                m = s % south % n
                if (ynorth == reflective) then
                    tm = d % ydel(m) / d % ydel(n)
                    tp = 1.
                    p1m = tm + 1.
                    p1p = tp + 1.
                    hp = 2. * p1m * p1p* (tm + tp + 1.)
                    Lmom2 = (p1p * (d % Sy(m,g) - d % Sy(n,g))) / hp
                else
                    Lmom2 = 0.
                end if
            elseif (.not. associated(s % south)) then
                p = s % north % n
                if (ysouth == reflective) then
                    tm = 1.
                    tp = d % ydel(p) / d % ydel(p)
                    p1m = tm + 1.
                    p1p = tp + 1.
                    hp = 2. * p1m * p1p* (tm + tp + 1.)
                    Lmom2 = (p1m * (d % Sy(p,g) - d % Sy(n,g))) / hp
                else
                    Lmom2 = 0.
                end if
            else
                m = s % south % n
                p = s % north % n
                tm = d % ydel(m) / d % ydel(n)
                tp = d % ydel(p) / d % ydel(n)
                p1m = tm + 1.
                p1p = tp + 1.
                hp = 2. * p1m * p1p * (tm + tp + 1.)
                Lmom2 = (p1m * (d % Sy(p,g) - d % Sy(n,g)) + &
                         p1p * (d % Sy(m,g) - d % Sy(n,g))) / hp
            endif

            Lmom2 = 0.25 * d % ydel(n)**2 / d % xs % D(n,g) * Lmom2

        else
            ! z-direction second moment transverse leakage
            if (.not. associated(s % top)) then
                m = s % bottom % n
                if (ztop == reflective) then
                    tm = d % zdel(m) / d % zdel(n)
                    tp = 1.
                    p1m = tm + 1.
                    p1p = tp + 1.
                    hp = 2. * p1m * p1p* (tm + tp + 1.)
                    Lmom2 = (p1p * (d % Sz(m,g) - d % Sz(n,g))) / hp
                else
                    Lmom2 = 0.
                end if
            elseif (.not. associated(s % bottom)) then
                p = s % top % n
                if (zbott == reflective) then
                    tm = 1.
                    tp = d % zdel(p) / d % zdel(p)
                    p1m = tm + 1.
                    p1p = tp + 1.
                    hp = 2. * p1m * p1p* (tm + tp + 1.)
                    Lmom2 = (p1m * (d % Sz(p,g) - d % Sz(n,g))) / hp
                else
                    Lmom2 = 0.
                end if
            else
                m = s % bottom % n
                p = s % top % n
                tm = d % zdel(m) / d % zdel(n)
                tp = d % zdel(p) / d % zdel(n)
                p1m = tm + 1.
                p1p = tp + 1.
                hp = 2. * p1m * p1p * (tm + tp + 1.)
                Lmom2 = (p1m * (d % Sz(p,g) - d % Sz(n,g)) + &
                         p1p * (d % Sz(m,g) - d % Sz(n,g))) / hp
            endif

            Lmom2 = 0.25 * d % zdel(n)**2 / d % xs % D(n,g) * Lmom2

        end if

    end subroutine

    
end module


    