module nodal

  use sdata, only: dp
  implicit none
  save

  integer :: cmode
  real(dp), dimension(:,:), allocatable  :: Bcn, Bcp          ! Buckling for left and right nodes
  real(dp), dimension(:), allocatable    :: An, Bn, En, Fn, Gn, Hn
  real(dp), dimension(:), allocatable    :: Ap, Bp, Ep, Fp, Gp, Hp
  real(dp), dimension(:), allocatable    :: Lm2                ! Second order source
  real(dp), dimension(:,:), allocatable :: S1, S2, S3          ! zeroth order source in x, y, and z directions

contains

  !****************************************************************************!

  subroutine nodal_update(cal_mode)

    !Purpose: to calculate flux expansion coefficients using SANM


    use sdata, only: ng, xyz, ystag, xstag, nyy, nzz, nxx, &
    a1n, a2n, a3n, a4n, a1p, a2p, a3p, a4p, ndmax, &
    Ln1, Lp1, xeast, xwest, ysouth, ynorth, zbott, ztop

    implicit none

    integer, intent(in)  :: cal_mode  !Calc. mode -> 0=adjoint, 1=forward

    integer  :: n, p
    integer  :: i, j, k
    logical  :: first =.true.

    if (first) then
      allocate(a1n(ng), a2n(ng), a3n(ng), a4n(ng))
      allocate(a1p(ng), a2p(ng), a3p(ng), a4p(ng))
      first = .false.
    end if

    allocate(Ln1(ng), Lp1(ng))
    allocate(Bcn(ng,ng), Bcp(ng,ng))
    allocate(An(ng), Bn(ng), En(ng), Fn(ng), Gn(ng), Hn(ng))
    allocate(Ap(ng), Bp(ng), Ep(ng), Fp(ng), Gp(ng), Hp(ng))
    allocate(Lm2(ng))

    cmode = cal_mode

    call get_source()


    !Node sweeps in x-direction
    do k = 1, nzz
      do j = 1, nyy
        p = xyz(ystag(j)%smin,j,k)
        Bcp = get_B(1,p)
        call get_ABEFGH (p,1,Ap,Bp,Ep,Fp,Gp,Hp)
        call get_coefs_first(xwest, 1, p)
        call nodal_coup_upd(1, a1p, a2p, a3p, a4p, p = p)
        do i = ystag(j)%smin, ystag(j)%smax-1
          n = xyz(i,j,k); p = xyz(i+1,j,k)
          Bcn = Bcp
          An = Ap; Bn = Bp; En = Ep; Fn = Fp; Gn = Gp; Hn = Hp
          Bcp = get_B(1,p)
          call get_ABEFGH (p,1,Ap,Bp,Ep,Fp,Gp,Hp)
          call get_coefs(1, n, p)
          call nodal_coup_upd(1, a1n, a2n, a3n, a4n, n, p)
        end do
        n = xyz(ystag(j)%smax,j,k)
        Bcn = Bcp
        An = Ap; Bn = Bp; En = Ep; Fn = Fp; Gn = Gp; Hn = Hp
        call get_coefs_last(xeast, 1, n)
        call nodal_coup_upd(1, a1n, a2n, a3n, a4n, n = n)
      end do
    end do

    !Node sweeps in y-direction
    do k = 1, nzz
      do i = 1, nxx
        p = xyz(i,xstag(i)%smin,k)
        Bcp = get_B(2,p)
        call get_ABEFGH (p,2,Ap,Bp,Ep,Fp,Gp,Hp)
        call get_coefs_first(ysouth, 2, p)
        call nodal_coup_upd(2, a1p, a2p, a3p, a4p, p = p)
        do j = xstag(i)%smin, xstag(i)%smax-1
          n = xyz(i,j,k); p = xyz(i,j+1,k)
          Bcn = Bcp
          An = Ap; Bn = Bp; En = Ep; Fn = Fp; Gn = Gp; Hn = Hp
          Bcp = get_B(2,p)
          call get_ABEFGH (p,2,Ap,Bp,Ep,Fp,Gp,Hp)
          call get_coefs(2, n, p)
          call nodal_coup_upd(2, a1n, a2n, a3n, a4n, n, p)
        end do
        n = xyz(i,xstag(i)%smax,k)
        Bcn = Bcp
        An = Ap; Bn = Bp; En = Ep; Fn = Fp; Gn = Gp; Hn = Hp
        call get_coefs_last(ynorth, 2, n)
        call nodal_coup_upd(2, a1n, a2n, a3n, a4n, n = n)
      end do
    end do

    !Node sweeps in z-direction
    do j = 1, nyy
      do i = ystag(j)%smin, ystag(j)%smax
        p = xyz(i,j,1)
        Bcp = get_B(3,p)
        call get_ABEFGH (p,3,Ap,Bp,Ep,Fp,Gp,Hp)
        call get_coefs_first(zbott, 3, p)
        call nodal_coup_upd(3, a1p, a2p, a3p, a4p, p = p)
        do k = 1, nzz-1
          n = xyz(i,j,k); p = xyz(i,j,k+1)
          Bcn = Bcp
          An = Ap; Bn = Bp; En = Ep; Fn = Fp; Gn = Gp; Hn = Hp
          Bcp = get_B(3,p)
          call get_ABEFGH (p,3,Ap,Bp,Ep,Fp,Gp,Hp)
          call get_coefs(3, n, p)
          call nodal_coup_upd(3, a1n, a2n, a3n, a4n, n, p)
        end do
        n = xyz(i,j,nzz)
        Bcn = Bcp
        An = Ap; Bn = Bp; En = Ep; Fn = Fp; Gn = Gp; Hn = Hp
        call get_coefs_last(ztop, 3, n)
        call nodal_coup_upd(3, a1n, a2n, a3n, a4n, n = n)
      end do
    end do

    deallocate(Ln1, Lp1,  Bcp, Bcn, Lm2)
    deallocate(An, Bn, En, Fn, Gn, Hn)
    deallocate(Ap, Bp, Ep, Fp, Gp, Hp)

    if (ndmax > 1.e3) then
      write(*,*)
      write(*,1236) ndmax
      write(*,*)"The two-node nonlinear iteration seems not stable."
      write(*,*)"Try to reduce time step size for transient..."
      write(*,*)"and/or change iteration control using %ITER card, ..."
      write(*,*)"perhaps by making nodal update less frequent."
      write(*,*)"Or increase number inner iteration per outer iteration."
      write(*,*)"If this error persists, contact me at makrus.imron@gmail.com"
      write(*,*)"Thank you for using ADPRES"
      stop
    end if

    1236 format(" Error: Max. change in nodal coupling coefficient = ", F10.1)


  end subroutine nodal_update

  !****************************************************************************!

  subroutine nodal_update_pnm(cal_mode)

    !Purpose: to calculate flux expansion coefficients using PNM


    use sdata, only: ng, xyz, ystag, xstag, nyy, nzz, nxx, &
    a1n, a2n, a3n, a4n, a1p, a2p, a3p, a4p, ndmax, &
    Ln1, Lp1, xeast, xwest, ysouth, ynorth, zbott, ztop

    implicit none

    integer, intent(in)  :: cal_mode  !Calc. mode -> 0=adjoint, 1=forward

    integer  :: n, p
    integer  :: i, j, k
    logical  :: first =.true.

    if (first) then
      allocate(a1n(ng), a2n(ng), a3n(ng), a4n(ng))
      allocate(a1p(ng), a2p(ng), a3p(ng), a4p(ng))
      first = .false.
    end if

    allocate(Ln1(ng), Lp1(ng))
    allocate(Bcn(ng,ng), Bcp(ng,ng))
    allocate(An(ng), Bn(ng), En(ng), Fn(ng), Gn(ng), Hn(ng))
    allocate(Ap(ng), Bp(ng), Ep(ng), Fp(ng), Gp(ng), Hp(ng))
    allocate(Lm2(ng))

    An = 1._dp / 15._dp; Bn = 1._dp / 35._dp; En = 2._dp / 7._dp
    Fn = 2._dp / 5._dp; Gn = 10._dp; Hn = 6._dp
    Ap = An; Bp = Bn; Ep = En; Fp = Fn; Gp = Gn; Hp = Hn

    cmode = cal_mode

    call get_source()


    !Node sweeps in x-direction
    do k = 1, nzz
      do j = 1, nyy
        p = xyz(ystag(j)%smin,j,k)
        Bcp = get_B(1,p)
        call get_coefs_first(xwest, 1, p)
        call nodal_coup_upd(1, a1p, a2p, a3p, a4p, p = p)
        do i = ystag(j)%smin, ystag(j)%smax-1
          n = xyz(i,j,k); p = xyz(i+1,j,k)
          Bcn = Bcp
          An = Ap; Bn = Bp; En = Ep; Fn = Fp; Gn = Gp; Hn = Hp
          Bcp = get_B(1,p)
          call get_coefs(1, n, p)
          call nodal_coup_upd(1, a1n, a2n, a3n, a4n, n, p)
        end do
        n = xyz(ystag(j)%smax,j,k)
        Bcn = Bcp
        An = Ap; Bn = Bp; En = Ep; Fn = Fp; Gn = Gp; Hn = Hp
        call get_coefs_last(xeast, 1, n)
        call nodal_coup_upd(1, a1n, a2n, a3n, a4n, n = n)
      end do
    end do

    !Node sweeps in y-direction
    do k = 1, nzz
      do i = 1, nxx
        p = xyz(i,xstag(i)%smin,k)
        Bcp = get_B(2,p)
        call get_coefs_first(ysouth, 2, p)
        call nodal_coup_upd(2, a1p, a2p, a3p, a4p, p = p)
        do j = xstag(i)%smin, xstag(i)%smax-1
          n = xyz(i,j,k); p = xyz(i,j+1,k)
          Bcn = Bcp
          An = Ap; Bn = Bp; En = Ep; Fn = Fp; Gn = Gp; Hn = Hp
          Bcp = get_B(2,p)
          call get_coefs(2, n, p)
          call nodal_coup_upd(2, a1n, a2n, a3n, a4n, n, p)
        end do
        n = xyz(i,xstag(i)%smax,k)
        Bcn = Bcp
        An = Ap; Bn = Bp; En = Ep; Fn = Fp; Gn = Gp; Hn = Hp
        call get_coefs_last(ynorth, 2, n)
        call nodal_coup_upd(2, a1n, a2n, a3n, a4n, n = n)
      end do
    end do

    !Node sweeps in z-direction
    do j = 1, nyy
      do i = ystag(j)%smin, ystag(j)%smax
        p = xyz(i,j,1)
        Bcp = get_B(3,p)
        call get_coefs_first(zbott, 3, p)
        call nodal_coup_upd(3, a1p, a2p, a3p, a4p, p = p)
        do k = 1, nzz-1
          n = xyz(i,j,k); p = xyz(i,j,k+1)
          Bcn = Bcp
          An = Ap; Bn = Bp; En = Ep; Fn = Fp; Gn = Gp; Hn = Hp
          Bcp = get_B(3,p)
          call get_coefs(3, n, p)
          call nodal_coup_upd(3, a1n, a2n, a3n, a4n, n, p)
        end do
        n = xyz(i,j,nzz)
        Bcn = Bcp
        An = Ap; Bn = Bp; En = Ep; Fn = Fp; Gn = Gp; Hn = Hp
        call get_coefs_last(ztop, 3, n)
        call nodal_coup_upd(3, a1n, a2n, a3n, a4n, n = n)
      end do
    end do

    deallocate(Ln1, Lp1,  Bcp, Bcn, Lm2)
    deallocate(An, Bn, En, Fn, Gn, Hn)
    deallocate(Ap, Bp, Ep, Fp, Gp, Hp)

    if (ndmax > 1.e3) then
      write(*,*)
      write(*,1236) ndmax
      write(*,*)"The two-node nonlinear iteration seems not stable."
      write(*,*)"Try to reduce time step size for transient..."
      write(*,*)"and/or change iteration control using %ITER card, ..."
      write(*,*)"perhaps by making nodal update less frequent."
      write(*,*)"Or increase number inner iteration per outer iteration."
      write(*,*)"If this error persists, contact me at makrus.imron@gmail.com"
      write(*,*)"Thank you for using ADPRES"
      stop
    end if

    1236 format(" Error: Max. change in nodal coupling coefficient = ", F10.1)


  end subroutine nodal_update_pnm

  !****************************************************************************!

  subroutine nodal_coup_upd(u, a1, a2, a3, a4, n, p)

    !Purpose: to update nodal coupling coefficients


    use sdata, only: ng, ix, iy, iz, D, xdel, ydel, zdel, nod, f0, &
    ndmax, im, jm, km

    implicit none

    integer, intent(in)                :: u
    integer, intent(in), optional      :: n, p
    real(dp), dimension(:), intent(in) :: a1, a2, a3, a4

    integer    :: g, sf
    real(dp)   :: dh, jp, nder, ndpr

    if (u == 1) then
      if (present(n)) then
        dh = xdel(ix(n))
      else
        dh = xdel(ix(p))
      end if
      sf = 1
    else if (u == 2) then
      if (present(n)) then
        dh = ydel(iy(n))
      else
        dh = ydel(iy(p))
      end if
      sf = 3
    else
      if (present(n)) then
        dh = zdel(iz(n))
      else
        dh = zdel(iz(p))
      end if
      sf = 5
    end if

    if (present(n) .and. present(p)) then
      do g = 1, ng
        jp = -2._dp*D(n,g)/dh * (a1(g)+3._dp*a2(g)+Hn(g)*a3(g)+Gn(g)*a4(g))

        ! Update nodal coupling
        ndpr = nod(n,g)%dn(sf)
        nod(n,g)%dn(sf) = (nod(n,g)%df(sf) * (f0(n,g) - f0(p,g)) - jp) &
        / (f0(n,g) + f0(p,g))
        nod(p,g)%dn(sf+1) = nod(n,g)%dn(sf)

        ! Check max difference on new nodal coupling coefficients
        nder = ABS(nod(n,g)%dn(sf) - ndpr)
        if (nder > ndmax) then
          ndmax = nder
          im = ix(n); jm = iy(n); km = iz(n)
        end if
      end do
    else if (present(p)) then
      do g = 1, ng
        jp = -2._dp*D(p,g)/dh*(a1(g)-3._dp*a2(g)+Hp(g)*a3(g)-Gp(g)*a4(g))

        ! Update nodal coupling
        ndpr = nod(p,g)%dn(sf+1)
        nod(p,g)%dn(sf+1) = -(jp / f0(p,g) + nod(p,g)%df(sf+1))

        ! Check max difference on new nodal coupling coefficients
        nder = ABS(nod(p,g)%dn(sf+1) - ndpr)
        if (nder > ndmax) then
          ndmax = nder
          im = ix(p); jm = iy(p); km = iz(p)
        end if
      end do
    else
      do g = 1, ng
        jp = -2._dp*D(n,g)/dh*(a1(g)+3._dp*a2(g)+Hn(g)*a3(g)+Gn(g)*a4(g))

        ! Update nodal coupling
        ndpr = nod(n,g)%dn(sf)
        nod(n,g)%dn(sf) = -(jp / f0(n,g) - nod(n,g)%df(sf))

        ! Check max difference on new nodal coupling coefficients
        nder = ABS(nod(n,g)%dn(sf) - ndpr)
        if (nder > ndmax) then
          ndmax = nder
          im = ix(n); jm = iy(n); km = iz(n)
        end if
      end do
    end if

  end subroutine nodal_coup_upd

  !****************************************************************************!

  subroutine get_coefs_first(bc, u, p)

    !Purpose: to calculate flux expansion coefficients


    use sdata, only: ng, a1p, a2p, a3p, a4p, Lp1

    implicit none

    integer, intent(in)                  :: bc
    integer, intent(in)                  :: u, p

    real(dp), dimension(ng,ng)     :: A           ! GxG Matrix
    real(dp), dimension(ng)        :: b           ! G vector

    !Setup GxG matrix and G vector to obtain a2(g) for node p
    call get_a2matvec(u,p,A,b)
    !calculate a2 expansion coefficients4
    a2p = LU_solve(p,ng,A,b)
    !calculate a4 expansion coefficients
    a4p = get_a4(a2p)

    !Setup GxG matrix and G vector to obtain a1(g) for left right node
    call get_a1matvec_first(bc, u, p, a2p, a4p, A, b)

    !calculate a1 expansion coefficients
    a1p = LU_solve(p,ng,A,b)

    !calculate a3 expansion coefficients
    a3p = get_a3(2,a1p,Lp1)

  end subroutine get_coefs_first

  !****************************************************************************!

  subroutine get_coefs_last(bc, u, n)

    !Purpose: to calculate flux expansion coefficients


    use sdata, only: ng, a1n, a2n, a3n, a4n, a2p, a4p, Ln1

    implicit none

    integer, intent(in)                  :: bc
    integer, intent(in)                  :: u, n

    real(dp), dimension(ng,ng)     :: A           ! GxG Matrix
    real(dp), dimension(ng)        :: b           ! G vector

    !get a2n and a4n expansion coefficients
    a2n = a2p
    a4n = a4p

    !Setup GxG matrix and G vector to obtain a1(g) for most right node
    call get_a1matvec_last(bc, u, n, a2n, a4n, A, b)

    !calculate a1 expansion coefficients
    a1n = LU_solve(n,ng,A,b)

    !calculate a3 expansion coefficients
    a3n = get_a3(1,a1n,Ln1)

  end subroutine get_coefs_last

  !****************************************************************************!

  subroutine get_coefs(u, n, p)

    !Purpose: to calculate flux expansion coefficients


    use sdata, only: ng, a1n, a2n, a3n, a4n, a2p, a4p, Ln1

    implicit none

    integer, intent(in)                  :: u, n, p

    real(dp), dimension(ng,ng)     :: A           ! GxG Matrix
    real(dp), dimension(ng)        :: b           ! G vector
    real(dp), dimension(2*ng,2*ng) :: R           ! 2Gx2G Matrix
    real(dp), dimension(2*ng)      :: s           ! 2G vector

    !get a2n and a4n expansion coefficients
    a2n = a2p
    a4n = a4p

    !Setup GxG matrix and G vector to obtain a2(g) for node p
    call get_a2matvec(u,p,A,b)
    !calculate a2 expansion coefficients
    a2p = LU_solve(p,ng,A,b)
    !calculate a4 expansion coefficients
    a4p = get_a4(a2p)

    !Setup 2Gx2G matrix and 2G vector to obtain a1(g) for node n
    call get_a1matvec(u, n, p, a2n, a4n, a2p, a4p, R, s)

    !calculate a1 expansion coefficients
    s = LU_solve(n,2*ng,R,s)
    a1n = s(1:ng)

    !calculate a3 expansion coefficients
    a3n = get_a3(1,a1n,Ln1)

  end subroutine get_coefs

  !****************************************************************************!

  subroutine get_a1matvec_first(bc, u, p, a2p, a4p, A, b)

    !Purpose: To get matrix vector to calculate a1 for most left node


    use sdata, only: ng, xdel, ydel, zdel, ix, iy, iz, f0, D, Lp1, dc

    implicit none

    integer, intent(in)                  :: bc
    integer, intent(in)                  :: u, p         ! Direction and node n umber
    real(dp), dimension(:), intent(in)   :: a2p, a4p     ! a2 and a4 expansion coefficients
    real(dp), dimension(:,:), intent(out):: A            ! 2Gx2G Materix
    real(dp), dimension(:), intent(out)  :: b            ! 2G vector

    real(dp)                   :: dn
    real(dp)                   :: Pp                       !dif. coef/dx
    integer                    :: g, h, sf

    ! define node size in direction u
    if (u == 1) then
      dn = xdel(ix(p))
      sf = 2
    else if (u == 2) then
      dn = ydel(iy(p))
      sf = 4
    else
      dn = zdel(iz(p))
      sf = 6
    end if

    do g = 1, ng
      !Calculate transverse leakage first moment
      call TLUpd1 (u,p,g,Lp1(g))
      Pp = 2._dp * D(p,g) / dn

      !Setup GxG matrix A and G vector B to obtain a2(g)
      if (bc == 2) then
        do h = 1, ng
          if (h == g) then
            A(g,g) = Pp * (Bcp(g,h)*Fp(g) + 1._dp)
          else
            A(g,h) = Pp*Bcp(g,h)*Fp(g)
          end if
        end do
        b(g) = Pp * (3._dp*a2p(g) + Gp(g)*a4p(g) - Fp(g)*Lp1(g))
      else if (bc == 1) then
        do h = 1, ng
          if (h == g) then
            A(g,g) = -dc(p,g,sf) * (1._dp + Ap(g)*Bcp(g,h)) &
                   - 2._dp*Pp*(Ap(g)*Bcp(g,h)*Hp(g)+1._dp)
          else
            A(g,h) = -dc(p,g,sf)*Ap(g)*Bcp(g,h) - 2._dp*Pp*Ap(g)*Bcp(g,h)*Hp(g)
          end if
        end do
        b(g) = 2._dp*Pp*(Ap(g)*Hp(g)*Lp1(g) - 3._dp*a2p(g) - Gp(g)*a4p(g)) &
             - dc(p,g,sf) * (a2p(g) + a4p(g) + f0(p,g) - Ap(g)*Lp1(g))
      else
        do h = 1, ng
          if (h == g) then
            A(g,g) = dc(p,g,sf) * (1._dp + Ap(g)*Bcp(g,h))
          else
            A(g,h) = dc(p,g,sf)*Ap(g)*Bcp(g,h)
          end if
        end do
        b(g) =  dc(p,g,sf) * (a2p(g) + a4p(g) + f0(p,g) - Ap(g) * Lp1(g))
      end if
    end do

  end subroutine get_a1matvec_first

  !****************************************************************************!

  subroutine get_a1matvec_last(bc, u, n, a2n, a4n, A, b)

    !Purpose: To get matrix vector to calculate a1 for most right node


    use sdata, only: ng, xdel, ydel, zdel, ix, iy, iz, f0, D, Lp1, Ln1, dc

    implicit none

    integer, intent(in)                  :: bc
    integer, intent(in)                  :: u, n         ! Direction and node number
    real(dp), dimension(:), intent(in)   :: a2n, a4n           ! a2 and a4 expansion coefficients
    real(dp), dimension(:,:), intent(out):: A           ! 2Gx2G Materix
    real(dp), dimension(:), intent(out)  :: b           ! 2G vector

    real(dp)                   :: dn
    real(dp)                   :: Pn                       !dif. coef/dx
    integer                    :: g, h, sf

    ! define node size in direction u
    if (u == 1) then
      dn = xdel(ix(n))
      sf = 1
    else if (u == 2) then
      dn = ydel(iy(n))
      sf = 3
    else
      dn = zdel(iz(n))
      sf = 5
    end if


    do g = 1, ng
      Ln1(g) = Lp1(g)
      Pn  = 2._dp * D(n,g) / dn

      !Setup 2Gx2G matrix A and 2G vector B to obtain a2(g)
      if (bc == 2) then
        do h = 1, ng
          if (h == g) then
            A(g,g) = -Pn * (Bcn(g,h)*Fn(g) + 1._dp)
          else
            A(g,h) = -Pn*Bcn(g,h)*Fn(g)
          end if
        end do
        b(g) = Pn * (3._dp*a2n(g) + Gn(g)*a4n(g) + Fn(g)*Ln1(g))
      else if (bc == 1) then
        do h = 1, ng
          if (h == g) then
            A(g,g) = dc(n,g,sf) * (1._dp + An(g)*Bcn(g,h)) &
                   + 2._dp*Pn*(An(g)*Bcn(g,h)*Hn(g)+1._dp)
          else
            A(g,h) = dc(n,g,sf)*An(g)*Bcn(g,h) + 2._dp*Pn*An(g)*Bcn(g,h)*Hn(g)
          end if
        end do
        b(g) = -2._dp*Pn*(An(g)*Hn(g)*Ln1(g) + 3._dp*a2n(g) + Gn(g)*a4n(g)) &
             - dc(n,g,sf) * (a2n(g) + a4n(g) + f0(n,g) + An(g)*Ln1(g))
      else
        do h = 1, ng
          if (h == g) then
            A(g,g) = dc(n,g,sf) * (1._dp + An(g)*Bcn(g,h))
          else
            A(g,h) = dc(n,g,sf)*An(g)*Bcn(g,h)
          end if
        end do
        b(g) = -dc(n,g,sf) * (a2n(g) + a4n(g) + f0(n,g) + An(g) * Ln1(g))
      end if
    end do


  end subroutine get_a1matvec_last

  !****************************************************************************!

  subroutine get_a1matvec(u, n, p, a2n, a4n, a2p, a4p, A,b)

    !Purpose: To setup 2Gx2G matrix and 2G vector to get a1 expansion coefficients


    use sdata, only: ng, xdel, ydel, zdel, ix, iy, iz, D, f0, Ln1, Lp1, dc

    implicit none

    integer, intent(in)                  :: u, n, p         ! Direction and node number
    real(dp), dimension(:), intent(in)   :: a2n, a4n           ! a2 and a4 expansion coefficients
    real(dp), dimension(:), intent(in)   :: a2p, a4p           ! a2 and a4 expansion coefficients
    real(dp), dimension(:,:), intent(out):: A           ! 2Gx2G Materix
    real(dp), dimension(:), intent(out)  :: b           ! 2G vector

    real(dp)                   :: hn, hp
    real(dp), dimension(ng)    :: Pn, Pp              !dif. coef/dx
    integer                    :: g, h, sf

    ! define node size in direction u
    if (u == 1) then
      hn = xdel(ix(n)); hp = xdel(ix(p))
      sf = 1
    else if (u == 2) then
      hn = ydel(iy(n)); hp = ydel(iy(p))
      sf = 3
    else
      hn = zdel(iz(n)); hp = zdel(iz(p))
      sf = 5
    end if

    do g = 1, ng
      ! Calculate transverse leakage moments
      Ln1(g) = Lp1(g)
      call TLUpd1 (u,p,g,Lp1(g))

      !Setup 2Gx2G matrix A and 2G vector B to obtain a2(g)
      Pn(g)  = 2._dp*D(n,g) / hn; Pp(g) = 2._dp*D(p,g) / hp
      do h = 1, ng
        if (h == g) then
          A(g,g)    = -Pn(g) * (Bcn(g,h)*Fn(g) + 1._dp)
          A(g,g+ng) =  Pp(g) * (Bcp(g,h)*Fp(g) + 1._dp)
        else
          A(g,h)    = -Pn(g)*Bcn(g,h)*Fn(g)
          A(g,h+ng) =  Pp(g)*Bcp(g,h)*Fp(g)
        end if
      end do
      b(g) = Pn(g) * (3._dp*a2n(g) + Gn(g)*a4n(g) + Fn(g)*Ln1(g)) &
           + Pp(g) * (3._dp*a2p(g) + Gp(g)*a4p(g) - Fp(g)*Lp1(g))
      ! write(*,*) Fn(g)*Ln1(g)
    end do
    ! stop
    do g = 1, ng
      do h = 1, ng
        if (h == g) then
          A(g+ng,g)    = dc(n,g,sf)   * (Bcn(g,h)*An(g) + 1._dp)
          A(g+ng,g+ng) = dc(p,g,sf+1) * (Bcp(g,h)*Ap(g) + 1._dp)
        else
          A(g+ng,h)    = dc(n,g,sf)   * Bcn(g,h) * An(g)
          A(g+ng,h+ng) = dc(p,g,sf+1) * Bcp(g,h) * Ap(g)
        end if
      end do
      ! Create vector b
        b(g+ng) = dc(p,g,sf+1) * (a2p(g) + a4p(g) + f0(p,g) - An(g)*Ln1(g)) &
                - dc(n,g,sf)   * (a2n(g) + a4n(g) + f0(n,g) + Ap(g)*Lp1(g))
    end do


  end subroutine get_a1matvec

  !****************************************************************************!

  function get_a3(cp,a1,Lmn1) result (a3)

    !Purpose: To  get a3 expansion coefficients


    use sdata, only: ng
    implicit none

    integer, intent(in)                 :: cp       ! to indicate left or right node
    real(dp), dimension(:), intent(in)  :: a1       ! a1 expansion coefficients
    real(dp), dimension(:), intent(in)  :: Lmn1     ! First transverse leakage moments
    real(dp), dimension(ng)             :: a3       ! a3 expansion coefficients

    real(dp) :: Bf
    real(dp), dimension(ng,ng)          :: Bc
    real(dp), dimension(ng)             :: Ac
    integer :: g, h

    if (cp == 1) then
      Bc = Bcn
      Ac = An
    else
      Bc = Bcp
      Ac = Ap
    end if

    do g = 1, ng
      Bf = 0._dp
      do h = 1, ng
        Bf = Bf + Bc(g,h)*a1(h)
      end do
      a3(g) = Ac(g) * (Bf + Lmn1(g))
    end do

  end function get_a3

  !****************************************************************************!

  function get_a4(a2) result (a4)

    !Purpose: To  get a4 expansion coefficients


    use sdata, only: ng


    implicit none

    real(dp), dimension(:), intent(in)  :: a2       ! a2 expansion coefficients
    real(dp), dimension(ng)             :: a4       ! a4 expansion coefficients

    real(dp) :: Bf
    integer :: g, h

    do g = 1, ng
      Bf = 0._dp
      do h = 1, ng
        Bf = Bf + Bcp(g,h) * a2(h)
      end do

      a4(g) = Bp(g) * (Bf + Lm2(g))
    end do

  end function get_a4

  !****************************************************************************!

  subroutine get_a2matvec(u,n,A,b)

    !Purpose: To setup GxG matrix and b vector to get a2 expansion coefficients


    use sdata, only: ng, D, f0, exsrc, xdel, ydel, zdel, ix, iy, iz

    implicit none

    integer, intent(in)                  :: u, n        ! Direction and node number
    real(dp), dimension(:,:), intent(out):: A           ! GxG Materix
    real(dp), dimension(:), intent(out)  :: b           ! G vector

    real(dp)                   :: S
    real(dp), dimension(ng)    :: Bf
    integer                    :: g, h

    !Setup GxG matrix and G vector to obtain a2(g)
    Bf = 0._dp
    do g = 1, ng
      !update zeroth source
      if (cmode == 2) then
        if (u == 1) then
          S = 0.25_dp * xdel(ix(n))**2 / D(n,g) * S1(n,g)
        else if (u == 2) then
          S = 0.25_dp * ydel(iy(n))**2 / D(n,g) * S2(n,g)
        else
          S = 0.25_dp * zdel(iz(n))**2 / D(n,g) * S3(n,g)
        end if
      else
        if (u == 1) then
          S = 0.25_dp * xdel(ix(n))**2 / D(n,g) * (S1(n,g) - exsrc(n,g))
        else if (u == 2) then
          S = 0.25_dp * ydel(iy(n))**2 / D(n,g) * (S2(n,g) - exsrc(n,g))
        else
          S = 0.25_dp * zdel(iz(n))**2 / D(n,g) * (S3(n,g) - exsrc(n,g))
        end if
      end if


      ! Create matrix A
      do h = 1, ng
        if (h == g) then
          A(g,g) = Bcp(g,h) * Ep(g) + 3._dp
        else
          A(g,h) = Bcp(g,h) * Ep(g)
        end if
        Bf(g) = Bf(g) + Bcp(g,h) * f0(n,h)
      end do

      ! Get second moment transverse leakage
      call TLUpd2 (u,n,g,Lm2(g))

      b(g) = Bf(g) - Ep(g)*Lm2(g) + S
    end do

  end subroutine get_a2matvec

  !******************************************************************************!

  function LU_solve(nt,msize,mat,b) result(x)

      !
      ! Purpose:
      !    To solve Ax=b by LU decomposition
      !

      USE io,      ONLY: ounit
      USE sdata,   ONLY: ix, iy, iz

      implicit none

      integer, intent(in)                  :: nt, msize  ! node and and matrix size
      real(dp), dimension(:,:), intent(in) :: mat           ! the matrix A
      real(dp), dimension(:), intent(in)   :: b            ! the vector b
      real(dp), dimension(msize)           :: x            ! the vector b

      real(dp), dimension(msize,msize) :: L, U
      real(dp), dimension(msize) :: y
      real(dp) :: piv, isum
      integer :: i, j, k

      U = mat
      L = 0._dp

      ! Start matrix decomposition
      do i= 1, msize
          if (ABS(mat(i,i)) < 10e-5) then
            write(ounit,*) 'ERROR IN MATRIX DECOMP: DIAGONAL ELEMENTS CLOSE TO ZERO'
            write(ounit,2001) ix(nt), iy(nt), iz(nt)
            write(*,*) 'ERROR IN MATRIX DECOMP: DIAGONAL ELEMENTS CLOSE TO ZERO'
            write(*,2001) ix(nt), iy(nt), iz(nt)
            STOP
          end if
          L(i,i) = 1._DP
          do j= i+1, msize
              piv = U(j,i)/U(i,i)
              L(j,i) = piv
              do k= i, msize
                  U(j,k) = U(j,k) - piv*U(i,k)
              end do
              U(j,i) = 0._dp
          end do
      end do


      !Solve y in Ly = b (Forward substitution)
      y(1) = b(1)
      do i=2,msize
          isum = 0._dp
          do k =1, i-1
              isum = isum + L(i,k)*y(k)
          end do
          y(i) = b(i)-isum
      end do

      ! Solve x in Ux=y(Backward substitution)
      x(msize) = y(msize)/U(msize,msize)
      do i = msize-1,1,-1
          isum = 0._dp
          do k =i+1,msize
              isum = isum + U(i,k)*x(k)
          end do
          x(i) = (y(i)-isum) / U(i,i)
      end do

      2001 FORMAT(2X, 'I = ', I2, ', J = ', I2, ', K = ', I2)

  end function LU_solve

  !******************************************************************************!

  SUBROUTINE Lxyz (n,g,L1, L2, L3)

  USE sdata, ONLY: nod, f0, xyz, ix, iy, iz, ystag, xstag, nzz, &
                   xeast, xwest, ysouth, ynorth, zbott, ztop, &
                   xdel, ydel, zdel

  ! Purpose:
     ! To update Transverse leakages for group g and nod n

  implicit none

  integer, intent(in)  :: n, g
  real(dp), intent(out) :: L1, L2, L3

  real(dp) :: jp, jm
  integer  :: p, m
  integer  :: i, j, k

  ! set i, j, k
  i = ix(n); j = iy(n); k = iz(n)

  ! x-direction zeroth transverse leakage
  if (i /= ystag(j)%smax) p = xyz(i+1,j,k)
  if (i /= ystag(j)%smin) m = xyz(i-1,j,k)

  if (i == ystag(j)%smax) then
    if (xeast == 2) then
      jp = 0._dp
    else
      jp = nod(n,g)%df(1)* f0(n,g) - nod(n,g)%dn(1)* f0(n,g)
    end if
  else
    jp = -nod(n,g)%df(1)*(f0(p,g) - f0(n,g)) - &
          nod(n,g)%dn(1)*(f0(p,g) + f0(n,g))
  end if
  if (i == ystag(j)%smin) then
    if (xwest == 2) then
      jm = 0._dp
    else
      jm = -nod(n,g)%df(2)*f0(n,g) - nod(n,g)%dn(2)* f0(n,g)
    end if
  else
    jm = -nod(n,g)%df(2)*(f0(n,g) - f0(m,g)) - &
          nod(n,g)%dn(2)*(f0(n,g) + f0(m,g))
  end if

  L1 = (jp - jm) / xdel(i)

  ! y-direction zeroth transverse leakage
  if (j /= xstag(i)%smax) p = xyz(i,j+1,k)
  if (j /= xstag(i)%smin) m = xyz(i,j-1,k)

  if (j == xstag(i)%smax) then
    if (ynorth == 2) then
      jp = 0._dp
    else
      jp = nod(n,g)%df(3)*f0(n,g) - nod(n,g)%dn(3)* f0(n,g)
    end if
  else
    jp = -nod(n,g)%df(3)*(f0(p,g) - f0(n,g)) - &
          nod(n,g)%dn(3)*(f0(p,g) + f0(n,g))
  end if
  if (j == xstag(i)%smin) then
    if (ysouth == 2) then
      jm = 0._dp
    else
      jm = -nod(n,g)%df(4)*f0(n,g) - nod(n,g)%dn(4)* f0(n,g)
    end if
  else
    jm = -nod(n,g)%df(4)*(f0(n,g) - f0(m,g)) - &
          nod(n,g)%dn(4)*(f0(n,g) + f0(m,g))
  end if


  L2 = (jp - jm) / ydel(j)

  ! z-direction zeroth transverse leakage
  if (k /= nzz)   p = xyz(i,j,k+1)
  if (k /= 1  )   m = xyz(i,j,k-1)

  if (k == nzz) then
    if (ztop == 2) then
      jp = 0._dp
    else
      jp = nod(n,g)%df(5)*f0(n,g) - nod(n,g)%dn(5)* f0(n,g)
    end if
  else
    jp = -nod(n,g)%df(5)*(f0(p,g) - f0(n,g)) - &
          nod(n,g)%dn(5)*(f0(p,g) + f0(n,g))
  end if
  if (k == 1) then
    if (zbott == 2) then
      jm = 0._dp
    else
      jm = -nod(n,g)%df(6)*f0(n,g) - nod(n,g)%dn(6)* f0(n,g)
    end if
  else
    jm = -nod(n,g)%df(6)*(f0(n,g) - f0(m,g)) - &
          nod(n,g)%dn(6)*(f0(n,g) + f0(m,g))
  end if

  L3 = (jp - jm) / zdel(k)


  end subroutine Lxyz

  !******************************************************************************!

  SUBROUTINE get_source ()

  USE sdata, ONLY: nnod, ng, exsrc

  ! Purpose:
     ! To update get source for the nodal update

  implicit none

  real(dp) :: L1, L2, L3
  integer  :: g, n
  logical  :: first = .true.

  if (first) then
    allocate(S1(nnod,ng), S2(nnod,ng), S3(nnod,ng))
    first = .false.
  end if

  do g = 1, ng
    do n = 1, nnod
      call Lxyz(n,g,L1,L2,L3)

      if (cmode == 2) then
        S1(n,g) =  L2 + L3 - exsrc(n,g)
        S2(n,g) =  L1 + L3 - exsrc(n,g)
        S3(n,g) =  L1 + L2 - exsrc(n,g)
      else
        S1(n,g) =  L2 + L3
        S2(n,g) =  L1 + L3
        S3(n,g) =  L1 + L2
      end if
    end do
  end do

  end subroutine get_source

  !******************************************************************************!

  SUBROUTINE TLUpd1 (u,n,g,Lmom1)

  USE sdata, ONLY: xdel, ydel, zdel, xstag, ystag, nzz, &
                   xwest, xeast, ynorth, ysouth, zbott, ztop, &
                   ix, iy, iz, xyz, D

  ! Purpose:
     ! To calaculate transverse leakage first moments

  implicit none

  integer, intent(in) :: u, n, g
  real(dp), intent(out) :: Lmom1

  real(dp) :: tm, tp
  real(dp) :: p1m, p2m, p1p, p2p, hp
  integer  :: p, m, i, j, k

  ! Set i, j, k
  i = ix(n); j = iy(n); k = iz(n)

  if (u==1) then
    ! Set paramaters for X-Direction Transverse leakage
    if (i /= ystag(j)%smax) p = xyz(i+1,j,k)
    if (i /= ystag(j)%smin) m = xyz(i-1,j,k)

    if (i == ystag(j)%smin) then
      if (xwest == 2) then
        tm = 1._dp
        tp = xdel(i+1)/xdel(i)
        p1m = tm + 1._dp; p2m = 2._dp*tm + 1._dp
        p1p = tp + 1._dp
        hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
        Lmom1 = (p1m * p2m * (S1(p,g) - S1(n,g))) / hp
      else
        tp = xdel(i+1)/xdel(i)
        p1p = tp + 1._dp
        Lmom1 = (S1(p,g) - S1(n,g)) / p1p
      end if
    else if (i == ystag(j)%smax) then
      if (xeast == 2) then
        tm = xdel(i-1)/xdel(i)
        tp = 1._dp
        p1m = tm + 1._dp
        p1p = tp + 1._dp; p2p = 2._dp*tp + 1._dp
        hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
        Lmom1 = (p1p * p2p * (S1(n,g) - S1(m,g))) / hp
      else
        tm = xdel(i-1)/xdel(i)
        p1m = tm + 1._dp
        Lmom1 = (S1(n,g) - S1(m,g)) / p1m
      end if
    else
      tm = xdel(i-1)/xdel(i)
      tp = xdel(i+1)/xdel(i)
      p1m = tm + 1._dp; p2m = 2._dp*tm + 1._dp
      p1p = tp + 1._dp; p2p = 2._dp*tp + 1._dp
      hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
      Lmom1 = (p1m * p2m * (S1(p,g) - S1(n,g)) &
            +  p1p * p2p * (S1(n,g) - S1(m,g))) / hp
    end if

    Lmom1 = 0.25_dp * xdel(i)**2 / D(n,g) * Lmom1

  else if (u == 2) then
    ! Set paramaters for Y-Direction Transverse leakage
    if (j /= xstag(i)%smax) p = xyz(i,j+1,k)
    if (j /= xstag(i)%smin) m = xyz(i,j-1,k)

    if (j == xstag(i)%smin) then
      if (ysouth == 2) then
        tm = 1._dp
        tp = ydel(j+1)/ydel(j)
        p1m = tm + 1._dp; p2m = 2._dp*tm + 1._dp
        p1p = tp + 1._dp
        hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
        Lmom1 = (p1m * p2m * (S2(p,g) - S2(n,g))) / hp
      else
        tp = ydel(j+1)/ydel(j)
        p1p = tp + 1._dp
        Lmom1 = (S2(p,g) - S2(n,g)) / p1p
      end if
    else if (j == xstag(i)%smax) then
      if (ynorth == 2) then
        tm = ydel(j-1)/ydel(j)
        tp = 1._dp
        p1m = tm + 1._dp
        p1p = tp + 1._dp; p2p = 2._dp*tp + 1._dp
        hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
        Lmom1 = (p1p * p2p * (S2(n,g) - S2(m,g))) / hp
      else
        tm = ydel(j-1)/ydel(j)
        p1m = tm + 1._dp
        Lmom1 = (S2(n,g) - S2(m,g)) / p1m
      end if
    else
      tm = ydel(j-1)/ydel(j)
      tp = ydel(j+1)/ydel(j)
      p1m = tm + 1._dp; p2m = 2._dp*tm + 1._dp
      p1p = tp + 1._dp; p2p = 2._dp*tp + 1._dp
      hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
      Lmom1 = (p1m * p2m * (S2(p,g) - S2(n,g)) &
            +  p1p * p2p * (S2(n,g) - S2(m,g))) / hp
    end if

    Lmom1 = 0.25_dp * ydel(j)**2 / D(n,g) * Lmom1

  else
    ! Set paramaters for Z-Direction Transverse leakage
    if (k /= nzz)   p = xyz(i,j,k+1)
    if (k /= 1  )   m = xyz(i,j,k-1)

    if (k == 1 ) then
      if (zbott == 2) then
        tm = 1._dp
        tp = zdel(k+1)/zdel(k)
        p1m = tm + 1._dp; p2m = 2._dp*tm + 1._dp
        p1p = tp + 1._dp
        hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
        Lmom1 = (p1m * p2m * (S3(p,g) - S3(n,g))) / hp
      else
        tp = zdel(k+1)/zdel(k)
        p1p = tp + 1._dp
        Lmom1 = (S3(p,g) - S3(n,g)) / p1p
      end if
    else if (k == nzz) then
      if (ztop == 2) then
        tm = zdel(k-1)/zdel(k)
        tp = 1._dp
        p1m = tm + 1._dp
        p1p = tp + 1._dp; p2p = 2._dp*tp + 1._dp
        hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
        Lmom1 = (p1p * p2p * (S3(n,g) - S3(m,g))) / hp
      else
        tm = zdel(k-1)/zdel(k)
        p1m = tm + 1._dp
        Lmom1 = (S3(n,g) - S3(m,g)) / p1m
      end if
    else
      tm = zdel(k-1)/zdel(k)
      tp = zdel(k+1)/zdel(k)
      p1m = tm + 1._dp; p2m = 2._dp*tm + 1._dp
      p1p = tp + 1._dp; p2p = 2._dp*tp + 1._dp
      hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
      Lmom1 = (p1m * p2m * (S3(p,g) - S3(n,g)) &
            +  p1p * p2p * (S3(n,g) - S3(m,g))) / hp
    end if

    Lmom1 = 0.25_dp * zdel(k)**2 / D(n,g) * Lmom1

  end if

  end subroutine TLUpd1

  !****************************************************************************!

  SUBROUTINE TLUpd2 (u,n,g,Lmom2)

  USE sdata, ONLY: xdel, ydel, zdel, xstag, ystag, nzz, &
                   xwest, xeast, ynorth, ysouth, zbott, ztop, &
                   ix, iy, iz, xyz, D

  ! Purpose:
     ! To calaculate transverse leakage second moments

  implicit none

  integer, intent(in) :: u, n, g
  real(dp), intent(out) :: Lmom2

  real(dp) :: tm, tp
  real(dp) :: p1m, p1p, hp
  integer  :: p, m, i, j, k

  ! Set i, j, k
  i = ix(n); j = iy(n); k = iz(n)

  if (u==1) then
    ! Set paramaters for X-Direction Transverse leakage
    if (i /= ystag(j)%smax) p = xyz(i+1,j,k)
    if (i /= ystag(j)%smin) m = xyz(i-1,j,k)

    if (i == ystag(j)%smin) then
      if (xwest == 2) then
        tm = 1._dp
        tp = xdel(i+1)/xdel(i)
        p1m = tm + 1._dp
        p1p = tp + 1._dp
        hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
        Lmom2 = (p1m * (S1(p,g) - S1(n,g))) / hp
      else
        Lmom2 = 0._dp
      end if
    else if (i == ystag(j)%smax) then
      if (xeast == 2) then
        tm = xdel(i-1)/xdel(i)
        tp = 1._dp
        p1m = tm + 1._dp
        p1p = tp + 1._dp
        hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
        Lmom2 = (p1p * (S1(m,g) - S1(n,g))) / hp
      else
        Lmom2 = 0._dp
      end if
    else
      tm = xdel(i-1)/xdel(i)
      tp = xdel(i+1)/xdel(i)
      p1m = tm + 1._dp
      p1p = tp + 1._dp
      hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
      Lmom2 = (p1m * (S1(p,g) - S1(n,g)) +  p1p * (S1(m,g) - S1(n,g))) / hp
    end if

    Lmom2 = 0.25_dp * xdel(i)**2 / D(n,g) * Lmom2

  else if (u == 2) then
    ! Set paramaters for Y-Direction Transverse leakage
    if (j /= xstag(i)%smax) p = xyz(i,j+1,k)
    if (j /= xstag(i)%smin) m = xyz(i,j-1,k)

    if (j == xstag(i)%smin) then
      if (ysouth == 2) then
        tm = 1._dp
        tp = ydel(j+1)/ydel(j)
        p1m = tm + 1._dp
        p1p = tp + 1._dp
        hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
        Lmom2 = (p1m * (S2(p,g) - S2(n,g))) / hp
      else
        Lmom2 = 0._dp
      end if
    else if (j == xstag(i)%smax) then
      if (ynorth == 2) then
        tm = ydel(j-1)/ydel(j)
        tp = 1._dp
        p1m = tm + 1._dp
        p1p = tp + 1._dp
        hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
        Lmom2 = (p1p * (S2(m,g) - S2(n,g))) / hp
      else
        Lmom2 = 0._dp
      end if
    else
      tm = ydel(j-1)/ydel(j)
      tp = ydel(j+1)/ydel(j)
      p1m = tm + 1._dp
      p1p = tp + 1._dp
      hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
      Lmom2 = (p1m * (S2(p,g) - S2(n,g)) +  p1p * (S2(m,g) - S2(n,g))) / hp
    end if

    Lmom2 = 0.25_dp * ydel(j)**2 / D(n,g) * Lmom2

  else
    ! Set paramaters for Z-Direction Transverse leakage
    if (k /= nzz) p = xyz(i,j,k+1)
    if (k /= 1  ) m = xyz(i,j,k-1)

    if (k == 1 ) then
      if (zbott == 2) then
        tm = 1._dp
        tp = zdel(k+1)/zdel(k)
        p1m = tm + 1._dp
        p1p = tp + 1._dp
        hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
        Lmom2 = (p1m * (S3(p,g) - S3(n,g))) / hp
      else
        Lmom2 = 0._dp
      end if
    else if (k == nzz) then
      if (ztop == 2) then
        tm = zdel(k-1)/zdel(k)
        tp = 1._dp
        p1m = tm + 1._dp
        p1p = tp + 1._dp
        hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
        Lmom2 = (p1p * (S3(m,g) - S3(n,g))) / hp
      else
        Lmom2 = 0._dp
      end if
    else
      tm = zdel(k-1)/zdel(k)
      tp = zdel(k+1)/zdel(k)
      p1m = tm + 1._dp
      p1p = tp + 1._dp
      hp = 2._dp * p1m * p1p* (tm + tp + 1._dp)
      Lmom2 = (p1m * (S3(p,g) - S3(n,g)) +  p1p * (S3(m,g) - S3(n,g))) / hp
    end if

    Lmom2 = 0.25_dp * zdel(k)**2 / D(n,g) * Lmom2

  end if


  end subroutine TLUpd2

  !****************************************************************************!

  function get_B(u,n) result (B)

    !Purpose: To  Buckling for node n and direction u


    use sdata, only: ng, xdel, ydel, zdel, ix, iy, iz, Ke, &
    sigr, D, chi, mat, nuf, sigs, Ke, tbeta, dfis
    implicit none

    integer, intent(in)                 :: u,n      ! direction and node number
    real(dp), dimension(ng,ng)          :: B        ! Buckling B2

    real(dp) :: dum, dn
    integer :: g, h

    if (u == 1) then
      dn = xdel(ix(n))
    else if (u == 2) then
      dn = ydel(iy(n))
    else
      dn = zdel(iz(n))
    end if

    if (cmode == 1) then
      do g = 1, ng
        do h = 1, ng
          if (g == h) then
            dum = sigr(n,g) - chi(mat(n),g)*nuf(n,h) / Ke
          else
            dum = -sigs(n,h,g) - chi(mat(n),g)*nuf(n,h) / Ke
          end if
          B(g,h) = 0.25_dp * dn**2 / D(n,g) * dum
        end do
      end do
    else if (cmode == 2) then
      do g = 1, ng
        do h = 1, ng
          if (g == h) then
            dum = sigr(n,g) - (1._dp - tbeta(mat(n)) + dfis(n)) &
                * chi(mat(n),g) * chi(mat(n),g)*nuf(n,h)
          else
            dum = -sigs(n,h,g) - (1._dp - tbeta(mat(n)) + dfis(n)) &
                * chi(mat(n),g)*nuf(n,h)
          end if
          B(g,h) = 0.25_dp * dn**2 / D(n,g) * dum
        end do
      end do
    else
      do g = 1, ng
        do h = 1, ng
          if (g == h) then
            dum = sigr(n,g) - chi(mat(n),g)*nuf(n,h) / Ke
          else
            dum = -sigs(n,g,h) - chi(mat(n),h)*nuf(n,g) / Ke
          end if
          B(g,h) = 0.25_dp * dn**2 / D(n,g) * dum
        end do
      end do
    end if

  end function get_B

  !****************************************************************************!

  SUBROUTINE get_ABEFGH (n,u,Aa,Bb,Ee,Ff,Gg,Hh)

  USE sdata, ONLY: ng, xdel, ydel, zdel, ix, iy, iz, D, sigr

  ! Purpose:
     ! To calaculate A,B,E,F,G,H paramters used to calculate matrix elements for
     ! nodal update

  implicit none

  integer, intent(in) :: u, n
  real(dp), dimension(:), intent(out) :: Aa,Bb,Ee,Ff,Gg,Hh

  real(dp) :: dn, alp, alp2
  real(dp) :: m1s, m0c, m2c
  integer  :: g

  if (u == 1) then
    dn = xdel(ix(n))
  else if (u == 2) then
    dn = ydel(iy(n))
  else
    dn = zdel(iz(n))
  end if

  do g = 1, ng
    !Calculate alpha and othe paramters to calculate A,B,E,F,G,H
    alp  = 0.5_dp * sqrt(sigr(n,g) / D(n,g)) * dn
    alp2 = alp**2
    m0c  = sinh(alp) / alp
    m1s  = 3._dp * (cosh(alp)/alp - sinh(alp)/alp2)
    m2c  = 5._dp * (sinh(alp)/alp - 3._dp*cosh(alp)/alp2 &
         + 3._dp*sinh(alp)/alp**3)

    Aa(g) = (sinh(alp) - m1s) / (alp2 * m1s)
    Bb(g) = (cosh(alp) - m0c - m2c) / (alp2 * m2c)
    Ee(g) = (m0c/m2c - 3._dp/alp2)
    Ff(g) = (alp*cosh(alp) - m1s) / (alp2 * m1s)
    Gg(g) = (alp*sinh(alp) - 3._dp*m2c) / (cosh(alp) - m0c - m2c)
    Hh(g) = (alp*cosh(alp) - m1s) / (sinh(alp) - m1s)
  end do

end subroutine get_ABEFGH

end module nodal
