module nodal

   use sdata, only: dp
   implicit none
   save

   integer :: cmode
   real(dp), dimension(:, :), allocatable :: bcn, bcp          ! Buckling for left and right nodes
   real(dp), dimension(:), allocatable :: an, bn, en, fn, gn, hn
   real(dp), dimension(:), allocatable :: ap, bp, ep, fp, gp, hp
   real(dp), dimension(:), allocatable :: lm2                ! Second order source
   real(dp), dimension(:, :), allocatable :: sx, sy, sz          ! zeroth order source in x, y, and z directions

   contains

   !****************************************************************************!

   subroutine nodal_update(cal_mode)

      !Purpose: to calculate flux expansion coefficients using SANM

      use sdata, only: ng, xyz, ystag, xstag, nyy, nzz, nxx, &
      a1n, a2n, a3n, a4n, a1p, a2p, a3p, a4p, ndmax, &
      ln1, lp1, xeast, xwest, ysouth, ynorth, zbott, ztop

      implicit none

      integer, intent(in) :: cal_mode  !Calc. mode -> 0=adjoint, 1=forward

      integer :: n, p
      integer :: i, j, k
      logical :: first = .true.

      if (first) then
         allocate(a1n(ng), a2n(ng), a3n(ng), a4n(ng))
         allocate(a1p(ng), a2p(ng), a3p(ng), a4p(ng))
         first = .false.
      end if

      allocate(ln1(ng), lp1(ng))
      allocate(bcn(ng, ng), bcp(ng, ng))
      allocate(an(ng), bn(ng), en(ng), fn(ng), gn(ng), hn(ng))
      allocate(ap(ng), bp(ng), ep(ng), fp(ng), gp(ng), hp(ng))
      allocate(lm2(ng))

      cmode = cal_mode

      call get_source()


      !Node sweeps in x-direction
      do k = 1, nzz
         do j = 1, nyy
            p = xyz(ystag(j)%smin, j, k)
            bcp = get_b(1, p)
            call get_abefgh (p, 1, ap, bp, ep, fp, gp, hp)
            call get_coefs_first(xwest, 1, p)
            call nodal_coup_upd(1, a1p, a2p, a3p, a4p, p = p)
            do i = ystag(j)%smin, ystag(j)%smax - 1
               n = xyz(i, j, k) ;
               p = xyz(i + 1, j, k)
               bcn = bcp
               an = ap;
               bn = bp;
               en = ep;
               fn = fp;
               gn = gp;
               hn = hp
               bcp = get_b(1, p)
               call get_abefgh (p, 1, ap, bp, ep, fp, gp, hp)
               call get_coefs(1, n, p)
               call nodal_coup_upd(1, a1n, a2n, a3n, a4n, n, p)
            end do
            n = xyz(ystag(j)%smax, j, k)
            bcn = bcp
            an = ap;
            bn = bp;
            en = ep;
            fn = fp;
            gn = gp;
            hn = hp
            call get_coefs_last(xeast, 1, n)
            call nodal_coup_upd(1, a1n, a2n, a3n, a4n, n = n)
         end do
      end do

      !Node sweeps in y-direction
      do k = 1, nzz
         do i = 1, nxx
            p = xyz(i, xstag(i)%smin, k)
            bcp = get_b(2, p)
            call get_abefgh (p, 2, ap, bp, ep, fp, gp, hp)
            call get_coefs_first(ysouth, 2, p)
            call nodal_coup_upd(2, a1p, a2p, a3p, a4p, p = p)
            do j = xstag(i)%smin, xstag(i)%smax - 1
               n = xyz(i, j, k) ;
               p = xyz(i, j + 1, k)
               bcn = bcp
               an = ap;
               bn = bp;
               en = ep;
               fn = fp;
               gn = gp;
               hn = hp
               bcp = get_b(2, p)
               call get_abefgh (p, 2, ap, bp, ep, fp, gp, hp)
               call get_coefs(2, n, p)
               call nodal_coup_upd(2, a1n, a2n, a3n, a4n, n, p)
            end do
            n = xyz(i, xstag(i)%smax, k)
            bcn = bcp
            an = ap;
            bn = bp;
            en = ep;
            fn = fp;
            gn = gp;
            hn = hp
            call get_coefs_last(ynorth, 2, n)
            call nodal_coup_upd(2, a1n, a2n, a3n, a4n, n = n)
         end do
      end do

      !Node sweeps in z-direction
      do j = 1, nyy
         do i = ystag(j)%smin, ystag(j)%smax
            p = xyz(i, j, 1)
            bcp = get_b(3, p)
            call get_abefgh (p, 3, ap, bp, ep, fp, gp, hp)
            call get_coefs_first(zbott, 3, p)
            call nodal_coup_upd(3, a1p, a2p, a3p, a4p, p = p)
            do k = 1, nzz - 1
               n = xyz(i, j, k) ;
               p = xyz(i, j, k + 1)
               bcn = bcp
               an = ap;
               bn = bp;
               en = ep;
               fn = fp;
               gn = gp;
               hn = hp
               bcp = get_b(3, p)
               call get_abefgh (p, 3, ap, bp, ep, fp, gp, hp)
               call get_coefs(3, n, p)
               call nodal_coup_upd(3, a1n, a2n, a3n, a4n, n, p)
            end do
            n = xyz(i, j, nzz)
            bcn = bcp
            an = ap;
            bn = bp;
            en = ep;
            fn = fp;
            gn = gp;
            hn = hp
            call get_coefs_last(ztop, 3, n)
            call nodal_coup_upd(3, a1n, a2n, a3n, a4n, n = n)
         end do
      end do

      call check_stability(ndmax)

      deallocate(ln1, lp1, bcp, bcn, lm2)
      deallocate(an, bn, en, fn, gn, hn)
      deallocate(ap, bp, ep, fp, gp, hp)

   end subroutine nodal_update

   !****************************************************************************!

   subroutine nodal_update_pnm(cal_mode)

      !Purpose: to calculate flux expansion coefficients using PNM

      use sdata, only: ng, xyz, ystag, xstag, nyy, nzz, nxx, &
      a1n, a2n, a3n, a4n, a1p, a2p, a3p, a4p, ndmax, &
      ln1, lp1, xeast, xwest, ysouth, ynorth, zbott, ztop

      implicit none

      integer, intent(in) :: cal_mode  !Calc. mode -> 0=adjoint, 1=forward

      integer :: n, p
      integer :: i, j, k
      logical :: first = .true.

      if (first) then
         allocate(a1n(ng), a2n(ng), a3n(ng), a4n(ng))
         allocate(a1p(ng), a2p(ng), a3p(ng), a4p(ng))
         first = .false.
      end if

      allocate(ln1(ng), lp1(ng))
      allocate(bcn(ng, ng), bcp(ng, ng))
      allocate(an(ng), bn(ng), en(ng), fn(ng), gn(ng), hn(ng))
      allocate(ap(ng), bp(ng), ep(ng), fp(ng), gp(ng), hp(ng))
      allocate(lm2(ng))

      an = 1._dp / 15._dp;
      bn = 1._dp / 35._dp;
      en = 2._dp / 7._dp
      fn = 2._dp / 5._dp;
      gn = 10._dp;
      hn = 6._dp
      ap = an;
      bp = bn;
      ep = en;
      fp = fn;
      gp = gn;
      hp = hn

      cmode = cal_mode

      call get_source()


      !Node sweeps in x-direction
      do k = 1, nzz
         do j = 1, nyy
            p = xyz(ystag(j)%smin, j, k)
            bcp = get_b(1, p)
            call get_coefs_first(xwest, 1, p)
            call nodal_coup_upd(1, a1p, a2p, a3p, a4p, p = p)
            do i = ystag(j)%smin, ystag(j)%smax - 1
               n = xyz(i, j, k) ;
               p = xyz(i + 1, j, k)
               bcn = bcp
               an = ap;
               bn = bp;
               en = ep;
               fn = fp;
               gn = gp;
               hn = hp
               bcp = get_b(1, p)
               call get_coefs(1, n, p)
               call nodal_coup_upd(1, a1n, a2n, a3n, a4n, n, p)
            end do
            n = xyz(ystag(j)%smax, j, k)
            bcn = bcp
            an = ap;
            bn = bp;
            en = ep;
            fn = fp;
            gn = gp;
            hn = hp
            call get_coefs_last(xeast, 1, n)
            call nodal_coup_upd(1, a1n, a2n, a3n, a4n, n = n)
         end do
      end do

      !Node sweeps in y-direction
      do k = 1, nzz
         do i = 1, nxx
            p = xyz(i, xstag(i)%smin, k)
            bcp = get_b(2, p)
            call get_coefs_first(ysouth, 2, p)
            call nodal_coup_upd(2, a1p, a2p, a3p, a4p, p = p)
            do j = xstag(i)%smin, xstag(i)%smax - 1
               n = xyz(i, j, k) ;
               p = xyz(i, j + 1, k)
               bcn = bcp
               an = ap;
               bn = bp;
               en = ep;
               fn = fp;
               gn = gp;
               hn = hp
               bcp = get_b(2, p)
               call get_coefs(2, n, p)
               call nodal_coup_upd(2, a1n, a2n, a3n, a4n, n, p)
            end do
            n = xyz(i, xstag(i)%smax, k)
            bcn = bcp
            an = ap;
            bn = bp;
            en = ep;
            fn = fp;
            gn = gp;
            hn = hp
            call get_coefs_last(ynorth, 2, n)
            call nodal_coup_upd(2, a1n, a2n, a3n, a4n, n = n)
         end do
      end do

      !Node sweeps in z-direction
      do j = 1, nyy
         do i = ystag(j)%smin, ystag(j)%smax
            p = xyz(i, j, 1)
            bcp = get_b(3, p)
            call get_coefs_first(zbott, 3, p)
            call nodal_coup_upd(3, a1p, a2p, a3p, a4p, p = p)
            do k = 1, nzz - 1
               n = xyz(i, j, k) ;
               p = xyz(i, j, k + 1)
               bcn = bcp
               an = ap;
               bn = bp;
               en = ep;
               fn = fp;
               gn = gp;
               hn = hp
               bcp = get_b(3, p)
               call get_coefs(3, n, p)
               call nodal_coup_upd(3, a1n, a2n, a3n, a4n, n, p)
            end do
            n = xyz(i, j, nzz)
            bcn = bcp
            an = ap;
            bn = bp;
            en = ep;
            fn = fp;
            gn = gp;
            hn = hp
            call get_coefs_last(ztop, 3, n)
            call nodal_coup_upd(3, a1n, a2n, a3n, a4n, n = n)
         end do
      end do

      call check_stability(ndmax)

      deallocate(ln1, lp1, bcp, bcn, lm2)
      deallocate(an, bn, en, fn, gn, hn)
      deallocate(ap, bp, ep, fp, gp, hp)

   end subroutine nodal_update_pnm

   !****************************************************************************!

   subroutine check_stability(ndmax)

      !
      ! Purpose:
      !    To check the stability of the two-node non linear iteration. If it's not
      !    stable provide an error message and how to overcome the error

      implicit none

      ! Maximum change in nodal coupling coefficients, if large means not stable
      real(dp), intent(in) :: ndmax

      if (ndmax > 1.e3) then
         write(*, *)
         write(*, 1236) ndmax
         write(*, *) "The two-node nonlinear iteration seems not stable. Try these:"
         write(*, *) "1. Change iteration control using %ITER card, "
         write(*, *) "   Perhaps by making nodal update less frequent,"
         write(*, *) "   or increase number inner iteration per outer iteration."
         write(*, *) "2. Ensure the node sizes in all directions are as uniform as possible. "
         write(*, *) "   Also try smaller node size."
         write(*, *) "3. For transient problem, try to reduce time step size."
         write(*, *)
         write(*, *) "If this error persists, contact me at makrus.imron@gmail.com"
         write(*, *) "Thank you for using KOMODO"
         stop
      end if

      1236 format(" Error: Max. change in nodal coupling coefficient = ", f10.1)

   end subroutine check_stability

   !******************************************************************************!

   subroutine nodal_coup_upd(u, a1, a2, a3, a4, n, p)

      !Purpose: to update nodal coupling coefficients

      use sdata, only: ng, ix, iy, iz, d, xdel, ydel, zdel, nod, f0, &
      ndmax, im, jm, km

      implicit none

      integer, intent(in) :: u
      integer, intent(in), optional :: n, p
      real(dp), dimension(:), intent(in) :: a1, a2, a3, a4

      integer :: g, sf
      real(dp) :: dh, jp, nder, ndpr

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
            jp = -2._dp * d(n, g) / dh * (a1(g) + 3._dp * a2(g) + hn(g) * a3(g) + gn(g) * a4(g))

            ! Update nodal coupling
            ndpr = nod(n, g)%dn(sf)
            nod(n, g)%dn(sf) = (nod(n, g)%df(sf) * (f0(n, g) - f0(p, g)) - jp) &
             / (f0(n, g) + f0(p, g))
            nod(p, g)%dn(sf + 1) = nod(n, g)%dn(sf)

            ! Check max difference on new nodal coupling coefficients
            nder = abs(nod(n, g)%dn(sf) - ndpr)
            if (nder > ndmax) then
               ndmax = nder
               im = ix(n) ;
               jm = iy(n) ;
               km = iz(n)
            end if
         end do
      else if (present(p)) then
         do g = 1, ng
            jp = -2._dp * d(p, g) / dh * (a1(g) - 3._dp * a2(g) + hp(g) * a3(g) - gp(g) * a4(g))

            ! Update nodal coupling
            ndpr = nod(p, g)%dn(sf + 1)
            nod(p, g)%dn(sf + 1) = -(jp / f0(p, g) + nod(p, g)%df(sf + 1))

            ! Check max difference on new nodal coupling coefficients
            nder = abs(nod(p, g)%dn(sf + 1) - ndpr)
            if (nder > ndmax) then
               ndmax = nder
               im = ix(p) ;
               jm = iy(p) ;
               km = iz(p)
            end if
         end do
      else
         do g = 1, ng
            jp = -2._dp * d(n, g) / dh * (a1(g) + 3._dp * a2(g) + hn(g) * a3(g) + gn(g) * a4(g))

            ! Update nodal coupling
            ndpr = nod(n, g)%dn(sf)
            nod(n, g)%dn(sf) = -(jp / f0(n, g) - nod(n, g)%df(sf))

            ! Check max difference on new nodal coupling coefficients
            nder = abs(nod(n, g)%dn(sf) - ndpr)
            if (nder > ndmax) then
               ndmax = nder
               im = ix(n) ;
               jm = iy(n) ;
               km = iz(n)
            end if
         end do
      end if

   end subroutine nodal_coup_upd

   !****************************************************************************!

   subroutine get_coefs_first(bc, u, p)

      !Purpose: to calculate flux expansion coefficients

      use sdata, only: ng, a1p, a2p, a3p, a4p, lp1

      implicit none

      integer, intent(in) :: bc
      integer, intent(in) :: u, p

      real(dp), dimension(ng, ng) :: a           ! GxG Matrix
      real(dp), dimension(ng) :: b           ! G vector

      !Setup GxG matrix and G vector to obtain a2(g) for node p
      call get_a2matvec(u, p, a, b)
      !calculate a2 expansion coefficients4
      a2p = lu_solve(p, ng, a, b)
      !calculate a4 expansion coefficients
      a4p = get_a4(a2p)

      !Setup GxG matrix and G vector to obtain a1(g) for left right node
      call get_a1matvec_first(bc, u, p, a2p, a4p, a, b)

      !calculate a1 expansion coefficients
      a1p = lu_solve(p, ng, a, b)

      !calculate a3 expansion coefficients
      a3p = get_a3(2, a1p, lp1)

   end subroutine get_coefs_first

   !****************************************************************************!

   subroutine get_coefs_last(bc, u, n)

      !Purpose: to calculate flux expansion coefficients

      use sdata, only: ng, a1n, a2n, a3n, a4n, a2p, a4p, ln1

      implicit none

      integer, intent(in) :: bc
      integer, intent(in) :: u, n

      real(dp), dimension(ng, ng) :: a           ! GxG Matrix
      real(dp), dimension(ng) :: b           ! G vector

      !get a2n and a4n expansion coefficients
      a2n = a2p
      a4n = a4p

      !Setup GxG matrix and G vector to obtain a1(g) for most right node
      call get_a1matvec_last(bc, u, n, a2n, a4n, a, b)

      !calculate a1 expansion coefficients
      a1n = lu_solve(n, ng, a, b)

      !calculate a3 expansion coefficients
      a3n = get_a3(1, a1n, ln1)

   end subroutine get_coefs_last

   !****************************************************************************!

   subroutine get_coefs(u, n, p)

      !Purpose: to calculate flux expansion coefficients

      use sdata, only: ng, a1n, a2n, a3n, a4n, a2p, a4p, ln1

      implicit none

      integer, intent(in) :: u, n, p

      real(dp), dimension(ng, ng) :: a           ! GxG Matrix
      real(dp), dimension(ng) :: b           ! G vector
      real(dp), dimension(2 * ng, 2 * ng) :: r           ! 2Gx2G Matrix
      real(dp), dimension(2 * ng) :: s           ! 2G vector

      !get a2n and a4n expansion coefficients
      a2n = a2p
      a4n = a4p

      !Setup GxG matrix and G vector to obtain a2(g) for node p
      call get_a2matvec(u, p, a, b)
      !calculate a2 expansion coefficients
      a2p = lu_solve(p, ng, a, b)
      !calculate a4 expansion coefficients
      a4p = get_a4(a2p)

      !Setup 2Gx2G matrix and 2G vector to obtain a1(g) for node n
      call get_a1matvec(u, n, p, a2n, a4n, a2p, a4p, r, s)

      !calculate a1 expansion coefficients
      s = lu_solve(n, 2 * ng, r, s)
      a1n = s(1:ng)

      !calculate a3 expansion coefficients
      a3n = get_a3(1, a1n, ln1)

   end subroutine get_coefs

   !****************************************************************************!

   subroutine get_a1matvec_first(bc, u, p, a2p, a4p, a, b)

      !Purpose: To get matrix vector to calculate a1 for most left node

      use sdata, only: ng, xdel, ydel, zdel, ix, iy, iz, f0, d, lp1, dc

      implicit none

      integer, intent(in) :: bc
      integer, intent(in) :: u, p         ! Direction and node n umber
      real(dp), dimension(:), intent(in) :: a2p, a4p     ! a2 and a4 expansion coefficients
      real(dp), dimension(:, :), intent(out) :: a            ! 2Gx2G Materix
      real(dp), dimension(:), intent(out) :: b            ! 2G vector

      real(dp) :: dn
      real(dp) :: pp                       !dif. coef/dx
      integer :: g, h, sf

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
         call tlupd1 (u, p, g, lp1(g))
         pp = 2._dp * d(p, g) / dn

         !Setup GxG matrix A and G vector B to obtain a2(g)
         if (bc == 2) then
            do h = 1, ng
               if (h == g) then
                  a(g, g) = pp * (bcp(g, h) * fp(g) + 1._dp)
               else
                  a(g, h) = pp * bcp(g, h) * fp(g)
               end if
            end do
            b(g) = pp * (3._dp * a2p(g) + gp(g) * a4p(g) - fp(g) * lp1(g))
         else if (bc == 1) then
            do h = 1, ng
               if (h == g) then
                  a(g, g) = -dc(p, g, sf) * (1._dp + ap(g) * bcp(g, h)) &
                  - 2._dp * pp * (ap(g) * bcp(g, h) * hp(g) + 1._dp)
               else
                  a(g, h) = -dc(p, g, sf) * ap(g) * bcp(g, h) - 2._dp * pp * ap(g) * bcp(g, h) * hp(g)
               end if
            end do
            b(g) = 2._dp * pp * (ap(g) * hp(g) * lp1(g) - 3._dp * a2p(g) - gp(g) * a4p(g)) &
            - dc(p, g, sf) * (a2p(g) + a4p(g) + f0(p, g) - ap(g) * lp1(g))
         else
            do h = 1, ng
               if (h == g) then
                  a(g, g) = dc(p, g, sf) * (1._dp + ap(g) * bcp(g, h))
               else
                  a(g, h) = dc(p, g, sf) * ap(g) * bcp(g, h)
               end if
            end do
            b(g) = dc(p, g, sf) * (a2p(g) + a4p(g) + f0(p, g) - ap(g) * lp1(g))
         end if
      end do

   end subroutine get_a1matvec_first

   !****************************************************************************!

   subroutine get_a1matvec_last(bc, u, n, a2n, a4n, a, b)

      !Purpose: To get matrix vector to calculate a1 for most right node

      use sdata, only: ng, xdel, ydel, zdel, ix, iy, iz, f0, d, lp1, ln1, dc

      implicit none

      integer, intent(in) :: bc
      integer, intent(in) :: u, n         ! Direction and node number
      real(dp), dimension(:), intent(in) :: a2n, a4n           ! a2 and a4 expansion coefficients
      real(dp), dimension(:, :), intent(out) :: a           ! 2Gx2G Materix
      real(dp), dimension(:), intent(out) :: b           ! 2G vector

      real(dp) :: dn
      real(dp) :: pn                       !dif. coef/dx
      integer :: g, h, sf

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
         ln1(g) = lp1(g)
         pn = 2._dp * d(n, g) / dn

         !Setup 2Gx2G matrix A and 2G vector B to obtain a2(g)
         if (bc == 2) then
            do h = 1, ng
               if (h == g) then
                  a(g, g) = -pn * (bcn(g, h) * fn(g) + 1._dp)
               else
                  a(g, h) = -pn * bcn(g, h) * fn(g)
               end if
            end do
            b(g) = pn * (3._dp * a2n(g) + gn(g) * a4n(g) + fn(g) * ln1(g))
         else if (bc == 1) then
            do h = 1, ng
               if (h == g) then
                  a(g, g) = dc(n, g, sf) * (1._dp + an(g) * bcn(g, h)) &
                  + 2._dp * pn * (an(g) * bcn(g, h) * hn(g) + 1._dp)
               else
                  a(g, h) = dc(n, g, sf) * an(g) * bcn(g, h) + 2._dp * pn * an(g) * bcn(g, h) * hn(g)
               end if
            end do
            b(g) = -2._dp * pn * (an(g) * hn(g) * ln1(g) + 3._dp * a2n(g) + gn(g) * a4n(g)) &
            - dc(n, g, sf) * (a2n(g) + a4n(g) + f0(n, g) + an(g) * ln1(g))
         else
            do h = 1, ng
               if (h == g) then
                  a(g, g) = dc(n, g, sf) * (1._dp + an(g) * bcn(g, h))
               else
                  a(g, h) = dc(n, g, sf) * an(g) * bcn(g, h)
               end if
            end do
            b(g) = -dc(n, g, sf) * (a2n(g) + a4n(g) + f0(n, g) + an(g) * ln1(g))
         end if
      end do


   end subroutine get_a1matvec_last

   !****************************************************************************!

   subroutine get_a1matvec(u, n, p, a2n, a4n, a2p, a4p, a, b)

      !Purpose: To setup 2Gx2G matrix and 2G vector to get a1 expansion coefficients

      use sdata, only: ng, xdel, ydel, zdel, ix, iy, iz, d, f0, ln1, lp1, dc

      implicit none

      integer, intent(in) :: u, n, p         ! Direction and node number
      real(dp), dimension(:), intent(in) :: a2n, a4n           ! a2 and a4 expansion coefficients
      real(dp), dimension(:), intent(in) :: a2p, a4p           ! a2 and a4 expansion coefficients
      real(dp), dimension(:, :), intent(out) :: a           ! 2Gx2G Materix
      real(dp), dimension(:), intent(out) :: b           ! 2G vector

      real(dp) :: hn, hp
      real(dp), dimension(ng) :: pn, pp              !dif. coef/dx
      integer :: g, h, sf

      ! define node size in direction u
      if (u == 1) then
         hn = xdel(ix(n)) ;
         hp = xdel(ix(p))
         sf = 1
      else if (u == 2) then
         hn = ydel(iy(n)) ;
         hp = ydel(iy(p))
         sf = 3
      else
         hn = zdel(iz(n)) ;
         hp = zdel(iz(p))
         sf = 5
      end if

      do g = 1, ng
         ! Calculate transverse leakage moments
         ln1(g) = lp1(g)
         call tlupd1 (u, p, g, lp1(g))

         !Setup 2Gx2G matrix A and 2G vector B to obtain a2(g)
         pn(g) = 2._dp * d(n, g) / hn;
         pp(g) = 2._dp * d(p, g) / hp
         do h = 1, ng
            if (h == g) then
               a(g, g) = -pn(g) * (bcn(g, h) * fn(g) + 1._dp)
               a(g, g + ng) = pp(g) * (bcp(g, h) * fp(g) + 1._dp)
            else
               a(g, h) = -pn(g) * bcn(g, h) * fn(g)
               a(g, h + ng) = pp(g) * bcp(g, h) * fp(g)
            end if
         end do
         b(g) = pn(g) * (3._dp * a2n(g) + gn(g) * a4n(g) + fn(g) * ln1(g)) &
         + pp(g) * (3._dp * a2p(g) + gp(g) * a4p(g) - fp(g) * lp1(g))
      end do
      do g = 1, ng
         do h = 1, ng
            if (h == g) then
               a(g + ng, g) = dc(n, g, sf) * (bcn(g, h) * an(g) + 1._dp)
               a(g + ng, g + ng) = dc(p, g, sf + 1) * (bcp(g, h) * ap(g) + 1._dp)
            else
               a(g + ng, h) = dc(n, g, sf) * bcn(g, h) * an(g)
               a(g + ng, h + ng) = dc(p, g, sf + 1) * bcp(g, h) * ap(g)
            end if
         end do
         ! Create vector b
         b(g + ng) = dc(p, g, sf + 1) * (a2p(g) + a4p(g) + f0(p, g) - an(g) * ln1(g)) &
         - dc(n, g, sf) * (a2n(g) + a4n(g) + f0(n, g) + ap(g) * lp1(g))
      end do


   end subroutine get_a1matvec

   !****************************************************************************!

   function get_a3(cp, a1, lmn1) result (a3)

      !Purpose: To  get a3 expansion coefficients

      use sdata, only: ng
      implicit none

      integer, intent(in) :: cp       ! to indicate whether it is first node (=2 means first node)
      real(dp), dimension(:), intent(in) :: a1       ! a1 expansion coefficients
      real(dp), dimension(:), intent(in) :: lmn1     ! First transverse leakage moments
      real(dp), dimension(ng) :: a3       ! a3 expansion coefficients

      real(dp) :: bf
      real(dp), dimension(ng, ng) :: bc
      real(dp), dimension(ng) :: ac
      integer :: g, h

      if (cp == 1) then
         bc = bcn
         ac = an
      else
         bc = bcp    ! Take buckling values from the left node
         ac = ap
      end if

      do g = 1, ng
         bf = 0._dp
         do h = 1, ng
            bf = bf + bc(g, h) * a1(h)
         end do
         a3(g) = ac(g) * (bf + lmn1(g))
      end do

   end function get_a3

   !****************************************************************************!

   function get_a4(a2) result (a4)

      !Purpose: To  get a4 expansion coefficients

      use sdata, only: ng

      implicit none

      real(dp), dimension(:), intent(in) :: a2       ! a2 expansion coefficients
      real(dp), dimension(ng) :: a4       ! a4 expansion coefficients

      real(dp) :: bf
      integer :: g, h

      do g = 1, ng
         bf = 0._dp
         do h = 1, ng
            bf = bf + bcp(g, h) * a2(h)
         end do

         a4(g) = bp(g) * (bf + lm2(g))
      end do

   end function get_a4

   !****************************************************************************!

   subroutine get_a2matvec(u, n, a, b)

      !Purpose: To setup GxG matrix and b vector to get a2 expansion coefficients

      use sdata, only: ng, d, f0, exsrc, xdel, ydel, zdel, ix, iy, iz

      implicit none

      integer, intent(in) :: u, n        ! Direction and node number
      real(dp), dimension(:, :), intent(out) :: a           ! GxG Materix
      real(dp), dimension(:), intent(out) :: b           ! G vector

      real(dp) :: s
      real(dp), dimension(ng) :: bf
      integer :: g, h

      !Setup GxG matrix and G vector to obtain a2(g)
      bf = 0._dp
      do g = 1, ng
         !update zeroth source
         if (cmode == 2) then
            if (u == 1) then
               s = 0.25_dp * xdel(ix(n))**2 / d(n, g) * sx(n, g)
            else if (u == 2) then
               s = 0.25_dp * ydel(iy(n))**2 / d(n, g) * sy(n, g)
            else
               s = 0.25_dp * zdel(iz(n))**2 / d(n, g) * sz(n, g)
            end if
         else
            if (u == 1) then
               s = 0.25_dp * xdel(ix(n))**2 / d(n, g) * (sx(n, g) - exsrc(n, g))
            else if (u == 2) then
               s = 0.25_dp * ydel(iy(n))**2 / d(n, g) * (sy(n, g) - exsrc(n, g))
            else
               s = 0.25_dp * zdel(iz(n))**2 / d(n, g) * (sz(n, g) - exsrc(n, g))
            end if
         end if


         ! Create matrix A
         do h = 1, ng
            if (h == g) then
               a(g, g) = bcp(g, h) * ep(g) + 3._dp
            else
               a(g, h) = bcp(g, h) * ep(g)
            end if
            bf(g) = bf(g) + bcp(g, h) * f0(n, h)
         end do

         ! Get second moment transverse leakage
         call tlupd2 (u, n, g, lm2(g))

         b(g) = bf(g) - ep(g) * lm2(g) + s
      end do

   end subroutine get_a2matvec

   !******************************************************************************!

   function lu_solve(nt, msize, mat, b) result(x)

      !
      ! Purpose:
      !    To solve Ax=b by LU decomposition
      !

      use io, only: ounit
      use sdata, only: ix, iy, iz

      implicit none

      integer, intent(in) :: nt, msize  ! node and and matrix size
      real(dp), dimension(:, :), intent(in) :: mat           ! the matrix A
      real(dp), dimension(:), intent(in) :: b            ! the vector b
      real(dp), dimension(msize) :: x            ! the vector b

      real(dp), dimension(msize, msize) :: l, u
      real(dp), dimension(msize) :: y
      real(dp) :: piv, isum
      integer :: i, j, k

      u = mat
      l = 0._dp

      ! Start matrix decomposition
      do i = 1, msize
         if (abs(mat(i, i)) < 10e-5) then
            write(ounit, *) 'ERROR IN MATRIX DECOMP: DIAGONAL ELEMENTS CLOSE TO ZERO'
            write(ounit, 2001) ix(nt), iy(nt), iz(nt)
            write(*, *) 'ERROR IN MATRIX DECOMP: DIAGONAL ELEMENTS CLOSE TO ZERO'
            write(*, 2001) ix(nt), iy(nt), iz(nt)
            stop
         end if
         l(i, i) = 1._dp
         do j = i + 1, msize
            piv = u(j, i) / u(i, i)
            l(j, i) = piv
            do k = i, msize
               u(j, k) = u(j, k) - piv * u(i, k)
            end do
            u(j, i) = 0._dp
         end do
      end do


      !Solve y in Ly = b (Forward substitution)
      y(1) = b(1)
      do i = 2, msize
         isum = 0._dp
         do k = 1, i - 1
            isum = isum + l(i, k) * y(k)
         end do
         y(i) = b(i) - isum
      end do

      ! Solve x in Ux=y(Backward substitution)
      x(msize) = y(msize) / u(msize, msize)
      do i = msize - 1, 1, -1
         isum = 0._dp
         do k = i + 1, msize
            isum = isum + u(i, k) * x(k)
         end do
         x(i) = (y(i) - isum) / u(i, i)
      end do

      2001 format(2x, 'I = ', i2, ', J = ', i2, ', K = ', i2)

   end function lu_solve

   !******************************************************************************!

   subroutine lxyz (n, g, l1, l2, l3)

      use sdata, only: nod, f0, xyz, ix, iy, iz, ystag, xstag, nzz, &
      xeast, xwest, ysouth, ynorth, zbott, ztop, &
      xdel, ydel, zdel

      ! Purpose:
      ! To update Transverse leakages for group g and nod n

      implicit none

      integer, intent(in) :: n, g
      real(dp), intent(out) :: l1, l2, l3

      real(dp) :: jp, jm
      integer :: p, m
      integer :: i, j, k

      ! set i, j, k
      i = ix(n) ;
      j = iy(n) ;
      k = iz(n)

      ! x-direction zeroth transverse leakage
      if (i /= ystag(j)%smax) p = xyz(i + 1, j, k)
      if (i /= ystag(j)%smin) m = xyz(i - 1, j, k)

      if (i == ystag(j)%smax) then
         if (xeast == 2) then
            jp = 0._dp
         else
            jp = nod(n, g)%df(1) * f0(n, g) - nod(n, g)%dn(1) * f0(n, g)
         end if
      else
         jp = -nod(n, g)%df(1) * (f0(p, g) - f0(n, g)) - &
         nod(n, g)%dn(1) * (f0(p, g) + f0(n, g))
      end if
      if (i == ystag(j)%smin) then
         if (xwest == 2) then
            jm = 0._dp
         else
            jm = -nod(n, g)%df(2) * f0(n, g) - nod(n, g)%dn(2) * f0(n, g)
         end if
      else
         jm = -nod(n, g)%df(2) * (f0(n, g) - f0(m, g)) - &
         nod(n, g)%dn(2) * (f0(n, g) + f0(m, g))
      end if

      l1 = (jp - jm) / xdel(i)

      ! y-direction zeroth transverse leakage
      if (j /= xstag(i)%smax) p = xyz(i, j + 1, k)
      if (j /= xstag(i)%smin) m = xyz(i, j - 1, k)

      if (j == xstag(i)%smax) then
         if (ynorth == 2) then
            jp = 0._dp
         else
            jp = nod(n, g)%df(3) * f0(n, g) - nod(n, g)%dn(3) * f0(n, g)
         end if
      else
         jp = -nod(n, g)%df(3) * (f0(p, g) - f0(n, g)) - &
         nod(n, g)%dn(3) * (f0(p, g) + f0(n, g))
      end if
      if (j == xstag(i)%smin) then
         if (ysouth == 2) then
            jm = 0._dp
         else
            jm = -nod(n, g)%df(4) * f0(n, g) - nod(n, g)%dn(4) * f0(n, g)
         end if
      else
         jm = -nod(n, g)%df(4) * (f0(n, g) - f0(m, g)) - &
         nod(n, g)%dn(4) * (f0(n, g) + f0(m, g))
      end if


      l2 = (jp - jm) / ydel(j)

      ! z-direction zeroth transverse leakage
      if (k /= nzz) p = xyz(i, j, k + 1)
      if (k /= 1) m = xyz(i, j, k - 1)

      if (k == nzz) then
         if (ztop == 2) then
            jp = 0._dp
         else
            jp = nod(n, g)%df(5) * f0(n, g) - nod(n, g)%dn(5) * f0(n, g)
         end if
      else
         jp = -nod(n, g)%df(5) * (f0(p, g) - f0(n, g)) - &
         nod(n, g)%dn(5) * (f0(p, g) + f0(n, g))
      end if
      if (k == 1) then
         if (zbott == 2) then
            jm = 0._dp
         else
            jm = -nod(n, g)%df(6) * f0(n, g) - nod(n, g)%dn(6) * f0(n, g)
         end if
      else
         jm = -nod(n, g)%df(6) * (f0(n, g) - f0(m, g)) - &
         nod(n, g)%dn(6) * (f0(n, g) + f0(m, g))
      end if

      l3 = (jp - jm) / zdel(k)


   end subroutine lxyz

   !******************************************************************************!

   subroutine get_source ()

      use sdata, only: nnod, ng, exsrc

      ! Purpose:
      ! To update get source for the nodal update

      implicit none

      real(dp) :: l1, l2, l3
      integer :: g, n
      logical :: first = .true.

      if (first) then
         allocate(sx(nnod, ng), sy(nnod, ng), sz(nnod, ng))
         first = .false.
      end if

      do g = 1, ng
         do n = 1, nnod
            call lxyz(n, g, l1, l2, l3)

            if (cmode == 2) then
               sx(n, g) = l2 + l3 - exsrc(n, g)
               sy(n, g) = l1 + l3 - exsrc(n, g)
               sz(n, g) = l1 + l2 - exsrc(n, g)
            else
               sx(n, g) = l2 + l3
               sy(n, g) = l1 + l3
               sz(n, g) = l1 + l2
            end if
         end do
      end do

   end subroutine get_source

   !******************************************************************************!

   subroutine tlupd1 (u, n, g, lmom1)

      use sdata, only: xdel, ydel, zdel, xstag, ystag, nzz, &
      xwest, xeast, ynorth, ysouth, zbott, ztop, &
      ix, iy, iz, xyz, d

      ! Purpose:
      ! To calaculate transverse leakage first moments

      implicit none

      integer, intent(in) :: u, n, g
      real(dp), intent(out) :: lmom1

      real(dp) :: tm, tp
      real(dp) :: p1m, p2m, p1p, p2p, hp
      integer :: p, m, i, j, k

      ! Set i, j, k
      i = ix(n) ;
      j = iy(n) ;
      k = iz(n)

      if (u == 1) then
         ! Set paramaters for X-Direction Transverse leakage
         if (i /= ystag(j)%smax) p = xyz(i + 1, j, k)
         if (i /= ystag(j)%smin) m = xyz(i - 1, j, k)

         if (i == ystag(j)%smin) then
            if (xwest == 2) then
               tm = 1._dp
               tp = xdel(i + 1) / xdel(i)
               p1m = tm + 1._dp;
               p2m = 2._dp * tm + 1._dp
               p1p = tp + 1._dp
               hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
               lmom1 = (p1m * p2m * (sx(p, g) - sx(n, g))) / hp
            else
               tp = xdel(i + 1) / xdel(i)
               p1p = tp + 1._dp
               lmom1 = (sx(p, g) - sx(n, g)) / p1p
            end if
         else if (i == ystag(j)%smax) then
            if (xeast == 2) then
               tm = xdel(i - 1) / xdel(i)
               tp = 1._dp
               p1m = tm + 1._dp
               p1p = tp + 1._dp;
               p2p = 2._dp * tp + 1._dp
               hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
               lmom1 = (p1p * p2p * (sx(n, g) - sx(m, g))) / hp
            else
               tm = xdel(i - 1) / xdel(i)
               p1m = tm + 1._dp
               lmom1 = (sx(n, g) - sx(m, g)) / p1m
            end if
         else
            tm = xdel(i - 1) / xdel(i)
            tp = xdel(i + 1) / xdel(i)
            p1m = tm + 1._dp;
            p2m = 2._dp * tm + 1._dp
            p1p = tp + 1._dp;
            p2p = 2._dp * tp + 1._dp
            hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
            lmom1 = (p1m * p2m * (sx(p, g) - sx(n, g)) &
            + p1p * p2p * (sx(n, g) - sx(m, g))) / hp
         end if

         lmom1 = 0.25_dp * xdel(i)**2 / d(n, g) * lmom1

      else if (u == 2) then
         ! Set paramaters for Y-Direction Transverse leakage
         if (j /= xstag(i)%smax) p = xyz(i, j + 1, k)
         if (j /= xstag(i)%smin) m = xyz(i, j - 1, k)

         if (j == xstag(i)%smin) then
            if (ysouth == 2) then
               tm = 1._dp
               tp = ydel(j + 1) / ydel(j)
               p1m = tm + 1._dp;
               p2m = 2._dp * tm + 1._dp
               p1p = tp + 1._dp
               hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
               lmom1 = (p1m * p2m * (sy(p, g) - sy(n, g))) / hp
            else
               tp = ydel(j + 1) / ydel(j)
               p1p = tp + 1._dp
               lmom1 = (sy(p, g) - sy(n, g)) / p1p
            end if
         else if (j == xstag(i)%smax) then
            if (ynorth == 2) then
               tm = ydel(j - 1) / ydel(j)
               tp = 1._dp
               p1m = tm + 1._dp
               p1p = tp + 1._dp;
               p2p = 2._dp * tp + 1._dp
               hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
               lmom1 = (p1p * p2p * (sy(n, g) - sy(m, g))) / hp
            else
               tm = ydel(j - 1) / ydel(j)
               p1m = tm + 1._dp
               lmom1 = (sy(n, g) - sy(m, g)) / p1m
            end if
         else
            tm = ydel(j - 1) / ydel(j)
            tp = ydel(j + 1) / ydel(j)
            p1m = tm + 1._dp;
            p2m = 2._dp * tm + 1._dp
            p1p = tp + 1._dp;
            p2p = 2._dp * tp + 1._dp
            hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
            lmom1 = (p1m * p2m * (sy(p, g) - sy(n, g)) &
            + p1p * p2p * (sy(n, g) - sy(m, g))) / hp
         end if

         lmom1 = 0.25_dp * ydel(j)**2 / d(n, g) * lmom1

      else
         ! Set paramaters for Z-Direction Transverse leakage
         if (k /= nzz) p = xyz(i, j, k + 1)
         if (k /= 1) m = xyz(i, j, k - 1)

         if (k == 1) then
            if (zbott == 2) then
               tm = 1._dp
               tp = zdel(k + 1) / zdel(k)
               p1m = tm + 1._dp;
               p2m = 2._dp * tm + 1._dp
               p1p = tp + 1._dp
               hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
               lmom1 = (p1m * p2m * (sz(p, g) - sz(n, g))) / hp
            else
               tp = zdel(k + 1) / zdel(k)
               p1p = tp + 1._dp
               lmom1 = (sz(p, g) - sz(n, g)) / p1p
            end if
         else if (k == nzz) then
            if (ztop == 2) then
               tm = zdel(k - 1) / zdel(k)
               tp = 1._dp
               p1m = tm + 1._dp
               p1p = tp + 1._dp;
               p2p = 2._dp * tp + 1._dp
               hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
               lmom1 = (p1p * p2p * (sz(n, g) - sz(m, g))) / hp
            else
               tm = zdel(k - 1) / zdel(k)
               p1m = tm + 1._dp
               lmom1 = (sz(n, g) - sz(m, g)) / p1m
            end if
         else
            tm = zdel(k - 1) / zdel(k)
            tp = zdel(k + 1) / zdel(k)
            p1m = tm + 1._dp;
            p2m = 2._dp * tm + 1._dp
            p1p = tp + 1._dp;
            p2p = 2._dp * tp + 1._dp
            hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
            lmom1 = (p1m * p2m * (sz(p, g) - sz(n, g)) &
            + p1p * p2p * (sz(n, g) - sz(m, g))) / hp
         end if

         lmom1 = 0.25_dp * zdel(k)**2 / d(n, g) * lmom1

      end if

   end subroutine tlupd1

   !****************************************************************************!

   subroutine tlupd2 (u, n, g, lmom2)

      use sdata, only: xdel, ydel, zdel, xstag, ystag, nzz, &
      xwest, xeast, ynorth, ysouth, zbott, ztop, &
      ix, iy, iz, xyz, d

      ! Purpose:
      ! To calaculate transverse leakage second moments

      implicit none

      integer, intent(in) :: u, n, g
      real(dp), intent(out) :: lmom2

      real(dp) :: tm, tp
      real(dp) :: p1m, p1p, hp
      integer :: p, m, i, j, k

      ! Set i, j, k
      i = ix(n) ;
      j = iy(n) ;
      k = iz(n)

      if (u == 1) then
         ! Set paramaters for X-Direction Transverse leakage
         if (i /= ystag(j)%smax) p = xyz(i + 1, j, k)
         if (i /= ystag(j)%smin) m = xyz(i - 1, j, k)

         if (i == ystag(j)%smin) then
            if (xwest == 2) then
               tm = 1._dp
               tp = xdel(i + 1) / xdel(i)
               p1m = tm + 1._dp
               p1p = tp + 1._dp
               hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
               lmom2 = (p1m * (sx(p, g) - sx(n, g))) / hp
            else
               lmom2 = 0._dp
            end if
         else if (i == ystag(j)%smax) then
            if (xeast == 2) then
               tm = xdel(i - 1) / xdel(i)
               tp = 1._dp
               p1m = tm + 1._dp
               p1p = tp + 1._dp
               hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
               lmom2 = (p1p * (sx(m, g) - sx(n, g))) / hp
            else
               lmom2 = 0._dp
            end if
         else
            tm = xdel(i - 1) / xdel(i)
            tp = xdel(i + 1) / xdel(i)
            p1m = tm + 1._dp
            p1p = tp + 1._dp
            hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
            lmom2 = (p1m * (sx(p, g) - sx(n, g)) + p1p * (sx(m, g) - sx(n, g))) / hp
         end if

         lmom2 = 0.25_dp * xdel(i)**2 / d(n, g) * lmom2

      else if (u == 2) then
         ! Set paramaters for Y-Direction Transverse leakage
         if (j /= xstag(i)%smax) p = xyz(i, j + 1, k)
         if (j /= xstag(i)%smin) m = xyz(i, j - 1, k)

         if (j == xstag(i)%smin) then
            if (ysouth == 2) then
               tm = 1._dp
               tp = ydel(j + 1) / ydel(j)
               p1m = tm + 1._dp
               p1p = tp + 1._dp
               hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
               lmom2 = (p1m * (sy(p, g) - sy(n, g))) / hp
            else
               lmom2 = 0._dp
            end if
         else if (j == xstag(i)%smax) then
            if (ynorth == 2) then
               tm = ydel(j - 1) / ydel(j)
               tp = 1._dp
               p1m = tm + 1._dp
               p1p = tp + 1._dp
               hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
               lmom2 = (p1p * (sy(m, g) - sy(n, g))) / hp
            else
               lmom2 = 0._dp
            end if
         else
            tm = ydel(j - 1) / ydel(j)
            tp = ydel(j + 1) / ydel(j)
            p1m = tm + 1._dp
            p1p = tp + 1._dp
            hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
            lmom2 = (p1m * (sy(p, g) - sy(n, g)) + p1p * (sy(m, g) - sy(n, g))) / hp
         end if

         lmom2 = 0.25_dp * ydel(j)**2 / d(n, g) * lmom2

      else
         ! Set paramaters for Z-Direction Transverse leakage
         if (k /= nzz) p = xyz(i, j, k + 1)
         if (k /= 1) m = xyz(i, j, k - 1)

         if (k == 1) then
            if (zbott == 2) then
               tm = 1._dp
               tp = zdel(k + 1) / zdel(k)
               p1m = tm + 1._dp
               p1p = tp + 1._dp
               hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
               lmom2 = (p1m * (sz(p, g) - sz(n, g))) / hp
            else
               lmom2 = 0._dp
            end if
         else if (k == nzz) then
            if (ztop == 2) then
               tm = zdel(k - 1) / zdel(k)
               tp = 1._dp
               p1m = tm + 1._dp
               p1p = tp + 1._dp
               hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
               lmom2 = (p1p * (sz(m, g) - sz(n, g))) / hp
            else
               lmom2 = 0._dp
            end if
         else
            tm = zdel(k - 1) / zdel(k)
            tp = zdel(k + 1) / zdel(k)
            p1m = tm + 1._dp
            p1p = tp + 1._dp
            hp = 2._dp * p1m * p1p * (tm + tp + 1._dp)
            lmom2 = (p1m * (sz(p, g) - sz(n, g)) + p1p * (sz(m, g) - sz(n, g))) / hp
         end if

         lmom2 = 0.25_dp * zdel(k)**2 / d(n, g) * lmom2

      end if


   end subroutine tlupd2

   !****************************************************************************!

   function get_b(u, n) result (b)

      !Purpose: To calculate Buckling for node n and direction u


      use sdata, only: ng, xdel, ydel, zdel, ix, iy, iz, ke, &
      sigr, d, chi, mat, nuf, sigs, ke, tbeta, dfis
      implicit none

      integer, intent(in) :: u, n      ! direction and node number
      real(dp), dimension(ng, ng) :: b        ! Buckling B2

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
                  dum = sigr(n, g) - chi(mat(n), g) * nuf(n, h) / ke
               else
                  dum = -sigs(n, h, g) - chi(mat(n), g) * nuf(n, h) / ke
               end if
               b(g, h) = 0.25_dp * dn**2 / d(n, g) * dum
            end do
         end do
      else if (cmode == 2) then
         do g = 1, ng
            do h = 1, ng
               if (g == h) then
                  dum = sigr(n, g) - (1._dp - tbeta(mat(n)) + dfis(n)) &
                  * chi(mat(n), g) * chi(mat(n), g) * nuf(n, h)
               else
                  dum = -sigs(n, h, g) - (1._dp - tbeta(mat(n)) + dfis(n)) &
                  * chi(mat(n), g) * nuf(n, h)
               end if
               b(g, h) = 0.25_dp * dn**2 / d(n, g) * dum
            end do
         end do
      else
         do g = 1, ng
            do h = 1, ng
               if (g == h) then
                  dum = sigr(n, g) - chi(mat(n), g) * nuf(n, h) / ke
               else
                  dum = -sigs(n, g, h) - chi(mat(n), h) * nuf(n, g) / ke
               end if
               b(g, h) = 0.25_dp * dn**2 / d(n, g) * dum
            end do
         end do
      end if

   end function get_b

   !****************************************************************************!

   subroutine get_abefgh (n, u, aa, bb, ee, ff, gg, hh)

      use sdata, only: ng, xdel, ydel, zdel, ix, iy, iz, d, sigr

      ! Purpose:
      ! To calaculate A,B,E,F,G,H paramters used to calculate matrix elements for
      ! semi-analytic nodal update

      implicit none

      integer, intent(in) :: u, n
      real(dp), dimension(:), intent(out) :: aa, bb, ee, ff, gg, hh

      real(dp) :: dn, alp, alp2
      real(dp) :: m1s, m0c, m2c
      integer :: g

      if (u == 1) then
         dn = xdel(ix(n))
      else if (u == 2) then
         dn = ydel(iy(n))
      else
         dn = zdel(iz(n))
      end if

      do g = 1, ng
         !Calculate alpha and othe paramters to calculate A,B,E,F,G,H
         alp = 0.5_dp * sqrt(sigr(n, g) / d(n, g)) * dn
         alp2 = alp**2
         m0c = sinh(alp) / alp
         m1s = 3._dp * (cosh(alp) / alp - sinh(alp) / alp2)
         m2c = 5._dp * (sinh(alp) / alp - 3._dp * cosh(alp) / alp2 &
         + 3._dp * sinh(alp) / alp**3)

         aa(g) = (sinh(alp) - m1s) / (alp2 * m1s)
         bb(g) = (cosh(alp) - m0c - m2c) / (alp2 * m2c)
         ee(g) = (m0c / m2c - 3._dp / alp2)
         ff(g) = (alp * cosh(alp) - m1s) / (alp2 * m1s)
         gg(g) = (alp * sinh(alp) - 3._dp * m2c) / (cosh(alp) - m0c - m2c)
         hh(g) = (alp * cosh(alp) - m1s) / (sinh(alp) - m1s)
      end do

   end subroutine get_abefgh

end module nodal
