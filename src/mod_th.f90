module th

   use sdata, only: dp

   implicit none

   save

   contains

   subroutine th_iter(maxiter, oupt)
      !
      ! Purpose:
      !    To do thermal-hydrailics iteration
      !

      use sdata, only: nnod, ftem, mtem, cden, bcon, bpos, npow, pow, ppow, &
      zdel, node_nf, ix, iy, iz, th_err, node_nf, ix, iy, iz, &
      th_niter, fer, ferc, ser, serc, get_time, th_time, ke
      use cmfd, only: outer_th
      use io, only: ounit, scr
      use xsec, only: xs_updt

      implicit none

      integer, intent(in) :: maxiter    ! maximum number of th iteration
      integer, intent(in) :: oupt       ! print ouput option
      real(dp), dimension(nnod) :: pline
      real(dp), dimension(nnod) :: otem
      integer :: n, i
      real(dp) :: st, fn

      th_err = 1.
      do i = 1, maxiter
         ! Save old fuel temp
         otem = ftem

         call xs_updt(bcon, ftem, mtem, cden, bpos)

         ! Perform outer inner iteration
         call outer_th()

         !Get start th_time
         st = get_time()

         ! Calculate power density
         call powdis(npow)

         ! Calculate linear power density for each nodes (W/cm)
         do n = 1, nnod
            pline(n) = npow(n) * pow * ppow * 0.01_dp &
            / (node_nf(ix(n), iy(n)) * zdel(iz(n)))     ! Linear power density (W/cm)
         end do

         ! Update fuel, moderator temp. and coolant density
         call th_upd(pline)

         ! Get fuel absolute difference from current and previous th iteration
         call abse(ftem, otem, th_err)

         !Get th_time
         fn = get_time()
         th_time = th_time + (fn - st)

         if (oupt > 0) then
            write(ounit, '(I5,F13.6,3ES15.5)') i, ke, ser, fer, th_err        ! Write outer iteration evolution
            if (scr) write(*, '(I5,F13.6,3ES15.5)') i, ke, ser, fer, th_err   ! Write outer iteration evolution
         end if

         ! If error is small enough
         if (th_err < 0.01 .and. fer < ferc .and. ser < serc) exit

      end do

      if (i - 1 == th_niter) then
         write(ounit, *) '  MAXIMUM TH ITERATION REACHED.'
         write(ounit, *) '  CALCULATION MIGHT BE NOT CONVERGED OR CHANGE ITERATION CONTROL'
         write(*, *) '  MAXIMUM TH ITERATION REACHED.'
         write(*, *) '  CALCULATION MIGHT BE NOT CONVERGED OR CHANGE ITERATION CONTROL'
         stop
      end if



   end subroutine th_iter

   !****************************************************************************!

   subroutine powdis (p)

      !
      ! Purpose:
      !    To calculate power for each nodes
      !


      use sdata, only: ng, nnod, sigf, f0, vdel, powtot, mode
      use io, only: ounit

      implicit none

      real(dp), dimension(:), intent(out) :: p
      integer :: g, n
      real(dp) :: pow

      p = 0._dp
      do g = 1, ng
         do n = 1, nnod
            pow = f0(n, g) * sigf(n, g) * vdel(n)
            if (pow < 0.) pow = 0.
            p(n) = p(n) + pow
         end do
      end do

      ! Normalize to 1._DP
      powtot = 0._dp
      do n = 1, nnod
         powtot = powtot + p(n)
      end do

      if (powtot <= 0 .and. mode /= 'FIXEDSRC') then
         write(ounit, *) '   ERROR: TOTAL NODES POWER IS ZERO OR LESS'
         write(ounit, *) '   STOP IN subroutine POWDIS'
         stop
      end if

      if (powtot > 0.0) then
         do n = 1, nnod
            p(n) = p(n) / powtot
         end do
      end if

   end subroutine powdis

   !****************************************************************************!

   subroutine abse(newf, oldf, rel)

      !
      ! Purpose:
      !    To calculate Max Relative error

      use sdata, only: nnod

      implicit none

      real(dp), dimension(:), intent(in) :: newf, oldf
      real(dp), intent(out) :: rel

      real(dp) :: error
      integer :: n

      rel = 0.

      do n = 1, nnod
         if (abs(newf(n)) > 1.e-10_dp) then
            error = abs(newf(n) - oldf(n))
            if (error > rel) rel = error
         end if
      end do

   end subroutine abse

   !****************************************************************************!

   subroutine par_ave(par, ave)
      !
      ! Purpose:
      !    To calculate average fuel temp (only for active core)
      !

      use sdata, only: vdel, nnod, ng, nuf

      implicit none

      real(dp), dimension(:), intent(in) :: par
      real(dp), intent(out) :: ave
      real(dp) :: dum, dum2
      integer :: n

      dum = 0. ;
      dum2 = 0.
      do n = 1, nnod
         if (nuf(n, ng) > 0.) then
            dum = dum + par(n) * vdel(n)
            dum2 = dum2 + vdel(n)
         end if
      end do

      ave = dum / dum2

   end subroutine par_ave

   !****************************************************************************!

   subroutine par_ave_out(par, ave)
      !
      ! Purpose:
      !    To calculate average fuel temp (only for active core)
      !

      use sdata, only: vdel, nnod, iz, nzz, nuf, ng, ix, iy, xyz, nxx, nyy, nzz, ystag

      implicit none

      real(dp), dimension(:), intent(in) :: par
      real(dp), intent(out) :: ave
      real(dp) :: dum, dum2
      integer, dimension(nxx, nyy) :: zmax
      integer :: n, i, j, k

      ! get number of nodex in axial direction from bottom -> fuel
      do j = 1, nyy
         do i = ystag(j)%smin, ystag(j)%smax
            zmax(i, j) = 0
            do k = 1, nzz / 2
               if (nuf(xyz(i, j, k), ng) < 1.e-5_dp) zmax(i, j) = zmax(i, j) + 1
            end do
         end do
      end do
      ! get number of nodex in axial direction from fuel -> top reflectors
      do j = 1, nyy
         do i = ystag(j)%smin, ystag(j)%smax
            do k = 1, nzz
               if (nuf(xyz(i, j, k), ng) > 1.e-5) zmax(i, j) = zmax(i, j) + 1
            end do
         end do
      end do

      dum = 0. ;
      dum2 = 0.
      do n = 1, nnod
         if (iz(n) == zmax(ix(n), iy(n)) .and. nuf(n, ng) > 1.e-5) then
            dum = dum + par(n) * vdel(n)
            dum2 = dum2 + vdel(n)
         end if
      end do

      ave = dum / dum2

   end subroutine par_ave_out

   !****************************************************************************!

   subroutine par_max(par, pmax)
      !
      ! Purpose:
      !    To calculate maximum fuel tem, coolant tem, and density
      !

      use sdata, only: nnod

      implicit none

      real(dp), dimension(:), intent(in) :: par
      real(dp), intent(out) :: pmax
      integer :: n

      pmax = 0.
      do n = 1, nnod
         if (par(n) > pmax) pmax = par(n)
      end do

   end subroutine par_max

   !****************************************************************************!

   subroutine getent(t, ent)
      !
      ! Purpose:
      !    To get enthalpy for given coolant temp. from steam table
      !

      use sdata, only: stab, ntem
      use io, only : ounit

      implicit none

      real(dp), intent(in) :: t
      real(dp), intent(out) :: ent
      real(dp) :: t1, ent1
      real(dp) :: t2, ent2
      integer :: i

      if ((t < stab(1, 1)) .or. (t > stab(ntem, 1))) then
         write(ounit, *) '  Coolant temp. : ', t, 'K'
         write(ounit, *) '  ERROR : MODERATOR TEMP. IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
         write(ounit, *) '  CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
         write(*, *) '  Coolant temp. : ', t, 'K'
         write(*, *) '  ERROR : MODERATOR TEMP. IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
         write(*, *) '  CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
         stop
      end if

      t2 = stab(1, 1) ;
      ent2 = stab(1, 3)
      do i = 2, ntem
         t1 = t2
         ent1 = ent2
         t2 = stab(i, 1) ;
         ent2 = stab(i, 3)
         if ((t >= t1) .and. (t <= t2)) then
            ent = ent1 + (t - t1) / (t2 - t1) * (ent2 - ent1)
            exit
         end if
      end do


   end subroutine getent

   !****************************************************************************!

   subroutine gettd(ent, t, rho, prx, kvx, tcx, rx)
      !
      ! Purpose:
      !    To get enthalpy for given coolant temp. from steam table
      !

      use sdata, only: stab, ntem
      use io, only : ounit

      implicit none

      real(dp), intent(in) :: ent
      real(dp), intent(out) :: t, rho, prx, kvx, tcx
      real(dp), intent(out), optional :: rx
      real(dp) :: ratx

      integer :: i, i1, i2

      ! Get two closest interpolation points
      if (ent >= stab(1, 3) .and. ent <= stab(ntem, 3)) then  !If enthalpy inside data range
         do i = 2, ntem
            if (ent >= stab(i - 1, 3) .and. ent <= stab(i, 3)) then
               i1 = i - 1
               i2 = i
               exit
            end if
         end do
      else if (ent < stab(1, 3) .and. (stab(1, 3) - ent) / stab(1, 3) < 0.1) then !If 10% lower than min. steam table data
         i1 = 1
         i2 = 2
      else if (ent > stab(ntem, 3) .and. (ent - stab(ntem, 3)) / stab(ntem, 3) < 0.1) then !if 10% higher than max. steam table data
         i1 = ntem - 1
         i2 = ntem
      else
         write(ounit, 1557) ent / 1000.
         write(*, 1557) ent / 1000.
         write(ounit, *) '   CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
         write(*, *) '   CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
         stop
         1557 format(2x, '  ERROR: ENTHALPY', f8.1, ' KJ/Kg &
         & is out of the range in the steam table')
      end if

      ! Interpolate
      ratx = (ent - stab(i1, 3)) / (stab(i2, 3) - stab(i1, 3))
      t = stab(i1, 1) + ratx * (stab(i2, 1) - stab(i1, 1))
      rho = stab(i1, 2) + ratx * (stab(i2, 2) - stab(i1, 2))
      prx = stab(i1, 4) + ratx * (stab(i2, 4) - stab(i1, 4))
      kvx = stab(i1, 5) + ratx * (stab(i2, 5) - stab(i1, 5))
      tcx = stab(i1, 6) + ratx * (stab(i2, 6) - stab(i1, 6))
      if (present(rx)) then
         rx = 1000._dp * (stab(i2, 2) - stab(i1, 2)) / (stab(i2, 3) - stab(i1, 3))
      end if

   end subroutine gettd

   !****************************************************************************!

   real(dp) function getkc(t)
      !
      ! Purpose:
      !    To calculate thermal conductivity of cladding
      !

      implicit none

      real(dp), intent(in) :: t

      getkc = 7.51_dp + 2.09e-2_dp * t - 1.45e-5_dp * t**2 + 7.67e-9_dp * t**3

   end function getkc

   !****************************************************************************!

   real(dp) function getkf(t)
      !
      ! Purpose:
      !    To calculate thermal conductivity of fuel
      !

      implicit none

      real(dp), intent(in) :: t

      getkf = 1.05_dp + 2150.0_dp / (t - 73.15_dp)

   end function getkf

   !****************************************************************************!

   real(dp) function getcpc(t)
      !
      ! Purpose:
      !    To calculate specific heat capacity of cladding
      !

      implicit none

      real(dp), intent(in) :: t

      getcpc = 252.54_dp + 0.11474_dp * t

   end function getcpc

   !****************************************************************************!

   real(dp) function getcpf(t)
      !
      ! Purpose:
      !    To calculate specific heat capacity of fuel
      !

      implicit none

      real(dp), intent(in) :: t

      getcpf = 162.3_dp + 0.3038_dp * t - 2.391e-4_dp * t**2 + 6.404e-8_dp * t**3

   end function getcpf

   !****************************************************************************!

   subroutine tridiasolve(a, b, c, d, x)
      !
      ! Purpose:
      !    To solve tridiagonal matrix
      !

      implicit none

      real(dp), dimension(:), intent(inout) :: a, b, c, d
      real(dp), dimension(:), intent(out) :: x

      integer :: i, n

      n = size(d)

      ! Gauss Elimination
      c(1) = c(1) / b(1)
      d(1) = d(1) / b(1)
      do i = 2, n
         c(i) = c(i) / (b(i) - a(i) * c(i - 1))
         d(i) = (d(i) - a(i) * d(i - 1)) / (b(i) - a(i) * c(i - 1))
      end do

      ! Back Substitution
      x(n) = d(n)
      do i = n - 1, 1, -1
         x(i) = d(i) - c(i) * x(i + 1)
      end do

   end subroutine tridiasolve

   !****************************************************************************!

   real(dp) function geths(xden, tc, kv, pr)
      !
      ! Purpose:
      !    To calculate heat transfer coef.
      !

      use sdata, only: dh, farea, cflow

      implicit none

      real(dp), intent(in) :: xden  ! coolant densisty
      real(dp), intent(in) :: tc  ! coolant thermal conductivity
      real(dp), intent(in) :: kv  ! kinematic viscosity
      real(dp), intent(in) :: pr  ! Prandtl Number

      real(dp) :: cvelo, nu, re

      cvelo = cflow / (farea * xden * 1000._dp)        ! Calculate flow velocity (m/s)
      re = cvelo * dh / (kv * 1.e-6_dp)                 ! Calculate Reynolds Number
      nu = 0.023_dp * (pr**0.4_dp) * (re**0.8_dp)                ! Calculate Nusselt Number
      geths = (tc / dh) * nu                        ! Calculate heat transfer coefficient


   end function geths

   !****************************************************************************!

   subroutine th_trans(xpline, h)

      !
      ! Purpose:
      !    To perform fuel pin thermal transient
      !

      use sdata, only: nnod, mtem, cden, ftem, tin, cflow, nxx, nyy, cf, ent, heatf, &
      tfm, nt, rpos, rdel, rf, rg, rc, farea, dia, pi, zdel, &
      ix, iy, iz, frate, th_time, get_time

      implicit none

      real(dp), dimension(:), intent(in) :: xpline    ! Linear Power Density (W/cm)
      real(dp), intent(in) :: h                       ! Time step

      integer :: i, j, k, n
      real(dp), dimension(nt + 1) :: a, b, c, d
      real(dp) :: hs, hg = 1.d4, kt, kt1, kt2          ! coolant heat transfer coef., gap heat transfer coef, and thermal conductivity
      real(dp) :: alpha = 0.7_dp
      real(dp) :: xa, xc
      real(dp) :: fdens = 10.412e3            ! UO2 density (kg/m3)
      real(dp) :: cdens = 6.6e3               ! Cladding density (kg/m3)
      real(dp) :: cp                          ! Specific heat capacity
      real(dp) :: eps, eta
      real(dp) :: mdens, vol                  ! Coolant density and channel volume
      real(dp), dimension(nnod) :: entp        ! previous enthalpy

      real(dp) :: pdens      ! power densisty  (W/m3)
      real(dp) :: enti       ! Coolant inlet enthalpy
      real(dp), dimension(nxx, nyy) :: entm, bfrate
      real(dp) :: cpline     ! Coolant Linear power densisty (W/m)
      real(dp) :: pr, kv, tcon ! Coolant Prandtl Number, Kinematic viscosity, and thermal conductivity
      real(dp) :: r
      logical :: first = .true.
      real(dp) :: st, fn

      !Get start th_time
      st = get_time()

      !set initial tridiagonal matrix element a, b
      a = 0._dp;
      b = 0._dp;
      c = 0._dp;

      if (first) then
         allocate(frate(nnod))
         frate = cflow
         first = .false.
      end if
      ! frate = cflow
      call getent(tin, enti)
      entp = ent

      do n = 1, nnod
         mdens = cden(n) * 1000._dp                                    ! Coolant density (kg/m3)
         cpline = heatf(n) * pi * dia + cf * xpline(n) * 100._dp       ! Coolant Linear power densisty (W/m)
         vol = farea * zdel(iz(n)) * 0.01_dp                             ! channel node volume
         i = ix(n) ;
         j = iy(n) ;
         k = iz(n)

         if (k == 1) then                                              ! Calculate coolant enthalpy
            eps = mdens * vol / h
            ent(n) = (cpline * zdel(iz(n)) * 0.01_dp + 2._dp * frate(n) * enti &
            + eps * entp(n)) / (eps + 2._dp * frate(n))                             ! Calculate enthalpy
            call gettd(ent(n), mtem(n), cden(n), pr, kv, tcon, r)        ! Get corresponding temp and density
            entm(i, j) = 2._dp * ent(n) - enti
            frate(n) = cflow - 0.5_dp * vol / h * r * (ent(n) - entp(n))
            bfrate(i, j) = 2._dp * frate(n) - cflow
         else
            eps = mdens * vol / h
            ent(n) = (cpline * zdel(iz(n)) * 0.01_dp + 2._dp * frate(n) &
            * entm(ix(n), iy(n)) + eps * entp(n)) / (eps + 2._dp * frate(n))
            call gettd(ent(n), mtem(n), cden(n), pr, kv, tcon, r)
            entm(i, j) = 2._dp * ent(n) - entm(i, j)
            frate(n) = bfrate(i, j) - 0.5_dp * vol / h * r * (ent(n) - entp(n))
            bfrate(i, j) = 2._dp * frate(n) - bfrate(i, j)
         end if

         hs = geths(cden(n), pr, kv, tcon)                                               ! Calculate heat transfer coef
         pdens = 100._dp * xpline(n) / (pi * rf**2)                ! Fuel pin Power Density (W/m3)

         ! Calculate tridiagonal matrix: a, b, c and source: d
         ! For nt=1 [FUEL CENTERLINE]
         kt1 = getkf(tfm(n, 1))                                                     ! Get thermal conductivity
         kt2 = getkf(tfm(n, 2))
         kt = 2._dp * kt1 * kt2 / (kt1 + kt2)
         cp = getcpf(tfm(n, 1))                                                           ! Get specific heat capacity
         eta = fdens * cp * rpos(1)**2 / (2._dp * h)
         xc = kt * rpos(1) / rdel(1)
         b(1) = xc + eta
         c(1) = -xc
         d(1) = pdens * 0.5_dp * rpos(1)**2 + eta * tfm(n, 1)

         do i = 2, nt - 2
            kt1 = kt2
            kt2 = getkf(tfm(n, i + 1))
            kt = 2._dp * kt1 * kt2 / (kt1 + kt2)
            cp = getcpf(tfm(n, i))
            eta = fdens * cp * (rpos(i)**2 - rpos(i - 1)**2) / (2. * h)
            xa = xc
            xc = kt * rpos(i) / rdel(i)
            a(i) = -xa
            b(i) = xa + xc + eta
            c(i) = -xc
            d(i) = pdens * 0.5_dp * (rpos(i)**2 - rpos(i - 1)**2) + eta * tfm(n, i)
         end do

         ! For nt-1 [FUEL-GAP INTERFACE]
         cp = getcpf(tfm(n, nt - 1))
         eta = fdens * cp * (rf**2 - rpos(nt - 2)**2) / (2. * h)
         xa = xc
         xc = rg * hg
         a(nt - 1) = -xa
         b(nt - 1) = xa + xc + eta
         c(nt - 1) = -xc
         d(nt - 1) = pdens * 0.5_dp * (rf**2 - rpos(nt - 2)**2) + eta * tfm(n, nt - 1)

         ! For nt [GAP-CLADDING INTERFACE]
         kt1 = getkc(tfm(n, nt))
         kt2 = getkc(tfm(n, nt + 1))
         kt = 2._dp * kt1 * kt2 / (kt1 + kt2)     ! For cladding
         cp = getcpc(tfm(n, nt))
         eta = cdens * cp * (rpos(nt)**2 - rg**2) / (2. * h)
         xa = xc
         xc = kt * rpos(nt) / rdel(nt)
         a(nt) = -xa
         b(nt) = xa + xc + eta
         c(nt) = -xc
         d(nt) = eta * tfm(n, nt)

         ! For nt+1  [CLADDING-COOLANT INTERFACE]
         cp = getcpc(tfm(n, nt + 1))
         eta = cdens * cp * (rc**2 - rpos(nt)**2) / (2. * h)
         xa = xc
         xc = rc * hs
         a(nt + 1) = -xa
         b(nt + 1) = xa + xc + eta
         d(nt + 1) = rc * hs * mtem(n) + eta * tfm(n, nt + 1)

         ! Solve tridiagonal matrix
         call tridiasolve(a, b, c, d, tfm(n, :))

         ! Get lumped fuel temp
         ftem(n) = (1. - alpha) * tfm(n, 1) + alpha * tfm(n, nt - 1)

         ! Calculate heat flux
         heatf(n) = hs * (tfm(n, nt + 1) - mtem(n))
      end do

      !Get th_time
      fn = get_time()
      th_time = th_time + (fn - st)

   end subroutine th_trans

   !****************************************************************************!

   subroutine th_upd(xpline)

      !
      ! Purpose:
      !    To update thermal parameters
      !

      use sdata, only: nnod, mtem, cden, ftem, tin, ix, iy, iz, nxx, nyy, cflow, cf, &
      ent, heatf, tfm, nt, rpos, rdel, rf, rg, rc, pi, zdel, dia

      implicit none

      real(dp), dimension(:), intent(in) :: xpline    ! Linear Power Density (W/cm)

      integer :: i, n
      real(dp), dimension(nt + 1) :: a, b, c, d
      real(dp) :: hs, hg = 1.d4, kt, kt1, kt2
      real(dp) :: alp = 0.7_dp
      real(dp) :: xa, xc
      real(dp) :: pdens      ! power densisty  (W/m3)
      real(dp) :: enti       ! Coolant inlet enthalpy
      real(dp), dimension(nxx, nyy) :: entm   ! enthalpy at node boundary
      real(dp) :: cpline     ! Coolant Linear power densisty (W/m)
      real(dp) :: pr, kv, tcon ! Coolant Prandtl Number, Kinematic viscosity, and thermal conductivity
      real(dp) :: zd  ! zdel in meter

      !set initial tridiagonal matrix element a, b
      a = 0._dp;
      b = 0._dp;
      c = 0._dp;

      call getent(tin, enti)

      do n = 1, nnod
         cpline = heatf(n) * pi * dia + cf * xpline(n) * 100._dp          ! Coolant Linear power densisty (W/m)
         zd = zdel(iz(n)) * 0.01_dp
         if (iz(n) == 1) then                                              ! For most bootom channel
            ent(n) = enti + 0.5_dp * cpline * zd / cflow    ! Calculate coolant enthalpy
            call gettd(ent(n), mtem(n), cden(n), pr, kv, tcon)             ! Get corresponding temp and density
            entm(ix(n), iy(n)) = 2._dp * ent(n) - enti                      ! Extrapolate enthalpy at node boundary
         else
            ent(n) = entm(ix(n), iy(n)) + 0.5_dp * cpline * zd / cflow
            call gettd(ent(n), mtem(n), cden(n), pr, kv, tcon)
            entm(ix(n), iy(n)) = 2._dp * ent(n) - entm(ix(n), iy(n))
         end if

         hs = geths(cden(n), pr, kv, tcon)
         pdens = (1._dp - cf) * 100._dp * xpline(n) / (pi * rf**2)        ! Fuel pin Power Density (W/m3)

         ! Calculate tridiagonal matrix: a, b, c and source: d
         ! For nt=1 [FUEL CENTERLINE]
         kt1 = getkf(tfm(n, 1))                                         ! Get thermal conductivity
         kt2 = getkf(tfm(n, 2))
         kt = 2._dp * kt1 * kt2 / (kt1 + kt2)
         xc = kt * rpos(1) / rdel(1)
         b(1) = xc
         c(1) = -xc
         d(1) = pdens * 0.5_dp * rpos(1)**2

         do i = 2, nt - 2
            kt1 = kt2
            kt2 = getkf(tfm(n, i + 1))
            kt = 2._dp * kt1 * kt2 / (kt1 + kt2)
            xa = xc
            xc = kt * rpos(i) / rdel(i)
            a(i) = -xa
            b(i) = xa + xc
            c(i) = -xc
            d(i) = pdens * 0.5_dp * (rpos(i)**2 - rpos(i - 1)**2)
         end do

         ! For nt-1 [FUEL-GAP INTERFACE]
         xa = xc
         xc = rg * hg
         a(nt - 1) = -xa
         b(nt - 1) = xa + xc
         c(nt - 1) = -xc
         d(nt - 1) = pdens * 0.5_dp * (rf**2 - rpos(nt - 2)**2)

         ! For nt [GAP-CLADDING INTERFACE]
         kt1 = getkc(tfm(n, nt))
         kt2 = getkc(tfm(n, nt + 1))
         kt = 2._dp * kt1 * kt2 / (kt1 + kt2)     ! For cladding
         xa = xc
         xc = kt * rpos(nt) / rdel(nt)
         a(nt) = -xa
         b(nt) = xa + xc
         c(nt) = -xc
         d(nt) = 0.

         ! For nt+1  [CLADDING-COOLANT INTERFACE]
         xa = xc
         a(nt + 1) = -xa
         b(nt + 1) = xa + hs * rc
         d(nt + 1) = rc * hs * mtem(n)

         ! Solve tridiagonal matrix
         call tridiasolve(a, b, c, d, tfm(n, :))

         ! Get lumped fuel temp
         ftem(n) = (1. - alp) * tfm(n, 1) + alp * tfm(n, nt - 1)

         ! Calculate heat flux
         heatf(n) = hs * (tfm(n, nt + 1) - mtem(n))
      end do


   end subroutine th_upd

end module th
