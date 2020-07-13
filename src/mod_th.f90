MODULE th

USE sdata, ONLY: DP

IMPLICIT NONE

SAVE

CONTAINS

SUBROUTINE th_iter(ind)
!
! Purpose:
!    To do thermal-hydrailics iteration
!

USE sdata, ONLY: nnod, ftem, mtem, cden, bcon, bpos, npow, pow, ppow,  &
                 zdel, node_nf, ix, iy, iz, th_err, node_nf, ix, iy, iz, &
                 th_niter, nth, fer, ferc, ser, serc, get_time, th_time
USE cmfd, ONLY: outer_th, PowDis
USE io, ONLY: bxtab, ounit
USE xsec, ONLY: XS_updt, XStab_updt

IMPLICIT NONE

INTEGER, INTENT(IN), OPTIONAL :: ind    ! if iteration reaching th_iter and ind = 0 then STOP
REAL(DP), DIMENSION(nnod) :: pline
REAL(DP), DIMENSION(nnod) :: otem
INTEGER :: mx_iter, n, l
REAL(DP) :: st, fn

! Determine maximum iteration
IF (PRESENT(ind)) THEN
  mx_iter = th_niter
ELSE
  mx_iter = 2
END IF

th_err = 1.
DO l = 1, mx_iter
    ! Save old fuel temp
    otem = ftem

    ! Update XS
    IF (bxtab == 1) THEN
      CALL XStab_updt(bcon, ftem, mtem, cden, bpos)
    ELSE
      CALL XS_updt(bcon, ftem, mtem, cden, bpos)
    END IF

    ! Perform outer inner iteration
    CALL outer_th(nth)

    !Get start th_time
    st = get_time()

    ! Calculate power density
    CALL PowDis(npow)

    ! Calculate linear power density for each nodes (W/cm)
    DO n = 1, nnod
        pline(n) = npow(n) * pow * ppow * 0.01_DP &
                 / (node_nf(ix(n),iy(n)) * zdel(iz(n)))     ! Linear power density (W/cm)
    END DO

    ! Update fuel, moderator temp. and coolant density
    CALL th_upd(pline)

    ! Get fuel absolute difference from current and previous th iteration
    CALL AbsE(ftem, otem, th_err)

    !Get th_time
    fn = get_time()
    th_time = th_time + (fn-st)

    ! If error is small enough
    IF (th_err < 0.01 .and. fer < ferc .and. ser < serc .and. present(ind)) EXIT

END DO

IF (PRESENT(ind) .AND. l-1 == th_niter) THEN
  WRITE(ounit,*) '  MAXIMUM TH ITERATION REACHED.'
  WRITE(ounit,*) '  CALCULATION MIGHT BE NOT CONVERGED OR CHANGE ITERATION CONTROL'
  WRITE(*,*) '  MAXIMUM TH ITERATION REACHED.'
  WRITE(*,*) '  CALCULATION MIGHT BE NOT CONVERGED OR CHANGE ITERATION CONTROL'
  STOP
END IF



END SUBROUTINE th_iter


SUBROUTINE AbsE(newF, oldF, rel)

  !
  ! Purpose:
  !    To calculate Max Relative error

USE sdata, ONLY: nnod

IMPLICIT NONE

REAL(DP), DIMENSION(:), INTENT(IN) :: newF, oldF
REAL(DP), INTENT(OUT) :: rel

REAL(DP) :: error
INTEGER :: n

rel = 0.

DO n= 1, nnod
    IF (ABS(newF(n)) > 1.e-10_DP) THEN
        error = ABS(newF(n) - oldF(n))
        IF (error > rel) rel = error
    END IF
END DO

END SUBROUTINE AbsE


SUBROUTINE par_ave(par, ave)
!
! Purpose:
!    To calculate average fuel temp (only for active core)
!

USE sdata, ONLY: vdel, nnod, ng, nuf

IMPLICIT NONE

REAL(DP), DIMENSION(:), INTENT(IN) :: par
REAL(DP), INTENT(OUT) :: ave
REAL(DP) :: dum, dum2
INTEGER :: n

dum = 0.; dum2 = 0.
DO n = 1, nnod
   IF (nuf(n,ng) > 0.) THEN
      dum = dum + par(n) * vdel(n)
      dum2 = dum2 + vdel(n)
   END IF
END DO

ave = dum / dum2

END SUBROUTINE par_ave


SUBROUTINE par_ave_out(par, ave)
!
! Purpose:
!    To calculate average fuel temp (only for active core)
!

USE sdata, ONLY: vdel, nnod, iz, nzz, nuf, ng, ix, iy, xyz, nxx, nyy, nzz, ystag

IMPLICIT NONE

REAL(DP), DIMENSION(:), INTENT(IN) :: par
REAL(DP), INTENT(OUT) :: ave
REAL(DP) :: dum, dum2
INTEGER, DIMENSION(nxx,nyy) :: zmax
INTEGER :: n, i, j, k

! get number of nodex in axial direction from bottom -> fuel
DO j = 1, nyy
  DO i = ystag(j)%smin, ystag(j)%smax
    zmax(i,j) = 0
    DO k = 1, nzz/2
      if (nuf(xyz(i,j,k),ng) < 1.e-5_dp) zmax(i,j) = zmax(i,j) + 1
    END DO
  END DO
END DO
! get number of nodex in axial direction from fuel -> top reflectors
DO j = 1, nyy
  DO i = ystag(j)%smin, ystag(j)%smax
    DO k = 1, nzz
      if (nuf(xyz(i,j,k),ng) > 1.e-5) zmax(i,j) = zmax(i,j) + 1
    END DO
  END DO
END DO

dum = 0.; dum2 = 0.
DO n = 1, nnod
   IF (iz(n) == zmax(ix(n),iy(n)) .AND. nuf(n,ng) > 1.e-5) THEN
      dum = dum + par(n) * vdel(n)
      dum2 = dum2 + vdel(n)
   END IF
END DO

ave = dum / dum2

END SUBROUTINE par_ave_out


SUBROUTINE par_max(par, pmax)
!
! Purpose:
!    To calculate maximum fuel tem, coolant tem, and density
!

USE sdata, ONLY: nnod

IMPLICIT NONE

REAL(DP), DIMENSION(:), INTENT(IN) :: par
REAL(DP), INTENT(OUT) :: pmax
INTEGER :: n

pmax = 0.
DO n = 1, nnod
   IF (par(n) > pmax) pmax = par(n)
END DO

END SUBROUTINE par_max


SUBROUTINE getent(t,ent)
!
! Purpose:
!    To get enthalpy for given coolant temp. from steam table
!

USE sdata, ONLY: stab, ntem
USE io, ONLY : ounit

IMPLICIT NONE

REAL(DP), INTENT(IN) :: t
REAL(DP), INTENT(OUT) :: ent
REAL(DP) :: t1, ent1
REAL(DP) :: t2, ent2
INTEGER :: i

IF ((t < stab(1,1)) .OR. (t > stab(ntem,1))) THEN
    WRITE(ounit,*) '  Coolant temp. : ', t, 'K'
    WRITE(ounit,*) '  ERROR : MODERATOR TEMP. IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
    WRITE(ounit,*) '  CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
    WRITE(*,*) '  Coolant temp. : ', t, 'K'
    WRITE(*,*) '  ERROR : MODERATOR TEMP. IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
    WRITE(*,*) '  CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
    STOP
END IF

t2 = stab(1,1); ent2 = stab(1,3)
DO i = 2, ntem
    t1 = t2
    ent1 = ent2
    t2 = stab(i,1); ent2 = stab(i,3)
    IF ((t >= t1) .AND. (t <= t2)) THEN
        ent = ent1 + (t - t1) / (t2 - t1) * (ent2 - ent1)
        EXIT
    END IF
END DO


END SUBROUTINE getent


SUBROUTINE gettd(ent,t,rho,prx,kvx,tcx,Rx)
!
! Purpose:
!    To get enthalpy for given coolant temp. from steam table
!

USE sdata, ONLY: stab, ntem
USE io, ONLY : ounit

IMPLICIT NONE

REAL(DP), INTENT(IN) :: ent
REAL(DP), INTENT(OUT) :: t, rho, prx, kvx, tcx
REAL(DP), INTENT(OUT), OPTIONAL :: Rx
REAL(DP) :: ratx

INTEGER :: i, i1, i2

! Get two closest interpolation points
IF (ent >= stab(1,3) .AND. ent <= stab(ntem,3)) THEN  !If enthalpy inside data range
  DO i = 2, ntem
    IF (ent >= stab(i-1,3) .AND. ent <= stab(i,3)) THEN
      i1 = i-1
      i2 = i
      EXIT
    END IF
  END DO
ELSE IF (ent < stab(1,3)  .AND. (stab(1,3) - ent) / stab(1,3) < 0.1) THEN !If 10% lower than min. steam table data
  i1 = 1
  i2 = 2
ELSE IF (ent > stab(ntem,3) .AND. (ent - stab(ntem,3)) / stab(ntem,3) < 0.1) THEN !if 10% higher than max. steam table data
  i1 = ntem-1
  i2 = ntem
ELSE
  WRITE(ounit,1557) ent/1000.
  WRITE(*,1557) ent/1000.
  WRITE(ounit,*) '   CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
  WRITE(*,*) '   CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
  STOP
  1557 FORMAT(2X, '  ERROR: ENTHALPY', F8.1 ,' KJ/Kg &
  & IS OUT OF THE RANGE IN THE STEAM TABLE')
END IF

! Interpolate
ratx = (ent - stab(i1,3)) / (stab(i2,3) - stab(i1,3))
t   = stab(i1,1) + ratx * (stab(i2,1) - stab(i1,1))
rho = stab(i1,2) + ratx * (stab(i2,2) - stab(i1,2))
prx = stab(i1,4) + ratx * (stab(i2,4) - stab(i1,4))
kvx = stab(i1,5) + ratx * (stab(i2,5) - stab(i1,5))
tcx = stab(i1,6) + ratx * (stab(i2,6) - stab(i1,6))
IF (PRESENT(Rx)) THEN
  Rx = 1000._DP * (stab(i2,2) - stab(i1,2)) / (stab(i2,3) - stab(i1,3))
END IF



END SUBROUTINE gettd


REAL(DP) FUNCTION getkc(t)
!
! Purpose:
!    To calculate thermal conductivity of cladding
!

IMPLICIT NONE

REAL(DP), INTENT(IN) :: t

getkc = 7.51_DP + 2.09e-2_DP*t - 1.45e-5_DP*t**2 + 7.67e-9_DP*t**3

END FUNCTION getkc


REAL(DP) FUNCTION getkf(t)
!
! Purpose:
!    To calculate thermal conductivity of fuel
!

IMPLICIT NONE

REAL(DP), INTENT(IN) :: t

getkf = 1.05_DP + 2150.0_DP / (t - 73.15_DP)

END FUNCTION getkf


REAL(DP) FUNCTION getcpc(t)
!
! Purpose:
!    To calculate specific heat capacity of cladding
!

IMPLICIT NONE

REAL(DP), INTENT(IN) :: t

getcpc = 252.54_DP + 0.11474_DP*t

END FUNCTION getcpc


REAL(DP) FUNCTION getcpf(t)
!
! Purpose:
!    To calculate specific heat capacity of fuel
!

IMPLICIT NONE

REAL(DP), INTENT(IN) :: t

getcpf = 162.3_DP + 0.3038_DP*t - 2.391e-4_DP*t**2 + 6.404e-8_DP*t**3

END FUNCTION getcpf


SUBROUTINE TridiaSolve(a,b,c,d,x)
!
! Purpose:
!    To solve tridiagonal matrix
!

IMPLICIT NONE

REAL(DP), DIMENSION(:), INTENT(INOUT) :: a, b, c, d
REAL(DP), DIMENSION(:), INTENT(OUT) :: x

INTEGER :: i, n

n = SIZE(d)

! Gauss Elimination
c(1) = c(1)/b(1)
d(1) = d(1)/b(1)
DO i = 2, n
    c(i) = c(i) / (b(i) - a(i) * c(i-1))
    d(i) = (d(i) - a(i) * d(i-1)) / (b(i) - a(i) * c(i-1))
END DO

! Back Substitution
x(n) = d(n)
DO i = n-1, 1, -1
    x(i) = d(i) - c(i) * x(i+1)
END DO

END SUBROUTINE TridiaSolve



REAL(DP) FUNCTION geths(xden, tc, kv, Pr)
!
! Purpose:
!    To calculate heat transfer coef.
!

USE sdata, ONLY: dh, farea, cflow

IMPLICIT NONE

REAL(DP), INTENT(IN) :: xden  ! coolant densisty
REAL(DP), INTENT(IN) :: tc  ! coolant thermal conductivity
REAL(DP), INTENT(IN) :: kv  ! kinematic viscosity
REAL(DP), INTENT(IN) :: Pr  ! Prandtl Number

REAL(DP) :: cvelo, Nu, Re

cvelo = cflow / (farea * xden * 1000._DP)        ! Calculate flow velocity (m/s)
Re = cvelo * dh / (kv * 1.e-6_DP)                 ! Calculate Reynolds Number
Nu = 0.023_DP*(Pr**0.4_DP)*(Re**0.8_DP)                ! Calculate Nusselt Number
geths = (tc / dh) * Nu                        ! Calculate heat transfer coefficient


END FUNCTION geths



SUBROUTINE th_trans(xpline, h)

!
! Purpose:
!    To perform fuel pin thermal transient
!

USE sdata, ONLY: nnod, mtem, cden, ftem, tin, cflow, nxx, nyy, cf, ent, heatf, &
                 tfm, nt, rpos, rdel, rf, rg, rc, farea, dia, pi, zdel, &
                 ix, iy, iz, frate, th_time, get_time

IMPLICIT NONE

REAL(DP), DIMENSION(:), INTENT(IN) :: xpline    ! Linear Power Density (W/cm)
REAL(DP), INTENT(IN) :: h                       ! Time step

INTEGER :: i, j, k, n
REAL(DP), DIMENSION(nt+1) :: a, b, c, d
REAL(DP) :: hs, hg = 1.d4, kt , kt1, kt2          ! coolant heat transfer coef., gap heat transfer coef, and thermal conductivity
REAL(DP) :: alpha = 0.7_DP
REAL(DP) :: xa, xc
REAL(DP) :: fdens = 10.412e3            ! UO2 density (kg/m3)
REAL(DP) :: cdens = 6.6e3               ! Cladding density (kg/m3)
REAL(DP) :: cp                          ! Specific heat capacity
REAL(DP) :: eps, eta
REAL(DP) :: mdens, vol                  ! Coolant density and channel volume
REAL(DP), DIMENSION(nnod) :: entp        ! previous enthalpy

REAL(DP) :: pdens      ! power densisty  (W/m3)
REAL(DP) :: enti       ! Coolant inlet enthalpy
REAL(DP), DIMENSION(nxx, nyy) :: entm, bfrate
REAL(DP) :: cpline     ! Coolant Linear power densisty (W/m)
REAL(DP) :: Pr, kv, tcon ! Coolant Prandtl Number, Kinematic viscosity, and thermal conductivity
real(dp) :: R
logical  :: first = .true.
real(dp) :: st, fn

!Get start th_time
st = get_time()

!set initial tridiagonal matrix element a, b
a = 0._dp; b = 0._dp; c = 0._dp;

if (first) then
  allocate(frate(nnod))
  frate = cflow
  first = .false.
end if
! frate = cflow
CALL getent(tin, enti)
entp = ent

DO n = 1, nnod
   mdens = cden(n) * 1000._DP                                    ! Coolant density (kg/m3)
   cpline = heatf(n) * pi * dia + cf * xpline(n) * 100._DP       ! Coolant Linear power densisty (W/m)
   vol   = farea * zdel(iz(n)) * 0.01_DP                             ! channel node volume
   i = ix(n); j = iy(n); k = iz(n)

   IF (k == 1) THEN                                              ! Calculate coolant enthalpy
       eps = mdens * vol / h
       ent(n) = (cpline * zdel(iz(n)) * 0.01_DP + 2._DP * frate(n) * enti &
              + eps * entp(n)) / (eps + 2._DP * frate(n))                             ! Calculate enthalpy
       CALL gettd(ent(n), mtem(n), cden(n), Pr, kv, tcon, R)        ! Get corresponding temp and density
       entm(i,j)   = 2._DP * ent(n) - enti
       frate(n)    = cflow - 0.5_dp * vol / h * R * (ent(n) - entp(n))
       bfrate(i,j) = 2._DP * frate(n) - cflow
   ELSE
       eps = mdens * vol / h
       ent(n) = (cpline * zdel(iz(n)) * 0.01_DP + 2._DP * frate(n) &
              * entm(ix(n),iy(n)) + eps * entp(n)) / (eps + 2._DP * frate(n))
       CALL gettd(ent(n), mtem(n), cden(n), Pr, kv, tcon, R)
       entm(i,j)   = 2._DP * ent(n) - entm(i,j)
       frate(n)    = bfrate(i,j) - 0.5_dp * vol / h * R * (ent(n) - entp(n))
       bfrate(i,j) = 2._DP * frate(n) - bfrate(i,j)
   END IF

   hs = geths(cden(n), Pr, kv, tcon)                                               ! Calculate heat transfer coef
   pdens = 100._DP * xpline(n) / (pi * rf**2)                ! Fuel pin Power Density (W/m3)

   ! Calculate tridiagonal matrix: a, b, c and source: d
   ! For nt=1 [FUEL CENTERLINE]
   kt1 = getkf(tfm(n,1))                                                     ! Get thermal conductivity
   kt2 = getkf(tfm(n,2))
   kt  = 2._DP * kt1 * kt2 / (kt1 + kt2)
   cp = getcpf(tfm(n,1))                                                           ! Get specific heat capacity
   eta = fdens * cp * rpos(1)**2 / (2._DP * h)
   xc  = kt * rpos(1) / rdel(1)
   b(1) =  xc + eta
   c(1) = -xc
   d(1) = pdens * 0.5_DP * rpos(1)**2 + eta * tfm(n,1)

   DO i = 2, nt-2
       kt1 = kt2
       kt2 = getkf(tfm(n,i+1))
       kt  = 2._DP * kt1 * kt2 / (kt1 + kt2)
       cp = getcpf(tfm(n,i))
       eta = fdens * cp * (rpos(i)**2 - rpos(i-1)**2) / (2. * h)
       xa = xc
       xc = kt * rpos(i) / rdel(i)
       a(i) = -xa
       b(i) =  xa + xc + eta
       c(i) = -xc
       d(i) = pdens * 0.5_DP * (rpos(i)**2 - rpos(i-1)**2) + eta * tfm(n,i)
   END DO

   ! For nt-1 [FUEL-GAP INTERFACE]
   cp = getcpf(tfm(n,nt-1))
   eta = fdens * cp * (rf**2 - rpos(nt-2)**2) / (2. * h)
   xa = xc
   xc = rg * hg
   a(nt-1) = -xa
   b(nt-1) =  xa + xc + eta
   c(nt-1) = -xc
   d(nt-1) = pdens * 0.5_DP * (rf**2 - rpos(nt-2)**2) + eta * tfm(n,nt-1)

   ! For nt [GAP-CLADDING INTERFACE]
   kt1 = getkc(tfm(n,nt))
   kt2 = getkc(tfm(n,nt+1))
   kt  = 2._DP * kt1 * kt2 / (kt1 + kt2)     ! For cladding
   cp = getcpc(tfm(n,nt))
   eta = cdens * cp * (rpos(nt)**2 - rg**2) / (2. * h)
   xa = xc
   xc = kt * rpos(nt) / rdel(nt)
   a(nt) = -xa
   b(nt) =  xa + xc + eta
   c(nt) = -xc
   d(nt) = eta * tfm(n,nt)

   ! For nt+1  [CLADDING-COOLANT INTERFACE]
   cp = getcpc(tfm(n,nt+1))
   eta = cdens * cp * (rc**2 - rpos(nt)**2) / (2. * h)
   xa = xc
   xc = rc * hs
   a(nt+1) = -xa
   b(nt+1) =  xa + xc + eta
   d(nt+1) = rc * hs * mtem(n) + eta * tfm(n,nt+1)

   ! Solve tridiagonal matrix
   CALL TridiaSolve(a, b, c, d, tfm(n, :))

   ! Get lumped fuel temp
   ftem(n) = (1.-alpha) * tfm(n, 1) + alpha * tfm(n, nt-1)

   ! Calculate heat flux
   heatf(n) = hs * (tfm(n, nt+1) - mtem(n))
END DO

!Get th_time
fn = get_time()
th_time = th_time + (fn-st)

END SUBROUTINE th_trans


SUBROUTINE th_upd(xpline)

!
! Purpose:
!    To update thermal parameters
!

USE sdata, ONLY: nnod, mtem, cden, ftem, tin, ix, iy, iz, nxx, nyy, cflow, cf,  &
                 ent, heatf, tfm, nt, rpos, rdel, rf, rg, rc, pi, zdel, dia

IMPLICIT NONE

REAL(DP), DIMENSION(:), INTENT(IN) :: xpline    ! Linear Power Density (W/cm)

INTEGER :: i, n
REAL(DP), DIMENSION(nt+1) :: a, b, c, d
REAL(DP) :: hs, Hg = 1.D4, kt, kt1, kt2
REAL(DP) :: alp = 0.7_DP
REAL(DP) :: xa, xc
REAL(DP) :: pdens      ! power densisty  (W/m3)
REAL(DP) :: enti       ! Coolant inlet enthalpy
REAL(DP), DIMENSION(nxx, nyy) :: entm   ! enthalpy at node boundary
REAL(DP) :: cpline     ! Coolant Linear power densisty (W/m)
REAL(DP) :: Pr, kv, tcon ! Coolant Prandtl Number, Kinematic viscosity, and thermal conductivity
REAL(DP) :: zd  ! zdel in meter

!set initial tridiagonal matrix element a, b
a = 0._dp; b = 0._dp; c = 0._dp;

CALL getent(tin, enti)

DO n = 1, nnod
  cpline = heatf(n) * pi * dia  + cf * xpline(n) * 100._DP          ! Coolant Linear power densisty (W/m)
  zd = zdel(iz(n)) * 0.01_DP
  IF (iz(n) == 1) THEN                                              ! For most bootom channel
    ent(n) = enti + 0.5_DP * cpline * zd / cflow    ! Calculate coolant enthalpy
    CALL gettd(ent(n), mtem(n), cden(n), Pr, kv, tcon)             ! Get corresponding temp and density
    entm(ix(n),iy(n)) = 2._DP * ent(n) - enti                      ! Extrapolate enthalpy at node boundary
  ELSE
    ent(n) = entm(ix(n),iy(n)) + 0.5_DP * cpline * zd / cflow
    CALL gettd(ent(n), mtem(n), cden(n), Pr, kv, tcon)
    entm(ix(n),iy(n)) = 2._DP * ent(n) - entm(ix(n),iy(n))
  END IF

  hs = geths(cden(n), Pr, kv, tcon)
  pdens = (1._DP - cf) * 100._DP * xpline(n) / (pi * rf**2)        ! Fuel pin Power Density (W/m3)

  ! Calculate tridiagonal matrix: a, b, c and source: d
  ! For nt=1 [FUEL CENTERLINE]
  kt1 = getkf(tfm(n,1))                                         ! Get thermal conductivity
  kt2 = getkf(tfm(n,2))
  kt  = 2._DP * kt1 * kt2 / (kt1 + kt2)
  xc  = kt * rpos(1) / rdel(1)
  b(1) =  xc
  c(1) = -xc
  d(1) = pdens * 0.5_DP * rpos(1)**2

  DO i = 2, nt-2
    kt1 = kt2
    kt2 = getkf(tfm(n,i+1))
    kt  = 2._DP * kt1 * kt2 / (kt1 + kt2)
    xa = xc
    xc = kt * rpos(i) / rdel(i)
    a(i) = -xa
    b(i) =  xa + xc
    c(i) = -xc
    d(i) = pdens * 0.5_DP * (rpos(i)**2 - rpos(i-1)**2)
  END DO

  ! For nt-1 [FUEL-GAP INTERFACE]
  xa = xc
  xc = rg * Hg
  a(nt-1) = -xa
  b(nt-1) =  xa + xc
  c(nt-1) = -xc
  d(nt-1) = pdens * 0.5_DP * (rf**2 - rpos(nt-2)**2)

  ! For nt [GAP-CLADDING INTERFACE]
  kt1 = getkc(tfm(n,nt))
  kt2 = getkc(tfm(n,nt+1))
  kt  = 2._DP * kt1 * kt2 / (kt1 + kt2)     ! For cladding
  xa = xc
  xc = kt * rpos(nt) / rdel(nt)
  a(nt) = -xa
  b(nt) =  xa + xc
  c(nt) = -xc
  d(nt) = 0.

  ! For nt+1  [CLADDING-COOLANT INTERFACE]
  xa = xc
  a(nt+1) = -xa
  b(nt+1) =  xa + hs * rc
  d(nt+1) = rc * hs * mtem(n)

  ! Solve tridiagonal matrix
  CALL TridiaSolve(a, b, c, d, tfm(n, :))

  ! Get lumped fuel temp
  ftem(n) = (1.-alp) * tfm(n, 1) + alp * tfm(n, nt-1)

  ! Calculate heat flux
  heatf(n) = hs * (tfm(n, nt+1) - mtem(n))
END DO


END SUBROUTINE th_upd


SUBROUTINE print_head()

!
! Purpose:
!    To print critical boron concentration header output
!

USE io, ONLY: ounit, scr, bther

IMPLICIT NONE

if (bther == 0) then
  ! File Output
  WRITE(ounit,*); WRITE(ounit,*)
  WRITE(ounit,2176); WRITE(ounit,2177); WRITE(ounit,2176)
  WRITE(ounit,*); WRITE(ounit,2178); WRITE(ounit,2179)
  if (scr) then
    ! Terminal Output
    WRITE(*,*); WRITE(*,*)
    WRITE(*,2176); WRITE(*,2177); WRITE(*,2176)
    WRITE(*,*); WRITE(*,2178); WRITE(*,2179)
  end if
else
  ! File Output
  WRITE(ounit,*); WRITE(ounit,*)
  WRITE(ounit,1176); WRITE(ounit,1177); WRITE(ounit,1176)
  WRITE(ounit,*); WRITE(ounit,1178); WRITE(ounit,1179)
  if (scr) then
    ! Terminal Output
    WRITE(*,*); WRITE(*,*)
    WRITE(*,1176); WRITE(*,1177); WRITE(*,1176)
    WRITE(*,*); WRITE(*,1178); WRITE(*,1179)
  end if
end if

2176 format(' ============================================================')
2177 format(12X,'CRITICAL BORON CONCENTRATION SEARCH')
2178 format('Itr  Boron Conc.   K-EFF     FLUX ERROR   FISS. SRC. ERROR')
2179 format(' -----------------------------------------------------------')
1176 format &
(' =========================================================================')
1177 format(19X,'CRITICAL BORON CONCENTRATION SEARCH')
1178 format &
('Itr  Boron Conc.   K-EFF     FLUX ERR.    FISS. SRC. ERR.  DOPPLER ERR.')
1179 format &
(' -----------------------------------------------------------------------')

END SUBROUTINE print_head


SUBROUTINE cbsearch()

!
! Purpose:
!    To search critical boron concentration
!

USE sdata, ONLY: Ke, rbcon, ftem, mtem, cden, bpos, nnod, f0, fer, ser, &
                 aprad, apaxi, afrad, npow
USE io, ONLY: ounit, AsmFlux, AsmPow, AxiPow
USE cmfd, ONLY: outer, powdis
USE xsec, ONLY: XS_updt

IMPLICIT NONE

REAL(DP)  :: bc1, bc2, bcon     ! Boron Concentration
REAL(DP) :: ke1, ke2
INTEGER :: n

call print_head()

bcon = rbcon
CALL XS_updt(bcon, ftem, mtem, cden, bpos)
CALL outer(0)
bc1 = bcon
ke1 = Ke

WRITE(ounit,1791) 1, bc1, Ke1, ser, fer
WRITE(*,1791) 1, bc1, Ke1, ser, fer

bcon = bcon + (Ke - 1.) * bcon   ! Guess next critical boron concentration
CALL XS_updt(bcon, ftem, mtem, cden, bpos)
CALL outer(0)
bc2 = bcon
ke2 = Ke

WRITE(ounit,1791) 2, bc2, Ke2, ser, fer
WRITE(*,1791) 2, bc2, Ke2, ser, fer

n = 3
DO
  bcon = bc2 + (1._DP - ke2) / (ke1 - ke2) * (bc1 - bc2)
  CALL XS_updt(bcon, ftem, mtem, cden, bpos)
  CALL outer(0)
  bc1 = bc2
  bc2 = bcon
  ke1 = ke2
  ke2 = ke
  WRITE(ounit,1791) n, bcon, Ke, ser, fer
  WRITE(*,1791) n, bcon, Ke, ser, fer
    IF ((ABS(Ke - 1._DP) < 1.e-5_DP) .AND. (ser < 1.e-5_DP) .AND. (fer < 1.e-5_DP)) EXIT
    n = n + 1
    IF (bcon > 3000.) THEN
        WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
        WRITE(ounit,*) '  ADPRES IS STOPPING'
        WRITE(*,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
        STOP
    END IF
    IF (bcon < 0.) THEN
        WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
        WRITE(ounit,*) '  ADPRES IS STOPPING'
        WRITE(*,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
        STOP
    END IF
    IF (n == 20) THEN
        WRITE(ounit,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
        WRITE(ounit,*) '  ADPRES IS STOPPING'
        WRITE(*,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
        STOP
    END IF
END DO

ALLOCATE(npow(nnod))
IF (aprad == 1 .OR. apaxi == 1) THEN
    CALL PowDis(npow)
END IF

IF (aprad == 1) CALL AsmPow(npow)

IF (apaxi == 1) CALL AxiPow(npow)

IF (afrad == 1) CALL AsmFlux(f0, 1._DP)

1791 format(I3, F10.2, F14.5, ES14.5, ES13.5)

END SUBROUTINE cbsearch


SUBROUTINE cbsearcht()

!
! Purpose:
!    To search critical boron concentration with thermal feedback
!

USE sdata, ONLY: Ke, ftem, mtem, cden, bcon, rbcon, npow, nnod, &
                 f0, ser, fer, tfm, aprad, apaxi, afrad, npow, th_err, &
                 serc, ferc
USE io, ONLY: ounit, AsmFlux, AsmPow, AxiPow, scr
USE cmfd, ONLY: powdis, outer

IMPLICIT NONE

REAL(DP)  :: bc1, bc2    ! Boron Concentration
REAL(DP) :: ke1, ke2
INTEGER :: n
REAL(DP) :: tf, tm, mtm, mtf, otm, cd, ocd

call print_head()

ALLOCATE(npow(nnod))

bcon = rbcon
CALL th_iter()  ! Start thermal hydarulic iteration with current paramters
bc1 = bcon
ke1 = Ke

WRITE(ounit,1792) 1, bc1, Ke1, ser, fer, th_err
WRITE(*,1792) 1, bc1, Ke1, ser, fer, th_err

IF (bcon < 1.e-5) THEN
  bcon = 500.
ELSE
  bcon = bcon + (Ke - 1.) * bcon   ! Guess next critical boron concentration
END IF
CALL th_iter()                 ! Perform second thermal hydarulic iteration with updated parameters
bc2 = bcon
ke2 = Ke

WRITE(ounit,1792) 2, bc2, Ke2, ser, fer, th_err
WRITE(*,1792) 2, bc2, Ke2, ser, fer, th_err

n = 3
DO
    bcon = bc2 + (1._DP - ke2) / (ke1 - ke2) * (bc1 - bc2)
    CALL th_iter()
    bc1 = bc2
    bc2 = bcon
    ke1 = ke2
    ke2 = ke
    WRITE(ounit,1792) n, bcon, Ke, ser, fer, th_err
    WRITE(*,1792) n, bcon, Ke, ser, fer, th_err
    IF ((ABS(Ke - 1._DP) < 1.e-5_DP) .AND. (ser < serc) .AND. (fer < ferc)) EXIT
    n = n + 1
    IF (bcon > 3000.) THEN
        WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
        WRITE(ounit,*) '  ADPRES IS STOPPING'
        WRITE(*,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
        STOP
    END IF
    IF (bcon < 0.) THEN
        WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
        WRITE(ounit,*) '  ADPRES IS STOPPING'
        WRITE(*,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
        STOP
    END IF
    IF (n == 30) THEN
        WRITE(ounit,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
        WRITE(ounit,*) '  ADPRES IS STOPPING'
        WRITE(*,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
        STOP
    END IF
END DO

IF (aprad == 1 .OR. apaxi == 1) THEN
    CALL PowDis(npow)
END IF

IF (aprad == 1) CALL AsmPow(npow)

IF (apaxi == 1) CALL AxiPow(npow)

IF (afrad == 1) CALL AsmFlux(f0, 1._DP)

CALL par_ave(ftem, tf)
CALL par_ave(mtem, tm)

CALL par_max(tfm(:,1), mtf)
CALL par_max(mtem, mtm)

CALL par_ave_out(mtem, otm)
CALL par_ave(cden, cd)
CALL par_ave_out(cden, ocd)

! Write Output
WRITE(ounit,*)
WRITE(ounit, 5001) tf, tf-273.15; WRITE(ounit, 5002)  mtf, mtf-273.15
WRITE(ounit, 5003) tm, tm-273.15; WRITE(ounit, 5004) mtm, mtm-273.15
WRITE(ounit, 5005) otm, otm-273.15; WRITE(ounit, 5006) cd * 1000., cd
WRITE(ounit, 5007) ocd * 1000., ocd
if (scr) then
  WRITE(*,*)
  WRITE(*, 5001) tf, tf-273.15; WRITE(*, 5002)  mtf, mtf-273.15
  WRITE(*, 5003) tm, tm-273.15; WRITE(*, 5004) mtm, mtm-273.15
  WRITE(*, 5005) otm, otm-273.15; WRITE(*, 5006) cd * 1000., cd
  WRITE(*, 5007) ocd * 1000., ocd
end if

5001 FORMAT(2X, 'AVERAGE DOPPLER TEMPERATURE     : ', F7.1, ' K (', F7.1, ' C)')
5002 FORMAT(2X, 'MAX FUEL CENTERLINE TEMPERATURE : ', F7.1, ' K (', F7.1, ' C)')
5003 FORMAT(2X, 'AVERAGE MODERATOR TEMPERATURE   : ', F7.1, ' K (', F7.1, ' C)')
5004 FORMAT(2X, 'MAXIMUM MODERATOR TEMPERATURE   : ', F7.1, ' K (', F7.1, ' C)')
5005 FORMAT(2X, 'OUTLET MODERATOR TEMPERATURE    : ', F7.1, ' K (', F7.1, ' C)')
5006 FORMAT(2X, 'AVERAGE MODERATOR DENSITY       : ', F7.1, ' kg/m3 (', F7.3, ' g/cc)')
5007 FORMAT(2X, 'OUTLET MODERATOR DENSITY        : ', F7.1, ' kg/m3 (', F7.3, ' g/cc)')
1792 format(I3, F9.2, F14.5, ES14.5, ES13.5, ES17.5)

END SUBROUTINE cbsearcht


END MODULE th
