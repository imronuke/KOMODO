MODULE trans

!=========================
! Transient Module to solve transient diffusion problems
! Using Fully Implicit method with (or without) exponetial transformation
! =======================

USE data, ONLY: DP

implicit none

SAVE

CONTAINS


subroutine rod_eject()

!
! Purpose:
!    To perform rod ejection simulation
!

USE data, ONLY: ng, nnod, sigr, nf, nmat, total_time, time_mid, time_step_1, time_step_2, Ke, &
                 bcon, ftem, mtem, cden, bpos, nbank, &
                 beta, flux, flux_prev, c0, tbeta, omega, core_beta, L
USE read, ONLY: ounit, bextr, scr
USE xsec, ONLY: XS_updt
USE cmfd, ONLY: outer_tr, outer, outer_ad

implicit none

real(dp), dimension(nnod, ng) :: af                                      ! adjoint flux

real(dp) :: rho
real(dp) :: t1, t2
real(dp) :: tpow1
integer :: n, i, j, imax, step

! Allocate precusor density
allocate (c0(nnod,nf))

! Allocate Frequency transformation constant
allocate (omega(nnod,ng))

! allocate total leakages
allocate(L(nnod,ng))

! Update xsec
CALL XS_updt(bcon, ftem, mtem, cden, bpos)

! Calculate forward flux at t=0 and check if keff=1
CALL outer(0)
if (scr) then
  write(*,*)
  write(*,*) ' steady state calculation ... done'
end if

! If K-EFF NOT EQUAL TO 1.0
if (ABS(Ke - 1._DP) > 1.e-5_DP) CALL KNE1

! Calculate Adjoint flux
! NOTE: This adjoint flux is approximation where
! the same nodal coupling coefficients in forward calculation are used
CALL outer_ad(0)
af = flux   ! Save adjoint flux to af
if (scr) then
  write(*,*)
  write(*,*) ' adjoint calculation ... done'
end if

! ! ReCalculate forward flux
CALL outer(0)
if (scr) then
  write(*,*)
  write(*,*) ' re-calculate steady state condition ... done'
end if

! Calculate Initial precursor density
CALL iPden()

! Calculate initial power
CALL PowTot(flux,tpow1)

! Total beta
tbeta = 0.
do n = 1, nmat
  do j = 1, nf
    tbeta(n) = tbeta(n) + beta(j)
  end do
  core_beta = tbeta(1)
end do

! Calculate reactivity
call reactivity(af, sigr, rho)

! File output
WRITE(*, *)
WRITE(*, *) " TRANSIENT RESULTS :"
WRITE(*, *)
WRITE(*, *) " Step  Time(s)  React.($)   Rel. Power   CR Bank Pos. (1-end)"
WRITE(*, *) "--------------------------------------------------------------"
WRITE(*,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') 0, 0., rho, &
1.0, (bpos(n), n = 1, nbank)

! Terminal output
WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit, *) " TRANSIENT RESULTS :"
WRITE(ounit, *)
WRITE(ounit, *) " Step  Time(s)  React.($)   Rel. Power"
WRITE(ounit, *) "-------------------------------------------"
WRITE(ounit,'(I4, F10.3, F10.4, ES15.4)') 0, 0., rho, 1.0

! Start transient calculation
step = 0
t2 = 0.
imax = NINT(time_mid/time_step_1)

! First Time Step
do i = 1, imax

    step = step + 1
    t1 = t2
    t2 = REAL(i)*time_step_1

    if (i > 1) THEN
       if (bextr == 0) then
         omega = 0._dp
       else
         omega = LOG(flux / flux_prev) / time_step_1
       end if
    else
       omega = 0._dp
    end if

    call trans_calc(0,time_step_1,af,tpow1,step,t2)

end do

! Second Time Step
imax = NINT((total_time-time_mid)/time_step_2)

do i = 1, imax

    step = step + 1
    t1 = t2
    t2 = time_mid + REAL(i)*time_step_2

    if (bextr == 0) then
      omega = 0._dp
    else
      omega = LOG(flux / flux_prev) / time_step_2
    end if

    call trans_calc(0,time_step_2,af,tpow1,step,t2)

end do

end subroutine rod_eject

!****************************************************************************!

subroutine rod_eject_th()

!
! Purpose:
!    To perform rod ejection simulation with TH feedback
!

USE data, ONLY: ng, nnod, sigr, nf, total_time, time_mid, time_step_1, time_step_2, Ke, &
                 tfm, ppow, m, ftem, mtem, bpos, nbank, beta, &
                 flux, flux_prev, c0, tbeta, omega, node_power, core_beta, nmat, L, n_th_iter
USE read, ONLY: ounit, bextr, bxtab, scr
USE xsec, ONLY: XS_updt
USE cmfd, ONLY: outer_tr, outer, outer_ad
use th, only : th_iter, par_ave, par_max

implicit none

real(dp), dimension(nnod, ng) :: af                                      ! adjoint flux

real(dp) :: rho
real(dp) :: t1, t2
real(dp) :: tpow1
integer :: n, i, j, imax, step
real(dp) :: tf, tm, mtf, mtm

! Allocate precusor density
allocate (c0(nnod,nf))

! Allocate Frequency transformation constant
allocate (omega(nnod,ng))

! Allocate node power distribution
allocate(node_power(nnod))

! allocate total leakages
allocate(L(nnod,ng))

! Determine th paramters distribution
CALL th_iter(n_th_iter, 0)
if (scr) then
  write(*,*)
  write(*,*) ' steady state calculation ... done'
end if

! If K-EFF NOT EQUAL TO 1.0
IF (ABS(Ke - 1._DP) > 1.e-5_DP .AND. bxtab == 0) CALL KNE1

! Calculate Adjoint flux
! NOTE: This adjoint flux is approximation where
! the same nodal coupling coefficients in forward calculation are used
CALL outer_ad(0)
af = flux   ! Save adjoint flux to af
if (scr) then
  write(*,*)
  write(*,*) ' adjoint calculation ... done'
end if

! Calculate forward flux
CALL outer(0)
if (scr) then
  write(*,*)
  write(*,*) ' re-calculate steady state condition ... done'
end if

! Calculate power
CALL PowTot(flux,tpow1)

! Calculate Initial precursor density
CALL iPden()

! Total beta
IF (bxtab == 0) then
  tbeta = 0.
  do n = 1, nmat
    do j = 1, nf
      tbeta(n) = tbeta(n) + beta(j)
    end do
  end do
  core_beta = tbeta(1)
else
  tbeta = 0.
  do n = 1, nmat
    do j = 1, nf
      tbeta(n) = tbeta(n) + m(n)%iBeta(j)
    end do
  end do
  call calc_beta(af)
  write(*,*)
  write(*,1324) core_beta*1.e5
end if

! Calculate reactivity
call reactivity(af, sigr, rho)

CALL par_ave(ftem, tf)
CALL par_max(tfm(:,1), mtf)
CALL par_ave(mtem, tm)
CALL par_max(mtem, mtm)

! File output
WRITE(ounit, *)
WRITE(ounit, *) " TRANSIENT RESULTS :"
WRITE(ounit, *)
WRITE(ounit, *) " Step  Time(s)  React.($)   Rel. Power   Avg. Tm   Max. Tm   Avg. Tf   Max. Tf"
WRITE(ounit, *) "--------------------------------------------------------------------------------"
WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 4F10.2)') 0, 0., rho, &
ppow*0.01, tm-273.15, mtm-273.15, tf-273.15, mtf-273.15

! Terminal output
WRITE(*, *)
WRITE(*, *) " TRANSIENT RESULTS :"
WRITE(*, *)
WRITE(*, *) " Step  Time(s)  React.($)   Rel. Power   CR Bank Pos. (1-end)"
WRITE(*, *) "--------------------------------------------------------------"
WRITE(*,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') 0, 0., rho, &
ppow*0.01, (bpos(n), n = 1, nbank)

! Start transient calculation
step = 0
t2 = 0.
imax = NINT(time_mid/time_step_1)

! First Time Step
do i = 1, imax

    step = step + 1
    t1 = t2
    t2 = REAL(i)*time_step_1

    if (i > 1) THEN
       if (bextr == 0) then
         omega = 0._dp
       else
         omega = LOG(flux / flux_prev) / time_step_1
       end if
    else
       omega = 0._dp
    end if

    call trans_calc(1,time_step_1,af,tpow1,step,t2)

end do

! Second Time Step
imax = NINT((total_time-time_mid)/time_step_2)

do i = 1, imax

    step = step + 1
    t1 = t2
    t2 = time_mid + REAL(i)*time_step_2

    if (bextr == 0) then
      omega = 0._dp
    else
      omega = LOG(flux / flux_prev) / time_step_2
    end if

    call trans_calc(1,time_step_2,af,tpow1,step,t2)

end do

1324 format(2X,'Core-averaged delayed neutron fraction :', F7.2, ' pcm')

end subroutine rod_eject_th

!****************************************************************************!

subroutine trans_calc(thc,ht,af,tpow1,step,t2)

!
! Purpose:
!    To perform transient calculation for given time step
!

USE data, ONLY: ng, nnod, sigr, bcon, ftem, mtem, cden, &
                 fbpos, bpos, tmove, bspeed, mdir, nbank, neutron_velo, node_power, &
                 flux, flux_prev, fsrc_prev, fsrc, omega, transient_warning, ix, iy, iz, pow, &
                 tfm, zdel, ppow, node_nf, m, mat, dfis, core_beta, small_theta, sigrp
USE read, ONLY: ounit, bxtab
USE xsec, ONLY: XS_updt
USE cmfd, ONLY: outer_tr, outer
use th, only : th_trans, par_ave, par_max, powdis


implicit none

integer, intent(in)                  :: thc            !T-H indicator
real(dp), intent(in)                 :: ht, tpow1, t2
integer, intent(in)                  :: step
real(dp), dimension(:,:), intent(in) :: af             ! adjoint flux

real(dp) :: rho
real(dp) :: tpow2
integer :: n, g
logical :: maxi   ! Maximum Outer Iteration Reached?

real(dp), dimension(nnod) :: pline       ! Linear power density
real(dp) :: xppow
real(dp) :: tf, tm, mtf, mtm
logical  :: first = .true.

if (first) then
  allocate(flux_prev(nnod,ng), fsrc_prev(nnod))
  allocate(dfis(nnod))
  allocate(sigrp(nnod,ng))
  first = .false.
end if

! Rod bank changes
do n = 1, nbank
    if (mdir(n) == 1) THEN   ! If CR moving down
        if (t2-tmove(n) > 1.e-5_DP .AND. fbpos(n)-bpos(n) < 1.e-5_DP) THEN
            bpos(n) = bpos(n) - ht *  bspeed(n)
            if (bpos(n) < fbpos(n)) bpos(n) = fbpos(n)  ! If bpos exceed, set to fbpos
        end if
    else if (mdir(n) == 2) THEN ! If CR moving up
        if (t2-tmove(n) > 1.e-5_DP .AND. fbpos(n)-bpos(n) > 1.e-5_DP) THEN
            bpos(n) = bpos(n) + ht *  bspeed(n)
            if (bpos(n) > fbpos(n)) bpos(n) = fbpos(n)  ! If bpos exceed, set to fbpos
        end if
    else
        CONTINUE
    end if
 end do

 ! Calculate xsec after pertubation
 CALL XS_updt(bcon, ftem, mtem, cden, bpos)

 ! Modify removal xsec
  sigrp = sigr    ! Save sigr to sigrp
  IF (bxtab == 0) THEN
    DO g = 1, ng
      DO n = 1, nnod
         sigr(n,g) = sigr(n,g) + 1._DP / (small_theta * neutron_velo(g) * ht) + omega(n,g) / neutron_velo(g)
      END DO
    END DO
  ELSE
    DO g = 1, ng
      DO n = 1, nnod
        sigr(n,g) = sigr(n,g) + 1._DP / (small_theta * m(mat(n))%velo(g) * ht) &
                  + omega(n,g) / m(mat(n))%velo(g)
      END DO
    END DO
  END IF

! Save the previous fluxes and fission source
flux_prev = flux
fsrc_prev = fsrc

! Transient calculation
CALL outer_tr(ht, maxi)

! Update precursor density
CALL uPden(ht)

! Calculate power
CALL PowTot(flux, tpow2)

! Calculate reactivity
call reactivity(af, sigrp, rho)

if (thc == 1) then
  ! Calculate node power distribution
  CALL PowDis(node_power)

  ! Power change
  xppow = ppow * tpow2/tpow1 * 0.01_DP

  ! Calculate linear power density for each nodes (W/cm)
  DO n = 1, nnod
     pline(n) = node_power(n) * pow * xppow &
     / (node_nf(ix(n),iy(n)) * zdel(iz(n)))     ! Linear power density (W/cm)
  END DO

  ! TH transient
  CALL th_trans(pline,ht)

  CALL par_ave(ftem, tf)
  CALL par_max(tfm(:,1), mtf)
  CALL par_ave(mtem, tm)
  CALL par_max(mtem, mtm)
end if

if (thc == 1) then
  WRITE(*,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, rho/core_beta, &
  xppow, (bpos(n), n = 1, nbank)

  IF (maxi) THEN
    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 4F10.2, A35)') step, t2, rho/core_beta, &
    xppow, tm-273.15, mtm-273.15, tf-273.15, mtf-273.15, 'OUTER ITERATION DID NOT CONVERGE'
  ELSE
    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 4F10.2)') step, t2, rho/core_beta, &
    xppow, tm-273.15, mtm-273.15, tf-273.15, mtf-273.15
  END IF
else
  WRITE(*,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, rho/core_beta, &
  tpow2/tpow1, (bpos(n), n = 1, nbank)

  if (maxi) THEN
      WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, A35)') step, t2, rho/core_beta, &
      tpow2/tpow1, 'OUTER ITERATION DID NOT CONVERGE'
  else
      WRITE(ounit,'(I4, F10.3, F10.4, ES15.4)') step, t2, rho/core_beta, &
      tpow2/tpow1
  end if
end if

if (maxi) transient_warning = .TRUE.


end subroutine trans_calc

!****************************************************************************!

subroutine KNE1()

!
! Purpose:
!    To adjuts the Keff to 1.0 if it is not equal to 1.0
!

USE data, ONLY: Ke, xnuf, dnuf, bcon, ftem, mtem, cden, bpos
USE read, ONLY: ounit, scr
USE xsec, ONLY: XS_updt
USE cmfd, ONLY: outer

implicit none

integer :: i

WRITE(ounit, *)
WRITE(ounit, '(A46,F9.6)') '  INITIAL MULTIPLICATION EFFECTIVE (K-EFF) = ', Ke
WRITE(ounit, *) '  WARNING: THE STEADY STATE K-EFF IS NOT EQUAL TO 1.0'
WRITE(ounit, *) '  AND NOW IT IS FORCED TO 1.0 BY MODIFYING THE nu*sigf CROSS SECTIONS '
WRITE(ounit, *)
if (scr) then
  write(*,*)
  write(*,*) ' steady state k-eff not equal to one, force it to one ... done'
end if
do i = 1, 10
   xnuf = xnuf / Ke
   dnuf = dnuf / Ke
   CALL XS_updt(bcon, ftem, mtem, cden, bpos)
   CALL outer(0)
   if (ABS(Ke-1._DP) < 1.e-5_DP) EXIT
end do
if (i == 10) STOP "K-EFF STILL NOT EQUAL TO ONE. KOMODO IS STOPPING"


end subroutine KNE1


!******************************************************************************!

subroutine PowTot (fx,tpow)

!
! Purpose:
!    To calculate total power density
!


USE data, ONLY: ng, nnod, sigf, vdel

implicit none

real(dp), dimension(:,:), INTENT(IN) :: fx
real(dp), INTENT(OUT) :: tpow

real(dp), dimension(nnod) :: p
integer :: g, n
real(dp) :: pow

p = 0._DP
do g= 1, ng
    do n= 1, nnod
      pow = fx(n,g) * sigf(n,g) * vdel(n)
      if (pow < 0.) pow = 0.
      p(n) = p(n) + pow
    end do
end do


tpow = 0._DP
do n = 1, nnod
    tpow = tpow + p(n)
end do

end subroutine PowTot

!****************************************************************************!

subroutine iPden()

!
! Purpose:
!    Calculate Initial precursor density
!

USE data, ONLY: nnod, nf, fsrc, c0, beta, lambda, m, mat, nuf, ng
USE read, ONLY: bxtab

implicit none

integer :: n, j
real(dp) :: blamb

if (bxtab == 1) THEN
  do n = 1, nnod
     do j = 1, nf
       if (nuf(n,ng) > 0.) THEN  !If it is fuel
         blamb = m(mat(n))%iBeta(j) / m(mat(n))%lamb(j)
         c0(n,j)  = blamb * fsrc(n)
       else
         c0(n,j) = 0.
       end if
     end do
  end do
else
  do n = 1, nnod
     do j = 1, nf
        blamb = beta(j) / lambda(j)
        c0(n,j)  = blamb * fsrc(n)
     end do
  end do
end if


end subroutine iPden

!****************************************************************************!

subroutine uPden(ht)

!
! Purpose:
!    To update precursor density
!

USE data, ONLY: nnod, nf, fsrc, fsrc_prev, c0, beta, lambda,  m, mat, nuf, ng
USE read, ONLY: bxtab

implicit none

real(dp), INTENT(IN) :: ht
real(dp) :: a1, a2, pxe
integer  :: n, i

if (bxtab == 1) THEN
  do i = 1, nf
    do n = 1, nnod
      if (nuf(n,ng) > 0.) THEN  !If it is fuel
        pxe  = exp(-m(mat(n))%lamb(i)*ht)
        a1   = (1._dp - pxe) / (m(mat(n))%lamb(i)*ht)
        a2   = 1._dp - a1
        a1   = a1 - pxe
        c0(n,i)  = c0(n,i)  * pxe + m(mat(n))%iBeta(i) / m(mat(n))%lamb(i) &
                 * (a1*fsrc_prev(n) + a2*fsrc(n))
      end if
    end do
  end do
else
  do i = 1, nf
    pxe  = exp(-lambda(i)*ht)
    a1   = (1._dp - pxe) / (lambda(i)*ht)
    a2   = 1._dp - a1
    a1   = a1 - pxe
    do n = 1, nnod
      c0(n,i)  = c0(n,i)  * pxe + beta(i) / lambda(i) &
               * (a1*fsrc_prev(n) + a2*fsrc(n))
    end do
  end do
end if


end subroutine uPden

!****************************************************************************!

subroutine reactivity(af,sigrp, rho)

!
! Purpose:
!    To calculate dynamic reactivity
!

USE data, ONLY: nnod, ng, flux, sigs, chi, mat, fsrc, vdel, L
use nodal, only: Lxyz

implicit none

real(dp), dimension(:,:), INTENT(IN) :: af
real(dp), dimension(:,:), INTENT(IN) :: sigrp
real(dp), INTENT(OUT) :: rho

integer :: n, g, h
real(dp), dimension(nnod) :: scg
real(dp) :: rem, lea, src, fde, L1, L2, L3

src = 0.; rem = 0.; lea = 0.; fde = 0.
do g = 1, ng
  scg = 0.
  do h = 1, ng
     do n = 1, nnod
        if (g /= h) scg(n) = scg(n) + sigs(n,h,g) * flux(n,h)
     end do
  end do
  do n = 1, nnod
    call Lxyz(n,g,L1,L2,L3)
    L(n,g) = L1 + L2 + L3
    src = src + af(n,g) * (scg(n) + chi(mat(n),g) * fsrc(n)) * vdel(n)
    rem = rem + af(n,g) * sigrp(n,g) * flux(n,g) * vdel(n)
    lea = lea + af(n,g) * L(n,g) * vdel(n)
    fde = fde + af(n,g) * chi(mat(n),g) * fsrc(n) * vdel(n)
   end do
end do

rho = (src - lea - rem) / fde

end subroutine reactivity

!****************************************************************************!

subroutine calc_beta(af)

!
! Purpose:
!    To calculate core-averaged delayed neutron fraction
!

USE data, ONLY: ng, nnod, m, mat, chi, nf, flux, nuf, core_beta
USE cmfd, ONLY: Integrate
USE read, ONLY: ounit

implicit none

real(dp), dimension(:,:), INTENT(IN) :: af

INTEGER :: n, i, g
real(dp), dimension(nnod) :: vdum, vdum2
real(dp) :: F

! Calculate F
vdum = 0.
DO g = 1, ng
    DO n = 1, nnod
        vdum(n) = vdum(n) + nuf(n,g) * flux(n,g)
    END DO
END DO

vdum2 = 0.
DO g = 1, ng
    DO n = 1, nnod
        vdum2(n) = vdum2(n) + chi(mat(n),g) * vdum(n) * af(n,g)
    END DO
END DO

F = Integrate(vdum2)

! Calculate Delayed neutron fraction (beta)
core_beta = 0._dp
DO i = 1, nf
    vdum2 = 0.
    DO g = 1, ng
        DO n = 1, nnod
            vdum2(n) = vdum2(n) + chi(mat(n),g) * m(mat(n))%iBeta(i) * vdum(n) * af(n,g)
        END DO
    END DO
    core_beta = core_beta + Integrate(vdum2) / F
END DO

write(ounit,*)
write(ounit,1344) core_beta*1.e5

1344 format ('  CORE AVERAGED DELAYED NEUTRON FRACTION: ', F7.2, ' PCM')


end subroutine calc_beta


end MODULE trans
