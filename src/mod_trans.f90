module trans

   !=========================
   ! Transient Module to solve transient diffusion problems
   ! Using Fully Implicit method with (or without) exponetial transformation
   ! =======================

   use sdata, only: dp

   implicit none

   save

   contains


   subroutine rod_eject()

      !
      ! Purpose:
      !    To perform rod ejection simulation
      !

      use sdata, only: ng, nnod, sigr, nf, nmat, ttot, tdiv, tstep1, tstep2, ke, &
      bcon, ftem, mtem, cden, bpos, nb, &
      ibeta, f0, ft, c0, tbeta, omeg, ctbeta, l, npow
      use io, only: ounit, bextr, scr, bvtk, print_vtk
      use hdf5_output, only: hdf5_is_active, hdf5_write_step
      use xsec, only: xs_updt
      use cmfd, only: outer_tr, outer, outer_ad
      use th, only : powdis

      implicit none

      real(dp), dimension(nnod, ng) :: af                                      ! adjoint flux

      real(dp) :: rho
      real(dp) :: t1, t2
      real(dp) :: tpow1
      integer :: n, i, j, imax, step

      ! Allocate precusor density
      allocate (c0(nnod, nf))

      ! Allocate Frequency transformation constant
      allocate (omeg(nnod, ng))

      ! allocate total leakages
      allocate(l(nnod, ng))

      ! Update xsec
      call xs_updt(bcon, ftem, mtem, cden, bpos)

      ! Calculate forward flux at t=0 and check if keff=1
      call outer(0)
      if (scr) then
         write(*, *)
         write(*, *) ' steady state calculation ... done'
      end if

      ! If K-EFF NOT EQUAL TO 1.0
      if (abs(ke - 1._dp) > 1.e-5_dp) call kne1

      ! Calculate Adjoint flux
      ! NOTE: This adjoint flux is approximation where
      ! the same nodal coupling coefficients in forward calculation are used
      call outer_ad(0)
      af = f0   ! Save adjoint flux to af
      if (scr) then
         write(*, *)
         write(*, *) ' adjoint calculation ... done'
      end if

      ! ! ReCalculate forward flux
      call outer(0)
      if (scr) then
         write(*, *)
         write(*, *) ' re-calculate steady state condition ... done'
      end if

      ! Calculate Initial precursor density
      call ipden()

      ! Calculate initial power
      call powtot(f0, tpow1)

      ! Total beta
      tbeta = 0.
      do n = 1, nmat
         do j = 1, nf
            tbeta(n) = tbeta(n) + ibeta(j)
         end do
         ctbeta = tbeta(1)
      end do

      ! Calculate reactivity
      call reactivity(af, sigr, rho)

      ! print vtk
      allocate(npow(nnod))
      call powdis(npow)
      if (hdf5_is_active) call hdf5_write_step(0, 0._dp, reactivity = rho)
      if (bvtk == 1) call print_vtk(0)

      ! File output
      write(*, *)
      write(*, *) " TRANSIENT RESULTS :"
      write(*, *)
      write(*, *) " Step  Time(s)  React.($)   Rel. Power   CR Bank Pos. (1-end)"
      write(*, *) "--------------------------------------------------------------"
      write(*, '(I4, F10.3, F10.4, ES15.4, 12F9.2)') 0, 0. , rho, &
      1.0, (bpos(n), n = 1, nb)

      ! Terminal output
      write(ounit, *)
      write(ounit, *)
      write(ounit, *) " TRANSIENT RESULTS :"
      write(ounit, *)
      write(ounit, *) " Step  Time(s)  React.($)   Rel. Power"
      write(ounit, *) "-------------------------------------------"
      write(ounit, '(I4, F10.3, F10.4, ES15.4)') 0, 0. , rho, 1.0

      ! Start transient calculation
      step = 0
      t2 = 0.
      imax = nint(tdiv / tstep1)

      ! First Time Step
      do i = 1, imax

         step = step + 1
         t1 = t2
         t2 = real(i) * tstep1

         if (i > 1) then
            if (bextr == 0) then
               omeg = 0._dp
            else
               omeg = log(f0 / ft) / tstep1
            end if
         else
            omeg = 0._dp
         end if

         call trans_calc(0, tstep1, af, tpow1, step, t2)

      end do

      ! Second Time Step
      imax = nint((ttot - tdiv) / tstep2)

      do i = 1, imax

         step = step + 1
         t1 = t2
         t2 = tdiv + real(i) * tstep2

         if (bextr == 0) then
            omeg = 0._dp
         else
            omeg = log(f0 / ft) / tstep2
         end if

         call trans_calc(0, tstep2, af, tpow1, step, t2)

      end do

   end subroutine rod_eject

   !****************************************************************************!

   subroutine rod_eject_th()

      !
      ! Purpose:
      !    To perform rod ejection simulation with TH feedback
      !

      use sdata, only: ng, nnod, sigr, nf, ttot, tdiv, tstep1, tstep2, ke, &
      tfm, ppow, m, ftem, mtem, bpos, nb, ibeta, &
      f0, ft, c0, tbeta, omeg, npow, ctbeta, nmat, l, th_niter
      use io, only: ounit, bextr, bxtab, scr, bvtk, print_vtk
      use hdf5_output, only: hdf5_is_active, hdf5_write_step
      use xsec, only: xs_updt
      use cmfd, only: outer_tr, outer, outer_ad
      use th, only : th_iter, par_ave, par_max
      use control, only: print_tail

      implicit none

      real(dp), dimension(nnod, ng) :: af                                      ! adjoint flux

      real(dp) :: rho
      real(dp) :: t1, t2
      real(dp) :: tpow1
      integer :: n, i, j, imax, step
      real(dp) :: tf, tm, mtf, mtm
      character(len = 20) :: steady_name

      ! Allocate precusor density
      allocate (c0(nnod, nf))

      ! Allocate Frequency transformation constant
      allocate (omeg(nnod, ng))

      ! Allocate node power distribution
      allocate(npow(nnod))

      ! allocate total leakages
      allocate(l(nnod, ng))

      ! Determine th paramters distribution
      call th_iter(th_niter, 0)

      if (scr) then
         write(*, *)
         write(*, *) ' steady state calculation ... done'
      end if

      ! If K-EFF NOT EQUAL TO 1.0
      if (abs(ke - 1._dp) > 1.e-5_dp .and. bxtab == 0) call kne1

      ! Calculate Adjoint flux
      ! NOTE: This adjoint flux is approximation where
      ! the same nodal coupling coefficients in forward calculation are used
      call outer_ad(0)
      af = f0   ! Save adjoint flux to af
      if (scr) then
         write(*, *)
         write(*, *) ' adjoint calculation ... done'
      end if

      ! Calculate forward flux
      call outer(0)
      if (scr) then
         write(*, *)
         write(*, *) ' re-calculate steady state condition ... done'
      end if

      ! Calculate power
      call powtot(f0, tpow1)

      ! Calculate Initial precursor density
      call ipden()

      ! Total beta
      if (bxtab == 0) then
         tbeta = 0.
         do n = 1, nmat
            do j = 1, nf
               tbeta(n) = tbeta(n) + ibeta(j)
            end do
         end do
         ctbeta = tbeta(1)
      else
         tbeta = 0.
         do n = 1, nmat
            do j = 1, nf
               tbeta(n) = tbeta(n) + m(n)%ibeta(j)
            end do
         end do
         call calc_beta(af)
         write(*, *)
         write(*, 1324) ctbeta * 1.e5
      end if

      ! Calculate reactivity
      call reactivity(af, sigr, rho)

      call par_ave(ftem, tf)
      call par_max(tfm(:, 1), mtf)
      call par_ave(mtem, tm)
      call par_max(mtem, mtm)

      if (hdf5_is_active) call hdf5_write_step(0, 0._dp, reactivity = rho)

      if (bvtk == 1) call print_vtk(0)

      ! File output
      write(ounit, *)
      write(ounit, *) " TRANSIENT RESULTS :"
      write(ounit, *)
      write(ounit, *) " Step  Time(s)  React.($)   Rel. Power   Avg. Tm   Max. Tm   Avg. Tf   Max. Tf"
      write(ounit, *) "--------------------------------------------------------------------------------"
      write(ounit, '(I4, F10.3, F10.4, ES15.4, 4F10.2)') 0, 0. , rho, &
      ppow * 0.01, tm - 273.15, mtm - 273.15, tf - 273.15, mtf - 273.15

      ! Terminal output
      write(*, *)
      write(*, *) " TRANSIENT RESULTS :"
      write(*, *)
      write(*, *) " Step  Time(s)  React.($)   Rel. Power   CR Bank Pos. (1-end)"
      write(*, *) "--------------------------------------------------------------"
      write(*, '(I4, F10.3, F10.4, ES15.4, 12F9.2)') 0, 0. , rho, &
      ppow * 0.01, (bpos(n), n = 1, nb)

      ! Start transient calculation
      step = 0
      t2 = 0.
      imax = nint(tdiv / tstep1)

      ! First Time Step
      do i = 1, imax

         step = step + 1
         t1 = t2
         t2 = real(i) * tstep1

         if (i > 1) then
            if (bextr == 0) then
               omeg = 0._dp
            else
               omeg = log(f0 / ft) / tstep1
            end if
         else
            omeg = 0._dp
         end if

         call trans_calc(1, tstep1, af, tpow1, step, t2)

      end do

      ! Second Time Step
      imax = nint((ttot - tdiv) / tstep2)

      do i = 1, imax

         step = step + 1
         t1 = t2
         t2 = tdiv + real(i) * tstep2

         if (bextr == 0) then
            omeg = 0._dp
         else
            omeg = log(f0 / ft) / tstep2
         end if

         call trans_calc(1, tstep2, af, tpow1, step, t2)

      end do

      call print_tail()

      1324 format(2x, 'Core-averaged delayed neutron fraction :', f7.2, ' pcm')

   end subroutine rod_eject_th

   !****************************************************************************!

   subroutine trans_calc(thc, ht, af, tpow1, step, t2)

      !
      ! Purpose:
      !    To perform transient calculation for given time step
      !

      use sdata, only: ng, nnod, sigr, bcon, ftem, mtem, cden, &
      fbpos, bpos, tmove, bspeed, mdir, nb, velo, npow, &
      f0, ft, fst, fs0, omeg, tranw, ix, iy, iz, pow, &
      tfm, zdel, ppow, node_nf, m, mat, dfis, ctbeta, sth, sigrp
      use io, only: ounit, bxtab, bvtk, print_vtk
      use hdf5_output, only: hdf5_is_active, hdf5_write_step
      use xsec, only: xs_updt
      use cmfd, only: outer_tr, outer
      use th, only : th_trans, par_ave, par_max, powdis


      implicit none

      integer, intent(in) :: thc            !T-H indicator
      real(dp), intent(in) :: ht, tpow1, t2
      integer, intent(in) :: step
      real(dp), dimension(:, :), intent(in) :: af             ! adjoint flux

      real(dp) :: rho
      real(dp) :: tpow2
      integer :: n, g
      logical :: maxi   ! Maximum Outer Iteration Reached?

      real(dp), dimension(nnod) :: pline       ! Linear power density
      real(dp) :: xppow
      real(dp) :: tf, tm, mtf, mtm
      logical :: first = .true.

      if (first) then
         allocate(ft(nnod, ng), fst(nnod))
         allocate(dfis(nnod))
         allocate(sigrp(nnod, ng))
         first = .false.
      end if

      ! Rod bank changes
      do n = 1, nb
         if (mdir(n) == 1) then   ! If CR moving down
            if (t2 - tmove(n) > 1.e-5_dp .and. fbpos(n) - bpos(n) < 1.e-5_dp) then
               bpos(n) = bpos(n) - ht * bspeed(n)
               if (bpos(n) < fbpos(n)) bpos(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            end if
         else if (mdir(n) == 2) then ! If CR moving up
            if (t2 - tmove(n) > 1.e-5_dp .and. fbpos(n) - bpos(n) > 1.e-5_dp) then
               bpos(n) = bpos(n) + ht * bspeed(n)
               if (bpos(n) > fbpos(n)) bpos(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            end if
         else
            continue
         end if
      end do

      ! Calculate xsec after pertubation
      call xs_updt(bcon, ftem, mtem, cden, bpos)

      ! Modify removal xsec
      sigrp = sigr    ! Save sigr to sigrp
      if (bxtab == 0) then
         do g = 1, ng
            do n = 1, nnod
               sigr(n, g) = sigr(n, g) + 1._dp / (sth * velo(g) * ht) + omeg(n, g) / velo(g)
            end do
         end do
      else
         do g = 1, ng
            do n = 1, nnod
               sigr(n, g) = sigr(n, g) + 1._dp / (sth * m(mat(n))%velo(g) * ht) &
               + omeg(n, g) / m(mat(n))%velo(g)
            end do
         end do
      end if

      ! Save the previous fluxes and fission source
      ft = f0
      fst = fs0

      ! Transient calculation
      call outer_tr(ht, maxi)

      ! Update precursor density
      call upden(ht)

      ! Calculate power
      call powtot(f0, tpow2)

      ! Calculate reactivity
      call reactivity(af, sigrp, rho)

      ! Calculate node power distribution
      call powdis(npow)

      if (thc == 1) then
         ! Power change
         xppow = ppow * tpow2 / tpow1 * 0.01_dp

         ! Calculate linear power density for each nodes (W/cm)
         do n = 1, nnod
            pline(n) = npow(n) * pow * xppow &
            / (node_nf(ix(n), iy(n)) * zdel(iz(n)))     ! Linear power density (W/cm)
         end do

         ! TH transient
         call th_trans(pline, ht)

         call par_ave(ftem, tf)
         call par_max(tfm(:, 1), mtf)
         call par_ave(mtem, tm)
         call par_max(mtem, mtm)
      end if

      if (thc == 1) then
         write(*, '(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, rho / ctbeta, &
         xppow, (bpos(n), n = 1, nb)

         if (maxi) then
            write(ounit, '(I4, F10.3, F10.4, ES15.4, 4F10.2, A35)') step, t2, rho / ctbeta, &
            xppow, tm - 273.15, mtm - 273.15, tf - 273.15, mtf - 273.15, 'OUTER ITERATION DID NOT CONVERGE'
         else
            write(ounit, '(I4, F10.3, F10.4, ES15.4, 4F10.2)') step, t2, rho / ctbeta, &
            xppow, tm - 273.15, mtm - 273.15, tf - 273.15, mtf - 273.15
         end if
      else
         write(*, '(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, rho / ctbeta, &
         tpow2 / tpow1, (bpos(n), n = 1, nb)

         if (maxi) then
            write(ounit, '(I4, F10.3, F10.4, ES15.4, A35)') step, t2, rho / ctbeta, &
            tpow2 / tpow1, 'OUTER ITERATION DID NOT CONVERGE'
         else
            write(ounit, '(I4, F10.3, F10.4, ES15.4)') step, t2, rho / ctbeta, &
            tpow2 / tpow1
         end if
      end if

      if (hdf5_is_active) then
         if (thc == 1) then
            call hdf5_write_step(step, t2, reactivity = rho, power_w = pow * xppow)
         else
            call hdf5_write_step(step, t2, reactivity = rho)
         end if
      end if

      if (bvtk == 1) then
         call print_vtk(step)
      end if

      if (maxi) tranw = .true.


   end subroutine trans_calc

   !****************************************************************************!

   subroutine kne1()

      !
      ! Purpose:
      !    To adjuts the Keff to 1.0 if it is not equal to 1.0
      !

      use sdata, only: ke, xnuf, dnuf, bcon, ftem, mtem, cden, bpos
      use io, only: ounit, scr
      use xsec, only: xs_updt
      use cmfd, only: outer

      implicit none

      integer :: i

      write(ounit, *)
      write(ounit, '(A46,F9.6)') '  INITIAL MULTIPLICATION EFFECTIVE (K-EFF) = ', ke
      write(ounit, *) '  WARNING: THE STEADY STATE K-EFF IS NOT EQUAL TO 1.0'
      write(ounit, *) '  AND NOW IT IS FORCED TO 1.0 BY MODIFYING THE nu*sigf CROSS SECTIONS '
      write(ounit, *)
      if (scr) then
         write(*, *)
         write(*, *) ' steady state k-eff not equal to one, force it to one ... done'
      end if
      do i = 1, 10
         xnuf = xnuf / ke
         dnuf = dnuf / ke
         call xs_updt(bcon, ftem, mtem, cden, bpos)
         call outer(0)
         if (abs(ke - 1._dp) < 1.e-5_dp) exit
      end do
      if (i == 10) stop "K-EFF STILL NOT EQUAL TO ONE. KOMODO IS STOPPING"


   end subroutine kne1


   !******************************************************************************!

   subroutine powtot (fx, tpow)

      !
      ! Purpose:
      !    To calculate total power density
      !


      use sdata, only: ng, nnod, sigf, vdel

      implicit none

      real(dp), dimension(:, :), intent(in) :: fx
      real(dp), intent(out) :: tpow

      real(dp), dimension(nnod) :: p
      integer :: g, n
      real(dp) :: pow

      p = 0._dp
      do g = 1, ng
         do n = 1, nnod
            pow = fx(n, g) * sigf(n, g) * vdel(n)
            if (pow < 0.) pow = 0.
            p(n) = p(n) + pow
         end do
      end do


      tpow = 0._dp
      do n = 1, nnod
         tpow = tpow + p(n)
      end do

   end subroutine powtot

   !****************************************************************************!

   subroutine ipden()

      !
      ! Purpose:
      !    Calculate Initial precursor density
      !

      use sdata, only: nnod, nf, fs0, c0, ibeta, lamb, m, mat, nuf, ng
      use io, only: bxtab

      implicit none

      integer :: n, j
      real(dp) :: blamb

      if (bxtab == 1) then
         do n = 1, nnod
            do j = 1, nf
               if (nuf(n, ng) > 0.) then  !If it is fuel
                  blamb = m(mat(n))%ibeta(j) / m(mat(n))%lamb(j)
                  c0(n, j) = blamb * fs0(n)
               else
                  c0(n, j) = 0.
               end if
            end do
         end do
      else
         do n = 1, nnod
            do j = 1, nf
               blamb = ibeta(j) / lamb(j)
               c0(n, j) = blamb * fs0(n)
            end do
         end do
      end if


   end subroutine ipden

   !****************************************************************************!

   subroutine upden(ht)

      !
      ! Purpose:
      !    To update precursor density
      !

      use sdata, only: nnod, nf, fs0, fst, c0, ibeta, lamb, m, mat, nuf, ng
      use io, only: bxtab

      implicit none

      real(dp), intent(in) :: ht
      real(dp) :: a1, a2, pxe
      integer :: n, i

      if (bxtab == 1) then
         do i = 1, nf
            do n = 1, nnod
               if (nuf(n, ng) > 0.) then  !If it is fuel
                  pxe = exp(-m(mat(n))%lamb(i) * ht)
                  a1 = (1._dp - pxe) / (m(mat(n))%lamb(i) * ht)
                  a2 = 1._dp - a1
                  a1 = a1 - pxe
                  c0(n, i) = c0(n, i) * pxe + m(mat(n))%ibeta(i) / m(mat(n))%lamb(i) &
                  * (a1 * fst(n) + a2 * fs0(n))
               end if
            end do
         end do
      else
         do i = 1, nf
            pxe = exp(-lamb(i) * ht)
            a1 = (1._dp - pxe) / (lamb(i) * ht)
            a2 = 1._dp - a1
            a1 = a1 - pxe
            do n = 1, nnod
               c0(n, i) = c0(n, i) * pxe + ibeta(i) / lamb(i) &
               * (a1 * fst(n) + a2 * fs0(n))
            end do
         end do
      end if


   end subroutine upden

   !****************************************************************************!

   subroutine reactivity(af, sigrp, rho)

      !
      ! Purpose:
      !    To calculate dynamic reactivity
      !

      use sdata, only: nnod, ng, f0, sigs, chi, mat, fs0, vdel, l
      use nodal, only: lxyz

      implicit none

      real(dp), dimension(:, :), intent(in) :: af
      real(dp), dimension(:, :), intent(in) :: sigrp
      real(dp), intent(out) :: rho

      integer :: n, g, h
      real(dp), dimension(nnod) :: scg
      real(dp) :: rem, lea, src, fde, l1, l2, l3

      src = 0. ;
      rem = 0. ;
      lea = 0. ;
      fde = 0.
      do g = 1, ng
         scg = 0.
         do h = 1, ng
            do n = 1, nnod
               if (g /= h) scg(n) = scg(n) + sigs(n, h, g) * f0(n, h)
            end do
         end do
         do n = 1, nnod
            call lxyz(n, g, l1, l2, l3)
            l(n, g) = l1 + l2 + l3
            src = src + af(n, g) * (scg(n) + chi(mat(n), g) * fs0(n)) * vdel(n)
            rem = rem + af(n, g) * sigrp(n, g) * f0(n, g) * vdel(n)
            lea = lea + af(n, g) * l(n, g) * vdel(n)
            fde = fde + af(n, g) * chi(mat(n), g) * fs0(n) * vdel(n)
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

      use sdata, only: ng, nnod, m, mat, chi, nf, f0, nuf, ctbeta
      use cmfd, only: integrate
      use io, only: ounit

      implicit none

      real(dp), dimension(:, :), intent(in) :: af

      integer :: n, i, g
      real(dp), dimension(nnod) :: vdum, vdum2
      real(dp) :: f

      ! Calculate F
      vdum = 0.
      do g = 1, ng
         do n = 1, nnod
            vdum(n) = vdum(n) + nuf(n, g) * f0(n, g)
         end do
      end do

      vdum2 = 0.
      do g = 1, ng
         do n = 1, nnod
            vdum2(n) = vdum2(n) + chi(mat(n), g) * vdum(n) * af(n, g)
         end do
      end do

      f = integrate(vdum2)

      ! Calculate Delayed neutron fraction (beta)
      ctbeta = 0._dp
      do i = 1, nf
         vdum2 = 0.
         do g = 1, ng
            do n = 1, nnod
               vdum2(n) = vdum2(n) + chi(mat(n), g) * m(mat(n))%ibeta(i) * vdum(n) * af(n, g)
            end do
         end do
         ctbeta = ctbeta + integrate(vdum2) / f
      end do

      write(ounit, *)
      write(ounit, 1344) ctbeta * 1.e5

      1344 format ('  CORE AVERAGED DELAYED NEUTRON FRACTION: ', f7.2, ' PCM')


   end subroutine calc_beta


end module trans
