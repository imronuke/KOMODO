module xsec

   use sdata, only: dp
   use io, only: bbcon, bftem, bmtem, bcden, bcrod, bcbcs, bxtab, ounit
   implicit none
   save

   !Define + and - operators for XBRANCH type addition and substitution respectively
   interface operator (+)
      module procedure bradd
   end interface
   interface operator (-)
      module procedure brsubst
   end interface
   interface operator (*)
      module procedure brrealmult
   end interface

   contains

   !******************************************************************************!

   subroutine xs_updt (xbcon, xftem, xmtem, xcden, xbpos)
      !
      ! Purpose:
      !    To update XS for given TH paramaters and rod position
      !

      use sdata, only : nnod, mat, sigtr, siga, nuf, sigf, sigs, dc, &
      get_time, xs_time


      real(dp), intent(in) :: xbcon  ! Provided Boron Concentration
      real(dp), dimension(:), intent(in) :: xftem  ! Provided fuel temperature
      real(dp), dimension(:), intent(in) :: xmtem  ! Provided moderator temperature
      real(dp), dimension(:), intent(in) :: xcden  ! Provided coolant density
      real(dp), dimension(:), intent(in) :: xbpos  ! Provided control rod bank position

      integer :: n

      real(dp) :: st, fn

      st = get_time()

      if (bxtab == 0) then
         call base_updt()
         if (bbcon == 1 .or. bcbcs == 1) call bcon_updt(xbcon)
         if (bftem == 1) call ftem_updt(xftem)
         if (bmtem == 1) call mtem_updt(xmtem)
         if (bcden == 1) call cden_updt(xcden)
         if (bcrod == 1) call crod_updt(xbpos)
         call dsigr_updt()
      else
         do n = 1, nnod
            call brinterp(0, mat(n), xcden(n), xbcon, xftem(n), xmtem(n), sigtr(n, :), &
            siga(n, :), nuf(n, :), sigf(n, :), sigs(n, :, :), dc(n, :, :))
         end do
         if (bcrod == 1) call crod_tab_updt(xbpos)
         call dsigr_updt()
      end if

      call check_xs()

      fn = get_time()
      xs_time = xs_time + (fn - st)

   end subroutine xs_updt

   !******************************************************************************!

   subroutine check_xs ()
      !
      ! Purpose:
      !    To check errors in xsec
      !

      use sdata, only : d, sigr, nuf, chi, sigs, nnod, ng, mat, nmat



      integer :: n, g, h

      do g = 1, ng
         do n = 1, nnod
            if (d(n, g) < 1.e-20) then
               write(*, *)
               write(*, 3541) mat(n), g
               write(*, *) ' DIFFUSION COEF. IS CLOSE TO ZERO OR NEGATIVE'
               write(ounit, 3541) mat(n), g
               write(ounit, *) ' DIFFUSION COEF. IS CLOSE TO ZERO OR NEGATIVE'
               stop
            end if
            if (sigr(n, g) < 0.) then
               write(*, *)
               write(*, 3542) mat(n), g
               write(*, *) ' REMOVAL XS IS  NEGATIVE'
               write(ounit, 3542) mat(n), g
               write(ounit, *) ' REMOVAL XS IS NEGATIVE'
               stop
            end if
            if (nuf(n, g) < 0.) then
               write(*, *)
               write(*, 3544) mat(n), g
               write(*, *) ' NU*FISSION XS IS  NEGATIVE'
               write(ounit, 3544) mat(n), g
               write(ounit, *) ' NU*FISSION XS IS NEGATIVE'
               stop
            end if
         end do
         do n = 1, nmat
            if (chi(n, g) < 0.) then
               write(*, *)
               write(*, 3543) n, g
               write(*, *) ' FISSION SPECTRUM IS  NEGATIVE'
               write(ounit, 3543) n, g
               write(ounit, *) ' FISSION SPECTRUM IS NEGATIVE'
               stop
            end if
         end do
         do h = 1, ng
            do n = 1, nnod
               if (sigs(n, g, h) < 0.) then
                  write(*, *)
                  write(*, 3545) mat(n), g, h
                  write(*, *) ' SCATTERING XS IS  NEGATIVE'
                  write(ounit, 3545) mat(n), g, h
                  write(ounit, *) ' SCATTERING XS IS NEGATIVE'
                  stop
               end if
            end do
         end do
      end do

      3541 format (2x, 'ERROR IN THE DIFFUSION COEFFICIENT ', &
      'FOR MATERIAL NUMBER : ', i3, ' ENERGY GROUP : ', i2)
      3542 format (2x, 'ERROR IN THE REMOVAL CROSS SECTION ', &
      'FOR MATERIAL NUMBER : ', i3, ' ENERGY GROUP : ', i2)
      3543 format (2x, 'ERROR IN THE FISSION SPECTRUM ', &
      'FOR MATERIAL NUMBER : ', i3, ' ENERGY GROUP : ', i2)
      3544 format (2x, 'ERROR IN THE NU*FISSION XSEC ', &
      'FOR MATERIAL NUMBER : ', i3, ' ENERGY GROUP : ', i2)
      3545 format (2x, 'ERROR IN THE SCATTERING CROSS SECTION ', &
      'FOR MATERIAL NUMBER : ', i3, ' ENERGY GROUP : ', i2, &
      ' TO ENERGY GROUP : ', i2)


   end subroutine check_xs

   !******************************************************************************!

   subroutine base_updt ()
      !
      ! Purpose:
      !    To update current XS to base XS
      !


      use sdata, only: nnod, sigtr, siga, nuf, sigf, sigs, mat, &
      xsigtr, xsiga, xnuf, xsigf, xsigs


      integer :: n

      do n = 1, nnod
         sigtr(n, 1:) = xsigtr(mat(n), 1:)
         siga (n, 1:) = xsiga (mat(n), 1:)
         nuf (n, 1:) = xnuf (mat(n), 1:)
         sigf (n, 1:) = xsigf (mat(n), 1:)
         sigs (n, 1:, 1:) = xsigs (mat(n), 1:, 1:)
      end do


   end subroutine base_updt

   !******************************************************************************!

   subroutine dsigr_updt ()
      !
      ! Purpose:
      !    To update diffusion constant and removal XS
      !


      use sdata, only: nnod, ng, sigtr, siga, &
      sigs, d, sigr


      integer :: i, g, h
      real(dp) :: dum

      do i = 1, nnod
         do g = 1, ng
            if (sigtr(i, g) < 1.e-5) stop "Negative diffusion coefficient encountered"
            d(i, g) = 1._dp / (3._dp * sigtr(i, g))
            dum = 0.
            do h = 1, ng
               if (g /= h) dum = dum + sigs(i, g, h)
            end do
            sigr(i, g) = siga(i, g) + dum
         end do
      end do

   end subroutine dsigr_updt

   !******************************************************************************!

   subroutine crod_updt (bpos)
      !
      ! Purpose: TO UPDATE AND CALCUALTE VOLUME WEIGHTED HOMOGENIZED CX FOR RODDED NODE
      !

      use sdata, only: ng, nxx, nyy, nzz, xyz, zdel, mat, &
      sigtr, siga, nuf, sigf, sigs, &
      dsigtr, dsiga, dnuf, dsigf, dsigs, &
      coreh, fbmap, pos0, ssize


      real(dp), dimension(:), intent(in) :: bpos

      integer :: i, j, k, g, h
      real(dp) :: rodh, vfrac
      real(dp) :: dum

      ! For each node
      do j = 1, nyy
         do i = 1, nxx
            if (fbmap(i, j) > 0) then
               !!!(rodh -> posistion the tip of the control rod the top of core)
               rodh = coreh - pos0 - bpos(fbmap(i, j)) * ssize
               dum = 0._dp
               do k = nzz, 1, -1
                  ! For partially rodded node, get volume weighted homogenized CX (0 < vfrac < 1._DP)
                  if (rodh >= dum .and. rodh <= dum + zdel(k)) then   ! If this node partially rodded
                     vfrac = (rodh - dum) / zdel(k)
                     sigtr(xyz(i, j, k), :) = sigtr(xyz(i, j, k), :) + &
                     vfrac * dsigtr(mat(xyz(i, j, k)), :)
                     siga(xyz(i, j, k), :) = siga(xyz(i, j, k), :) + &
                     vfrac * dsiga(mat(xyz(i, j, k)), :)
                     nuf(xyz(i, j, k), :) = nuf(xyz(i, j, k), :) + &
                     vfrac * dnuf(mat(xyz(i, j, k)), :)
                     sigf(xyz(i, j, k), :) = sigf(xyz(i, j, k), :) + &
                     vfrac * dsigf(mat(xyz(i, j, k)), :)
                     sigs(xyz(i, j, k), :, :) = sigs(xyz(i, j, k), :, :) + &
                     vfrac * dsigs(mat(xyz(i, j, k)), :, :)
                     exit
                  end if
                  ! For fully rodded node, vfrac = 1.
                  sigtr(xyz(i, j, k), :) = sigtr(xyz(i, j, k), :) + dsigtr(mat(xyz(i, j, k)), :)
                  siga(xyz(i, j, k), :) = siga(xyz(i, j, k), :) + dsiga(mat(xyz(i, j, k)), :)
                  nuf(xyz(i, j, k), :) = nuf(xyz(i, j, k), :) + dnuf(mat(xyz(i, j, k)), :)
                  sigf(xyz(i, j, k), :) = sigf(xyz(i, j, k), :) + dsigf(mat(xyz(i, j, k)), :)
                  sigs(xyz(i, j, k), :, :) = sigs(xyz(i, j, k), :, :) + dsigs(mat(xyz(i, j, k)), :, :)

                  dum = dum + zdel(k)
               end do
               ! if negative CX found, Surpress CX to zero and calculate D and sigr
               do k = nzz, 1, -1
                  do g = 1, ng
                     if (siga(xyz(i, j, k), g) < 0.) siga(xyz(i, j, k), g) = 0.
                     if (nuf(xyz(i, j, k), g) < 0.) nuf(xyz(i, j, k), g) = 0.
                     if (sigf(xyz(i, j, k), g) < 0.) sigf(xyz(i, j, k), g) = 0.
                     do h = 1, ng
                        if (sigs(xyz(i, j, k), g, h) < 0.) sigs(xyz(i, j, k), g, h) = 0.
                     end do
                  end do
               end do
            end if
         end do
      end do


   end subroutine crod_updt

   !******************************************************************************!

   subroutine crod_tab_updt (bpos)
      !
      ! Purpose: TO UPDATE AND CALCUALTE VOLUME WEIGHTED HOMOGENIZED CX FOR RODDED NODE
      !

      use sdata, only: ng, nxx, nyy, nzz, xyz, zdel, mat, &
      sigtr, siga, nuf, sigf, sigs, dc, &
      coreh, fbmap, pos0, ssize, m, &
      nnod, cden, ftem, mtem, bcon


      real(dp), dimension(:), intent(in) :: bpos

      integer :: i, j, k, g, h
      real(dp) :: rodh, vfrac
      real(dp) :: dum

      integer :: nd

      real(dp), dimension(nnod, ng) :: rsigtr, rsiga, rnuf, rsigf
      real(dp), dimension(nnod, ng, ng) :: rsigs
      real(dp), dimension(nnod, ng, 6) :: rdc

      do j = 1, nyy
         do i = 1, nxx
            if (fbmap(i, j) > 0) then
               !(rodh -> posistion the tip of the control rod the top of core)
               rodh = coreh - pos0 - bpos(fbmap(i, j)) * ssize
               dum = 0._dp
               do k = nzz, 1, -1
                  nd = xyz(i, j, k)                                  ! Node number
                  if (m(mat(nd))%trod == 1) then
                     call brinterp(1, mat(nd), cden(nd), bcon, ftem(nd), &
                     mtem(nd), rsigtr(nd, :), rsiga(nd, :), rnuf(nd, :), &
                     rsigf(nd, :), rsigs(nd, :, :), rdc(nd, :, :))
                  else
                     write(101, 1671) fbmap(i, j), mat(nd)
                     write(*, 1671) fbmap(i, j), mat(nd)
                     stop
                     1671 format (2x, 'CONTROL ROD BANK NUMBER ', i4, &
                     ' COINCIDES WITH MATERIAL NUMBER ', i4, &
                     ' THAT DOES NOT HAVE CONTROL ROD DATA IN XTAB FILE')
                  end if
                  ! For partially rodded node, get volume weighted homogenized CX (0 < vfrac < 1._DP)
                  if (rodh >= dum .and. rodh <= dum + zdel(k)) then   ! If this node partially rodded
                     vfrac = (rodh - dum) / zdel(k)                    ! Rodded fraction
                     sigtr(nd, :) = (1._dp - vfrac) * sigtr(nd, :) + &
                     vfrac * rsigtr(nd, :)
                     siga(nd, :) = (1._dp - vfrac) * siga(nd, :) + &
                     vfrac * rsiga(nd, :)
                     nuf(nd, :) = (1._dp - vfrac) * nuf(nd, :) + &
                     vfrac * rnuf(nd, :)
                     sigf(nd, :) = (1._dp - vfrac) * sigf(nd, :) + &
                     vfrac * rsigf(nd, :)
                     sigs(nd, :, :) = (1._dp - vfrac) * sigs(nd, :, :) + &
                     vfrac * rsigs(nd, :, :)
                     dc(nd, :, :) = (1._dp - vfrac) * dc(nd, :, :) + &
                     vfrac * rdc(nd, :, :)
                     exit
                  end if
                  ! For fully rodded node, vfrac = 1.
                  sigtr(nd, :) = rsigtr(nd, :)
                  siga(nd, :) = rsiga(nd, :)
                  nuf(nd, :) = rnuf(nd, :)
                  sigf(nd, :) = rsigf(nd, :)
                  sigs(nd, :, :) = rsigs(nd, :, :)
                  dc(nd, :, :) = rdc(nd, :, :)

                  dum = dum + zdel(k)
               end do
               ! if negative CX found, Surpress CX to zero  and calculate D and sigr
               do k = nzz, 1, -1
                  nd = xyz(i, j, k)
                  do g = 1, ng
                     if (siga(nd, g) < 0.) siga(nd, g) = 0.
                     if (nuf(nd, g) < 0.) nuf(nd, g) = 0.
                     if (sigf(nd, g) < 0.) sigf(nd, g) = 0.
                     do h = 1, ng
                        if (sigs(nd, g, h) < 0.) sigs(nd, g, h) = 0.
                     end do
                     do h = 1, 6
                        if (dc(nd, g, h) < 0.) dc(nd, g, h) = 0.
                     end do
                  end do
               end do
            end if
         end do
      end do

   end subroutine crod_tab_updt

   !******************************************************************************!

   subroutine bcon_updt (bcon)

      !
      ! Purpose:
      !    To update CX for given boron concentration

      use sdata, only: nnod, ng, sigtr, siga, nuf, sigf, sigs, mat, &
      csigtr, csiga, cnuf, csigf, csigs, rbcon


      real(dp), intent(in) :: bcon
      integer :: n, g, h

      do n = 1, nnod
         do g = 1, ng
            sigtr(n, g) = sigtr(n, g) + csigtr(mat(n), g) * (bcon - rbcon)
            siga(n, g) = siga(n, g) + csiga(mat(n), g) * (bcon - rbcon)
            nuf(n, g) = nuf(n, g) + cnuf(mat(n), g) * (bcon - rbcon)
            sigf(n, g) = sigf(n, g) + csigf(mat(n), g) * (bcon - rbcon)
            do h = 1, ng
               sigs(n, g, h) = sigs(n, g, h) + csigs(mat(n), g, h) * (bcon - rbcon)
            end do
         end do
      end do


   end subroutine bcon_updt

   !******************************************************************************!

   subroutine ftem_updt (ftem)

      !
      ! Purpose:
      !    To update CX for given fuel temp

      use sdata, only: nnod, ng, sigtr, siga, nuf, sigf, sigs, mat, &
      fsigtr, fsiga, fnuf, fsigf, fsigs, rftem


      real(dp), dimension(:), intent(in) :: ftem
      integer :: n, g, h


      do n = 1, nnod
         do g = 1, ng
            sigtr(n, g) = sigtr(n, g) + fsigtr(mat(n), g) * (sqrt(ftem(n))  - sqrt(rftem))
            siga(n, g) = siga(n, g) + fsiga(mat(n), g) * (sqrt(ftem(n)) - sqrt(rftem))
            nuf(n, g) = nuf(n, g) + fnuf(mat(n), g) * (sqrt(ftem(n)) - sqrt(rftem))
            sigf(n, g) = sigf(n, g) + fsigf(mat(n), g) * (sqrt(ftem(n)) - sqrt(rftem))
            do h = 1, ng
               sigs(n, g, h) = sigs(n, g, h) + fsigs(mat(n), g, h) * (sqrt(ftem(n)) - sqrt(rftem))
            end do
         end do
      end do

   end subroutine ftem_updt

   !******************************************************************************!

   subroutine mtem_updt (mtem)

      !
      ! Purpose:
      !    To update CX for given moderator temperature

      use sdata, only: nnod, ng, sigtr, siga, nuf, sigf, sigs, mat, &
      msigtr, msiga, mnuf, msigf, msigs, rmtem


      real(dp), dimension(:), intent(in) :: mtem
      integer :: n, g, h


      do n = 1, nnod
         do g = 1, ng
            sigtr(n, g) = sigtr(n, g) + msigtr(mat(n), g) * (mtem(n) - rmtem)
            siga(n, g) = siga(n, g) + msiga(mat(n), g) * (mtem(n) - rmtem)
            nuf(n, g) = nuf(n, g) + mnuf(mat(n), g) * (mtem(n) - rmtem)
            sigf(n, g) = sigf(n, g) + msigf(mat(n), g) * (mtem(n) - rmtem)
            do h = 1, ng
               sigs(n, g, h) = sigs(n, g, h) + msigs(mat(n), g, h) * (mtem(n) - rmtem)
            end do
         end do
      end do


   end subroutine mtem_updt

   !******************************************************************************!

   subroutine cden_updt (cden)

      !
      ! Purpose:
      !    To update CX for given coolant density

      use sdata, only: nnod, ng, sigtr, siga, nuf, sigf, sigs, mat, &
      lsigtr, lsiga, lnuf, lsigf, lsigs, rcden


      real(dp), dimension(:), intent(in) :: cden
      integer :: n, g, h


      do n = 1, nnod
         do g = 1, ng
            sigtr(n, g) = sigtr(n, g) + lsigtr(mat(n), g) * (cden(n) - rcden)
            siga(n, g) = siga(n, g) + lsiga(mat(n), g) * (cden(n) - rcden)
            nuf(n, g) = nuf(n, g) + lnuf(mat(n), g) * (cden(n) - rcden)
            sigf(n, g) = sigf(n, g) + lsigf(mat(n), g) * (cden(n) - rcden)
            do h = 1, ng
               sigs(n, g, h) = sigs(n, g, h) + lsigs(mat(n), g, h) * (cden(n) - rcden)
            end do
         end do
      end do


   end subroutine cden_updt

   !******************************************************************************!

   subroutine brinterp (rod, mn, xcden, xbcon, xftem, xmtem, sigtr, siga, nuf, &
   sigf, sigs, dc)
      !Purpose: To interpolate the xsec data from xtab file for given bcon,
      ! ftem, mtem and cden

      use sdata, only: m, xbranch, ng


      integer, intent(in) :: rod, mn  ! CR indicator and material number
      real(dp), intent(in) :: xbcon, xftem, xmtem, xcden  ! TH Parameters
      real(dp), dimension(:), intent(out) :: sigtr, siga, nuf, sigf
      real(dp), dimension(:, :), intent(out) :: dc, sigs

      integer, parameter :: nx = 8
      type(xbranch), dimension(nx) :: xs   !Temporary xsec for interpolation
      integer :: s, t, u, v, mx
      integer :: s1, s2, t1, t2, u1, u2, v1, v2  ! Dimnesion Position of the given parameters
      integer :: i
      real(dp) :: radx

      ! Set to default
      s1 = 1;
      s2 = 1;
      t1 = 1;
      t2 = 1;
      u1 = 1;
      u2 = 1;
      v1 = 1;
      v2 = 1

      ! Get 2 closest points for interpolation
      ! FOR COOLANT DENSITY
      if (m(mn)%nd > 1) then
         mx = m(mn)%nd
         if (xcden >= m(mn)%pd(1) .and. xcden <= m(mn)%pd(mx)) then
            do s = 2, mx
               if (xcden >= m(mn)%pd(s - 1) .and. xcden <= m(mn)%pd(s)) then
                  s1 = s - 1
                  s2 = s
                  exit
               end if
            end do
         else if (xcden < m(mn)%pd(1) .and. (m(mn)%pd(1) - xcden) / m(mn)%pd(1) < 0.2) then
            s1 = 1
            s2 = 2
         else if (xcden > m(mn)%pd(mx) .and. (xcden - m(mn)%pd(mx)) / m(mn)%pd(mx) < 0.2) then
            s1 = mx - 1
            s2 = mx
         else
            write(ounit, 1567) xcden
            write(*, 1567) xcden
            stop
            1567 format(2x, '  ERROR: COOLANT DENSITY ', f7.3, ' g/cm3 &
            & is out of the range of the branch parameter')
         end if
      end if

      ! FOR BORON CONCENTRATION
      if (m(mn)%nb > 1) then
         mx = m(mn)%nb
         if (xbcon >= m(mn)%pb(1) .and. xbcon <= m(mn)%pb(mx)) then
            do t = 2, mx
               if (xbcon >= m(mn)%pb(t - 1) .and. xbcon <= m(mn)%pb(t)) then
                  t1 = t - 1
                  t2 = t
                  exit
               end if
            end do
         else if (xbcon < m(mn)%pb(1) .and. (m(mn)%pb(1) - xbcon) < 100.) then
            t1 = 1
            t2 = 2
         else if (xbcon > m(mn)%pb(mx) .and. (xbcon - m(mn)%pb(mx)) < 100.) then
            t1 = mx - 1
            t2 = mx
         else
            write(ounit, 1568) xbcon
            write(*, 1568) xbcon
            stop
            1568 format(2x, '  ERROR: BORON CONCENTRATION ', f8.1, ' PPM &
            & is out of the range of the branch parameter')
         end if
      end if

      ! FOR FUEL TEMPERATURE
      if (m(mn)%nf > 1) then
         mx = m(mn)%nf
         if (xftem >= m(mn)%pf(1) .and. xftem <= m(mn)%pf(mx)) then
            do u = 2, mx
               if (xftem >= m(mn)%pf(u - 1) .and. xftem <= m(mn)%pf(u)) then
                  u1 = u - 1
                  u2 = u
                  exit
               end if
            end do
         else if (xftem < m(mn)%pf(1) .and. (m(mn)%pf(1) - xftem) / m(mn)%pf(1) < 0.2) then
            u1 = 1
            u2 = 2
         else if (xftem > m(mn)%pf(mx) .and. (xftem - m(mn)%pf(mx)) / m(mn)%pf(mx) < 0.2) then
            u1 = mx - 1
            u2 = mx
         else
            write(ounit, 1570) xftem
            write(*, 1570) xftem
            stop
            1570 format(2x, '  ERROR: FUEL TEMPERATURE ', f7.1, ' K &
            & is out of the range of the branch parameter')
         end if
      end if

      ! FOR MODERATOR TEMPERATURE
      if (m(mn)%nm > 1) then
         mx = m(mn)%nm
         if (xmtem >= m(mn)%pm(1) .and. xmtem <= m(mn)%pm(mx)) then
            do v = 2, mx
               if (xmtem >= m(mn)%pm(v - 1) .and. xmtem <= m(mn)%pm(v)) then
                  v1 = v - 1
                  v2 = v
                  exit
               end if
            end do
         else if (xmtem < m(mn)%pm(1) .and. (m(mn)%pm(1) - xmtem) / m(mn)%pm(1) < 0.2) then
            v1 = 1
            v2 = 2
         else if (xmtem > m(mn)%pm(mx) .and. (xmtem - m(mn)%pm(mx)) / m(mn)%pm(mx) < 0.2) then
            v1 = mx - 1
            v2 = mx
         else
            write(ounit, 1569) xmtem
            write(*, 1569) xmtem
            stop
            1569 format(2x, '  ERROR: MODERATOR TEMPERATURE ', f7.1, ' K &
            & is out of the range of the branch parameter')
         end if
      end if

      !Start doing interpolation
      do i = 1, nx   !Allocate memory to xs
         allocate(xs(i)%sigtr(ng))
         allocate(xs(i)%siga(ng))
         allocate(xs(i)%nuf(ng))
         allocate(xs(i)%sigf(ng))
         allocate(xs(i)%dc(ng, 6))
         allocate(xs(i)%sigs(ng, ng))
      end do

      if (rod == 0) then   ! For Unrodded XSEC
         !interpolation on Moderator Temperature
         if (m(mn)%nm > 1) then
            radx = (xmtem - m(mn)%pm(v1)) / (m(mn)%pm(v2) - m(mn)%pm(v1))
            xs(1) = m(mn)%xsec(s1, t1, u1, v1) + &
            (radx * (m(mn)%xsec(s1, t1, u1, v2) - m(mn)%xsec(s1, t1, u1, v1)))
            xs(2) = m(mn)%xsec(s1, t1, u2, v1) + &
            (radx * (m(mn)%xsec(s1, t1, u2, v2) - m(mn)%xsec(s1, t1, u2, v1)))
            xs(3) = m(mn)%xsec(s1, t2, u1, v1) + &
            (radx * (m(mn)%xsec(s1, t2, u1, v2) - m(mn)%xsec(s1, t2, u1, v1)))
            xs(4) = m(mn)%xsec(s1, t2, u2, v1) + &
            (radx * (m(mn)%xsec(s1, t2, u2, v2) - m(mn)%xsec(s1, t2, u2, v1)))
            xs(5) = m(mn)%xsec(s2, t1, u1, v1) + &
            (radx * (m(mn)%xsec(s2, t1, u1, v2) - m(mn)%xsec(s2, t1, u1, v1)))
            xs(6) = m(mn)%xsec(s2, t1, u2, v1) + &
            (radx * (m(mn)%xsec(s2, t1, u2, v2) - m(mn)%xsec(s2, t1, u2, v1)))
            xs(7) = m(mn)%xsec(s2, t2, u1, v1) + &
            (radx * (m(mn)%xsec(s2, t2, u1, v2) - m(mn)%xsec(s2, t2, u1, v1)))
            xs(8) = m(mn)%xsec(s2, t2, u2, v1) + &
            (radx * (m(mn)%xsec(s2, t2, u2, v2) - m(mn)%xsec(s2, t2, u2, v1)))
         else
            xs(1) = m(mn)%xsec(s1, t1, u1, v1)
            xs(2) = m(mn)%xsec(s1, t1, u2, v1)
            xs(3) = m(mn)%xsec(s1, t2, u1, v1)
            xs(4) = m(mn)%xsec(s1, t2, u2, v1)
            xs(5) = m(mn)%xsec(s2, t1, u1, v1)
            xs(6) = m(mn)%xsec(s2, t1, u2, v1)
            xs(7) = m(mn)%xsec(s2, t2, u1, v1)
            xs(8) = m(mn)%xsec(s2, t2, u2, v1)
         end if
         !interpolation on Fuel Temperature
         if (m(mn)%nf > 1) then
            radx = (xftem - m(mn)%pf(u1)) / (m(mn)%pf(u2) - m(mn)%pf(u1))
            xs(1) = xs(1) + (radx * (xs(2) - xs(1)))
            xs(3) = xs(3) + (radx * (xs(4) - xs(3)))
            xs(5) = xs(5) + (radx * (xs(6) - xs(5)))
            xs(7) = xs(7) + (radx * (xs(8) - xs(7)))
         end if
         !interpolation on Boron concentration
         if (m(mn)%nb > 1) then
            radx = (xbcon - m(mn)%pb(t1)) / (m(mn)%pb(t2) - m(mn)%pb(t1))
            xs(1) = xs(1) + (radx * (xs(3) - xs(1)))
            xs(5) = xs(5) + (radx * (xs(7) - xs(5)))
         end if
         !interpolation on coolant density
         if (m(mn)%nd > 1) then
            xs(1) = xs(1) + ((xcden - m(mn)%pd(s1)) / (m(mn)%pd(s2) - m(mn)%pd(s1)) * &
            (xs(5) - xs(1)))
         end if
      else   ! For Rodded XSEC
         !interpolation on Moderator Temperature
         if (m(mn)%nm > 1) then
            radx = (xmtem - m(mn)%pm(v1)) / (m(mn)%pm(v2) - m(mn)%pm(v1))
            xs(1) = m(mn)%rxsec(s1, t1, u1, v1) + &
            (radx * (m(mn)%rxsec(s1, t1, u1, v2) - m(mn)%rxsec(s1, t1, u1, v1)))
            xs(2) = m(mn)%rxsec(s1, t1, u2, v1) + &
            (radx * (m(mn)%rxsec(s1, t1, u2, v2) - m(mn)%rxsec(s1, t1, u2, v1)))
            xs(3) = m(mn)%rxsec(s1, t2, u1, v1) + &
            (radx * (m(mn)%rxsec(s1, t2, u1, v2) - m(mn)%rxsec(s1, t2, u1, v1)))
            xs(4) = m(mn)%rxsec(s1, t2, u2, v1) + &
            (radx * (m(mn)%rxsec(s1, t2, u2, v2) - m(mn)%rxsec(s1, t2, u2, v1)))
            xs(5) = m(mn)%rxsec(s2, t1, u1, v1) + &
            (radx * (m(mn)%rxsec(s2, t1, u1, v2) - m(mn)%rxsec(s2, t1, u1, v1)))
            xs(6) = m(mn)%rxsec(s2, t1, u2, v1) + &
            (radx * (m(mn)%rxsec(s2, t1, u2, v2) - m(mn)%rxsec(s2, t1, u2, v1)))
            xs(7) = m(mn)%rxsec(s2, t2, u1, v1) + &
            (radx * (m(mn)%rxsec(s2, t2, u1, v2) - m(mn)%rxsec(s2, t2, u1, v1)))
            xs(8) = m(mn)%rxsec(s2, t2, u2, v1) + &
            (radx * (m(mn)%rxsec(s2, t2, u2, v2) - m(mn)%rxsec(s2, t2, u2, v1)))
         else
            xs(1) = m(mn)%rxsec(s1, t1, u1, v1)
            xs(2) = m(mn)%rxsec(s1, t1, u2, v1)
            xs(3) = m(mn)%rxsec(s1, t2, u1, v1)
            xs(4) = m(mn)%rxsec(s1, t2, u2, v1)
            xs(5) = m(mn)%rxsec(s2, t1, u1, v1)
            xs(6) = m(mn)%rxsec(s2, t1, u2, v1)
            xs(7) = m(mn)%rxsec(s2, t2, u1, v1)
            xs(8) = m(mn)%rxsec(s2, t2, u2, v1)
         end if
         !interpolation on Fuel Temperature
         if (m(mn)%nf > 1) then
            radx = (xftem - m(mn)%pf(u1)) / (m(mn)%pf(u2) - m(mn)%pf(u1))
            xs(1) = xs(1) + (radx * (xs(2) - xs(1)))
            xs(3) = xs(3) + (radx * (xs(4) - xs(3)))
            xs(5) = xs(5) + (radx * (xs(6) - xs(5)))
            xs(7) = xs(7) + (radx * (xs(8) - xs(7)))
         end if
         !interpolation on Boron concentration
         if (m(mn)%nb > 1) then
            radx = (xbcon - m(mn)%pb(t1)) / (m(mn)%pb(t2) - m(mn)%pb(t1))
            xs(1) = xs(1) + (radx * (xs(3) - xs(1)))
            xs(5) = xs(5) + (radx * (xs(7) - xs(5)))
         end if

         !interpolation on coolant density
         if (m(mn)%nd > 1) then
            xs(1) = xs(1) + ((xcden - m(mn)%pd(s1)) / (m(mn)%pd(s2) - m(mn)%pd(s1)) * &
            (xs(5) - xs(1)))
         end if
      end if
      sigtr = xs(1)%sigtr
      siga = xs(1)%siga
      nuf = xs(1)%nuf
      sigf = xs(1)%sigf
      sigs = xs(1)%sigs
      dc = xs(1)%dc

      do i = 1, nx   !DeAllocate memory to xs
         deallocate(xs(i)%sigtr)
         deallocate(xs(i)%siga)
         deallocate(xs(i)%nuf)
         deallocate(xs(i)%sigf)
         deallocate(xs(i)%dc)
         deallocate(xs(i)%sigs)
      end do


   end subroutine brinterp

   !******************************************************************************!

   function bradd(a, b) result (c)

      ! To perform XBRANCH data type addition

      use sdata, only: xbranch


      type(xbranch), intent(in) :: a, b
      type(xbranch) :: c

      c%sigtr = a%sigtr + b%sigtr
      c%siga = a%siga + b%siga
      c%nuf = a%nuf + b%nuf
      c%sigf = a%sigf + b%sigf
      c%sigs = a%sigs + b%sigs
      c%dc = a%dc + b%dc

   end function bradd

   !******************************************************************************!

   function brsubst(a, b) result (c)

      ! To perform XBRANCH data type substraction

      use sdata, only: xbranch


      type(xbranch), intent(in) :: a, b
      type(xbranch) :: c

      c%sigtr = a%sigtr - b%sigtr
      c%siga = a%siga - b%siga
      c%nuf = a%nuf - b%nuf
      c%sigf = a%sigf - b%sigf
      c%sigs = a%sigs - b%sigs
      c%dc = a%dc - b%dc

   end function brsubst

   !******************************************************************************!

   function brrealmult(re, a) result (b)

      ! To perform XBRANCH data type substraction

      use sdata, only: xbranch, dp


      real(dp), intent(in) :: re
      type(xbranch), intent(in) :: a
      type(xbranch) :: b

      b%sigtr = re * a%sigtr
      b%siga = re * a%siga
      b%nuf = re * a%nuf
      b%sigf = re * a%sigf
      b%sigs = re * a%sigs
      b%dc = re * a%dc

   end function brrealmult

   !******************************************************************************!

end module
