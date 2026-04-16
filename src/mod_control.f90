module control

   use sdata, only: dp
   implicit none
   save

   contains

   !******************************************************************************!

   subroutine forward()

      !
      ! Purpose:
      !    To solve forward (normal) problems
      !

      use sdata, only: nnod, aprad, apaxi, afrad, ftem, mtem, cden, &
      bcon, bpos, npow, th_niter, npow
      use io, only: asmpow, axipow, asmflux, inp_read, bther, &
      boutp, print_outp, bvtk, print_vtk
      use xsec, only: xs_updt
      use cmfd, only: outer, print_keff
      use th, only: th_iter


      !Update xsec
      call xs_updt(bcon, ftem, mtem, cden, bpos)

      call print_head()

      !Outer iteration
      if (bther == 0) then
         call outer(1)
      else
         allocate(npow(nnod))
         call th_iter(th_niter, 1)
         call print_tail()
      end if

      call print_keff()

      if (aprad == 1 .or. apaxi == 1) then
         if (.not. allocated(npow)) allocate(npow(nnod))
         call get_power_dist(npow)
      end if

      if (aprad == 1) call asmpow(npow)

      if (apaxi == 1) call axipow(npow)

      if (afrad == 1) call asmflux(1.e0_dp)

      if (boutp == 1) call print_outp(npow)

      if (bvtk == 1) call print_vtk(0)

   end subroutine forward

   !******************************************************************************!

   subroutine adjoint()

      !
      ! Purpose:
      !    To solve adjoint problems
      !

      use sdata, only: nnod, aprad, apaxi, afrad, ftem, mtem, &
      cden, bcon, bpos, npow
      use io, only: asmpow, axipow, asmflux, inp_read, &
      bvtk, print_vtk
      use xsec, only: xs_updt
      use cmfd, only: outer_ad, print_keff


      !Update xsec
      call xs_updt(bcon, ftem, mtem, cden, bpos)

      call print_head()

      !Outer iteration
      call outer_ad(1)

      call print_keff()

      if (aprad == 1 .or. apaxi == 1) then
         allocate(npow(nnod))
         call get_power_dist(npow)
      end if

      if (aprad == 1) call asmpow(npow)

      if (apaxi == 1) call axipow(npow)

      if (afrad == 1) call asmflux(1.e0_dp)

      if (bvtk == 1) call print_vtk(0)

   end subroutine adjoint

   !******************************************************************************!

   subroutine fixedsrc()

      !
      ! Purpose:
      !    To solve fixed source problems
      !

      use sdata, only: nnod, aprad, apaxi, afrad, ftem, mtem, cden, &
      bcon, bpos, powtot, npow
      use io, only: asmpow, axipow, asmflux, inp_read, bvtk, print_vtk
      use xsec, only: xs_updt
      use cmfd, only: outer_fs


      !Update xsec
      call xs_updt(bcon, ftem, mtem, cden, bpos)

      call print_head()

      !Outer iteration
      call outer_fs(1)

      if (aprad == 1 .or. apaxi == 1) then
         allocate(npow(nnod))
         call get_power_dist(npow)
      end if

      if (powtot > 0.0) then
         if (aprad == 1) call asmpow(npow)
         if (apaxi == 1) call axipow(npow)
      end if

      if (afrad == 1) call asmflux()

      if (bvtk == 1) call print_vtk(0)

   end subroutine fixedsrc

   !****************************************************************************!

   subroutine cbsearch()

      !
      ! Purpose:
      !    To search critical boron concentration
      !

      use sdata, only: ke, rbcon, ftem, mtem, cden, bpos, nnod, fer, ser, &
      aprad, apaxi, afrad, npow
      use io, only: ounit, asmflux, asmpow, axipow, bvtk, print_vtk, &
      print_outp, boutp
      use cmfd, only: outer
      use xsec, only: xs_updt


      real(dp) :: bc1, bc2, bcon     ! Boron Concentration
      real(dp) :: ke1, ke2
      integer :: n

      call print_head()

      bcon = rbcon
      call xs_updt(bcon, ftem, mtem, cden, bpos)
      call outer(0)
      bc1 = bcon
      ke1 = ke

      write(ounit, 1791) 1, bc1, ke1, ser, fer
      write(*, 1791) 1, bc1, ke1, ser, fer

      bcon = bcon + (ke - 1.) * bcon   ! Guess next critical boron concentration
      call xs_updt(bcon, ftem, mtem, cden, bpos)
      call outer(0)
      bc2 = bcon
      ke2 = ke

      write(ounit, 1791) 2, bc2, ke2, ser, fer
      write(*, 1791) 2, bc2, ke2, ser, fer

      n = 3
      do
         bcon = bc2 + (1._dp - ke2) / (ke1 - ke2) * (bc1 - bc2)
         call xs_updt(bcon, ftem, mtem, cden, bpos)
         call outer(0)
         bc1 = bc2
         bc2 = bcon
         ke1 = ke2
         ke2 = ke
         write(ounit, 1791) n, bcon, ke, ser, fer
         write(*, 1791) n, bcon, ke, ser, fer
         if ((abs(ke - 1._dp) < 1.e-5_dp) .and. (ser < 1.e-5_dp) .and. (fer < 1.e-5_dp)) exit
         n = n + 1
         call check_ppm(n, bcon)
      end do

      allocate(npow(nnod))
      if (aprad == 1 .or. apaxi == 1) then
         call get_power_dist(npow)
      end if

      if (aprad == 1) call asmpow(npow)

      if (apaxi == 1) call axipow(npow)

      if (afrad == 1) call asmflux(1._dp)

      if (bvtk == 1) call print_vtk(0)

      if (boutp == 1) call print_outp(npow)

      1791 format(i3, f10.2, f14.5, es14.5, es13.5)

   end subroutine cbsearch

   !****************************************************************************!

   subroutine cbsearcht()

      !
      ! Purpose:
      !    To search critical boron concentration with thermal feedback
      !

      use sdata, only: ke, bcon, rbcon, npow, nnod, &
      ser, fer, aprad, apaxi, afrad, npow, th_err, &
      serc, ferc
      use io, only: ounit, asmflux, asmpow, axipow, bvtk, print_vtk, &
      print_outp, boutp
      use cmfd, only: outer
      use th, only : th_iter
      ! use trans, only: vtk_out


      real(dp) :: bc1, bc2    ! Boron Concentration
      real(dp) :: ke1, ke2
      integer :: n
      character(len = 20) :: steady_name

      call print_head()

      allocate(npow(nnod))

      bcon = rbcon
      call th_iter(2, 0)  ! Start thermal hydarulic iteration with current paramters
      bc1 = bcon
      ke1 = ke

      write(ounit, 1792) 1, bc1, ke1, ser, fer, th_err
      write(*, 1792) 1, bc1, ke1, ser, fer, th_err

      if (bcon < 1.e-5) then
         bcon = 500.
      else
         bcon = bcon + (ke - 1.) * bcon   ! Guess next critical boron concentration
      end if
      call th_iter(2, 0)                 ! Perform second thermal hydarulic iteration with updated parameters
      bc2 = bcon
      ke2 = ke

      write(ounit, 1792) 2, bc2, ke2, ser, fer, th_err
      write(*, 1792) 2, bc2, ke2, ser, fer, th_err

      n = 3
      do
         bcon = bc2 + (1._dp - ke2) / (ke1 - ke2) * (bc1 - bc2)
         call th_iter(2, 0)
         bc1 = bc2
         bc2 = bcon
         ke1 = ke2
         ke2 = ke
         write(ounit, 1792) n, bcon, ke, ser, fer, th_err
         write(*, 1792) n, bcon, ke, ser, fer, th_err
         if ((abs(ke - 1._dp) < 1.e-5_dp) .and. (ser < serc) .and. (fer < ferc)) exit
         n = n + 1
         call check_ppm(n, bcon)
      end do

      if (aprad == 1 .or. apaxi == 1) then
         call get_power_dist(npow)
      end if

      if (aprad == 1) call asmpow(npow)

      if (apaxi == 1) call axipow(npow)

      if (afrad == 1) call asmflux(1._dp)

      if (bvtk == 1) call print_vtk(0)

      if (boutp == 1) call print_outp(npow)

      call print_tail()

      1792 format(i3, f9.2, f14.5, es14.5, es13.5, es17.5)

   end subroutine cbsearcht

   !****************************************************************************!

   subroutine check_ppm(n, bcon)

      !
      ! Purpose:
      !    To check critical boron concentration search
      !

      use io, only: ounit


      integer, intent(in) :: n
      real(dp), intent(in) :: bcon

      if (bcon > 3000.) then
         write(ounit, *) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
         write(ounit, *) '  KOMODO IS STOPPING'
         write(*, *) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
         stop
      end if
      if (bcon < 0.) then
         write(ounit, *) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
         write(ounit, *) '  KOMODO IS STOPPING'
         write(*, *) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
         stop
      end if
      if (n == 30) then
         write(ounit, *) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
         write(ounit, *) '  KOMODO IS STOPPING'
         write(*, *) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
         stop
      end if

   end subroutine check_ppm

   !****************************************************************************!

   subroutine get_power_dist (p)

      !
      ! Purpose:
      !    To calculate power for each nodes
      !


      use sdata, only: ng, nnod, sigf, f0, vdel, powtot, mode
      use io, only: ounit


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

      ! Normalize to 1._dp
      powtot = 0._dp
      do n = 1, nnod
         powtot = powtot + p(n)
      end do

      if (powtot <= 0 .and. mode /= 'FIXEDSRC') then
         write(ounit, *) '   ERROR: TOTAL NODES POWER IS ZERO OR LESS'
         write(ounit, *) '   stop IN subroutine get_power_dist'
         stop
      end if

      if (powtot > 0.0) then
         do n = 1, nnod
            p(n) = p(n) / powtot
         end do
      end if

   end subroutine

   !****************************************************************************!

   subroutine print_head()

      !
      ! Purpose:
      !    To print header
      !

      use sdata, only: mode
      use io, only: ounit, scr, bther

      if (mode == 'FORWARD' .and. bther == 1) then
         write(ounit, *) ;
         write(ounit, *)
         write(ounit, 3245) ;
         write(ounit, 3247) mode;
         write(ounit, 3245)
         write(ounit, *)
         write(ounit, 3251) ;
         write(ounit, 1179)
         if (scr) then
            write(*, *) ;
            write(*, *)
            write(*, 3245) ;
            write(*, 3247) mode;
            write(*, 3245)
            write(*, *)
            write(*, 3251) ;
            write(*, 1179)
         end if
      else if (mode == 'FORWARD' .or. mode == 'ADJOINT') then
         write(ounit, *) ;
         write(ounit, *)
         write(ounit, 3245) ;
         write(ounit, 3247) mode;
         write(ounit, 3245)
         write(ounit, *)
         write(ounit, 3248) ;
         write(ounit, 3249)
         if (scr) then
            write(*, *) ;
            write(*, *)
            write(*, 3245) ;
            write(*, 3247) mode;
            write(*, 3245)
            write(*, *)
            write(*, 3248) ;
            write(*, 3249)
         end if
      else if (mode == 'FIXEDSRC') then
         write(ounit, *) ;
         write(ounit, *)
         write(ounit, 3245) ;
         write(ounit, 3247) mode;
         write(ounit, 3245)
         write(ounit, *)
         write(ounit, 3248) ;
         write(ounit, 3249)
         if (scr) then
            write(*, *) ;
            write(*, *)
            write(*, 3245) ;
            write(*, 3247) mode;
            write(*, 3245)
            write(*, *)
            write(*, 3250) ;
            write(*, 3249)
         end if
      else
         if (bther == 0) then
            ! File Output
            write(ounit, *) ;
            write(ounit, *)
            write(ounit, 2176) ;
            write(ounit, 2177) ;
            write(ounit, 2176)
            write(ounit, *) ;
            write(ounit, 2178) ;
            write(ounit, 2179)
            if (scr) then
               ! Terminal Output
               write(*, *) ;
               write(*, *)
               write(*, 2176) ;
               write(*, 2177) ;
               write(*, 2176)
               write(*, *) ;
               write(*, 2178) ;
               write(*, 2179)
            end if
         else
            ! File Output
            write(ounit, *) ;
            write(ounit, *)
            write(ounit, 1176) ;
            write(ounit, 1177) ;
            write(ounit, 1176)
            write(ounit, *) ;
            write(ounit, 1178) ;
            write(ounit, 1179)
            if (scr) then
               ! Terminal Output
               write(*, *) ;
               write(*, *)
               write(*, 1176) ;
               write(*, 1177) ;
               write(*, 1176)
               write(*, *) ;
               write(*, 1178) ;
               write(*, 1179)
            end if
         end if
      end if



      3245 format (1x, ' ==============================================', &
      '================================')
      3247 format(23x, a8, ' CALCULATION RESULTS')
      3248 format(2x, 'Itr     k-eff     Fis.Src Error     Flux error')
      3249 format(1x, '----------------------------------------------------')
      3250 format(2x, 'Itr   Fis.Src Error     Flux error')
      2176 format(' ============================================================')
      2177 format(12x, 'CRITICAL BORON CONCENTRATION SEARCH')
      2178 format('Itr  Boron Conc.   K-eff     Flux Error   Fiss. Src. Error')
      2179 format(' -----------------------------------------------------------')
      1176 format &
      (' =========================================================================')
      1177 format(19x, 'CRITICAL BORON CONCENTRATION SEARCH')
      1178 format &
      ('Itr  Boron Conc.   K-eff     Flux Err.    Fiss. Src. Err.  Fuel Temp. Error.')
      1179 format &
      (' -----------------------------------------------------------------------')
      3251 format(2x, 'Itr     k-eff     Fis.Src Error     Flux error    Fuel Temp. Error')


   end subroutine print_head

   !****************************************************************************!

   subroutine print_tail()

      !
      ! Purpose:
      !    To print final th paramters
      !

      use sdata, only: ftem, mtem, cden, tfm
      use io, only: ounit, scr
      use th, only : par_ave, par_max, par_ave_out


      real(dp) :: tf, tm, mtm, mtf, otm, cd, ocd

      call par_ave(ftem, tf)
      call par_ave(mtem, tm)

      call par_max(tfm(:, 1), mtf)
      call par_max(mtem, mtm)

      call par_ave_out(mtem, otm)
      call par_ave(cden, cd)
      call par_ave_out(cden, ocd)

      ! Write Output
      write(ounit, *)
      write(ounit, 5001) tf, tf - 273.15;
      write(ounit, 5002) mtf, mtf - 273.15
      write(ounit, 5003) tm, tm - 273.15;
      write(ounit, 5004) mtm, mtm - 273.15
      write(ounit, 5005) otm, otm - 273.15;
      write(ounit, 5006) cd * 1000. , cd
      write(ounit, 5007) ocd * 1000. , ocd
      if (scr) then
         write(*, *)
         write(*, 5001) tf, tf - 273.15;
         write(*, 5002) mtf, mtf - 273.15
         write(*, 5003) tm, tm - 273.15;
         write(*, 5004) mtm, mtm - 273.15
         write(*, 5005) otm, otm - 273.15;
         write(*, 5006) cd * 1000. , cd
         write(*, 5007) ocd * 1000. , ocd
      end if

      5001 format(2x, 'AVERAGE FUEL TEMPERATURE        : ', f7.1, ' K (', f7.1, ' C)')
      5002 format(2x, 'MAX FUEL CENTERLINE TEMPERATURE : ', f7.1, ' K (', f7.1, ' C)')
      5003 format(2x, 'AVERAGE MODERATOR TEMPERATURE   : ', f7.1, ' K (', f7.1, ' C)')
      5004 format(2x, 'MAXIMUM MODERATOR TEMPERATURE   : ', f7.1, ' K (', f7.1, ' C)')
      5005 format(2x, 'OUTLET MODERATOR TEMPERATURE    : ', f7.1, ' K (', f7.1, ' C)')
      5006 format(2x, 'AVERAGE MODERATOR DENSITY       : ', f7.1, ' kg/m3 (', f7.3, ' g/cc)')
      5007 format(2x, 'OUTLET MODERATOR DENSITY        : ', f7.1, ' kg/m3 (', f7.3, ' g/cc)')

   end subroutine print_tail



end module
