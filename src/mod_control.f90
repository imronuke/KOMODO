module control

  use sdata, only: dp
  implicit none
  save

contains

  !******************************************************************************!

  SUBROUTINE forward()

    !
    ! Purpose:
    !    To solve forward (normal) problems
    !

    use sdata, only: nnod, aprad, apaxi, afrad, ftem, mtem, cden, &
    bcon, bpos, npow, th_niter
    use io,    only: AsmPow, AxiPow, AsmFlux, inp_read, bther, boutp, print_outp
    use xsec,  only: XS_updt
    use cmfd,  only: outer,print_keff
    use th,  only: th_iter

    IMPLICIT NONE

    REAL(DP), DIMENSION(:), ALLOCATABLE :: pow

    !Update xsec
    CALL XS_updt(bcon, ftem, mtem, cden, bpos)

    call print_head()

    !Outer iteration
    if (bther == 0) then
      CALL outer(1)
    else
      allocate(npow(nnod))
      call th_iter(th_niter, 1)
      call print_tail()
    end if

    call print_keff()

    IF (aprad == 1 .OR. apaxi == 1) THEN
        ALLOCATE(pow(nnod))
        CALL get_power_dist(pow)
    END IF

    IF (aprad == 1) CALL AsmPow(pow)

    IF (apaxi == 1) CALL AxiPow(pow)

    IF (afrad == 1) CALL AsmFlux(1.e0_DP)

    IF (boutp == 1) CALL print_outp(pow)


  END SUBROUTINE forward

  !******************************************************************************!

  SUBROUTINE adjoint()

    !
    ! Purpose:
    !    To solve adjoint problems
    !

    use sdata, only: nnod, aprad, apaxi, afrad, ftem, mtem, cden, bcon, bpos
    use io,    only: AsmPow, AxiPow, AsmFlux, inp_read
    use xsec,  only: XS_updt
    use cmfd,  only: outer_ad,print_keff

    IMPLICIT NONE

    REAL(DP), DIMENSION(:), ALLOCATABLE :: pow

    !Update xsec
    CALL XS_updt(bcon, ftem, mtem, cden, bpos)

    call print_head()

    !Outer iteration
    CALL outer_ad(1)

    call print_keff()

    IF (aprad == 1 .OR. apaxi == 1) THEN
        ALLOCATE(pow(nnod))
        CALL get_power_dist(pow)
    END IF

    IF (aprad == 1) CALL AsmPow(pow)

    IF (apaxi == 1) CALL AxiPow(pow)

    IF (afrad == 1) CALL AsmFlux(1.e0_DP)

  END SUBROUTINE adjoint

  !******************************************************************************!

  SUBROUTINE fixedsrc()

    !
    ! Purpose:
    !    To solve fixed source problems
    !

    use sdata, only: nnod, aprad, apaxi, afrad, ftem, mtem, cden, bcon, bpos, powtot
    use io,    only: AsmPow, AxiPow, AsmFlux, inp_read
    use xsec,  only: XS_updt
    use cmfd,  only: outer_fs

    IMPLICIT NONE

    REAL(DP), DIMENSION(:), ALLOCATABLE :: pow

    !Update xsec
    CALL XS_updt(bcon, ftem, mtem, cden, bpos)

    call print_head()

    !Outer iteration
    CALL outer_fs(1)

    IF (aprad == 1 .OR. apaxi == 1) THEN
        ALLOCATE(pow(nnod))
        CALL get_power_dist(pow)
    END IF

    IF (powtot > 0.0) THEN
      IF (aprad == 1) CALL AsmPow(pow)
      IF (apaxi == 1) CALL AxiPow(pow)
    END IF

    IF (afrad == 1) CALL AsmFlux()

  END SUBROUTINE fixedsrc

 !****************************************************************************!

  SUBROUTINE cbsearch()

  !
  ! Purpose:
  !    To search critical boron concentration
  !

  USE sdata, ONLY: Ke, rbcon, ftem, mtem, cden, bpos, nnod, fer, ser, &
                   aprad, apaxi, afrad, npow
  USE io, ONLY: ounit, AsmFlux, AsmPow, AxiPow
  USE cmfd, ONLY: outer
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
      call check_ppm(n, bcon)
  END DO

  ALLOCATE(npow(nnod))
  IF (aprad == 1 .OR. apaxi == 1) THEN
      CALL get_power_dist(npow)
  END IF

  IF (aprad == 1) CALL AsmPow(npow)

  IF (apaxi == 1) CALL AxiPow(npow)

  IF (afrad == 1) CALL AsmFlux(1._DP)

  1791 format(I3, F10.2, F14.5, ES14.5, ES13.5)

  END SUBROUTINE cbsearch

  !****************************************************************************!

  SUBROUTINE cbsearcht()

  !
  ! Purpose:
  !    To search critical boron concentration with thermal feedback
  !

  USE sdata, ONLY: Ke, bcon, rbcon, npow, nnod, &
                   ser, fer, aprad, apaxi, afrad, npow, th_err, &
                   serc, ferc
  USE io, ONLY: ounit, AsmFlux, AsmPow, AxiPow
  USE cmfd, ONLY: outer
  use th, only : th_iter

  IMPLICIT NONE

  REAL(DP)  :: bc1, bc2    ! Boron Concentration
  REAL(DP) :: ke1, ke2
  INTEGER :: n

  call print_head()

  ALLOCATE(npow(nnod))

  bcon = rbcon
  CALL th_iter(2, 0)  ! Start thermal hydarulic iteration with current paramters
  bc1 = bcon
  ke1 = Ke

  WRITE(ounit,1792) 1, bc1, Ke1, ser, fer, th_err
  WRITE(*,1792) 1, bc1, Ke1, ser, fer, th_err

  IF (bcon < 1.e-5) THEN
    bcon = 500.
  ELSE
    bcon = bcon + (Ke - 1.) * bcon   ! Guess next critical boron concentration
  END IF
  CALL th_iter(2, 0)                 ! Perform second thermal hydarulic iteration with updated parameters
  bc2 = bcon
  ke2 = Ke

  WRITE(ounit,1792) 2, bc2, Ke2, ser, fer, th_err
  WRITE(*,1792) 2, bc2, Ke2, ser, fer, th_err

  n = 3
  DO
      bcon = bc2 + (1._DP - ke2) / (ke1 - ke2) * (bc1 - bc2)
      CALL th_iter(2, 0)
      bc1 = bc2
      bc2 = bcon
      ke1 = ke2
      ke2 = ke
      WRITE(ounit,1792) n, bcon, Ke, ser, fer, th_err
      WRITE(*,1792) n, bcon, Ke, ser, fer, th_err
      IF ((ABS(Ke - 1._DP) < 1.e-5_DP) .AND. (ser < serc) .AND. (fer < ferc)) EXIT
      n = n + 1
      call check_ppm(n, bcon)
  END DO

  IF (aprad == 1 .OR. apaxi == 1) THEN
      CALL get_power_dist(npow)
  END IF

  IF (aprad == 1) CALL AsmPow(npow)

  IF (apaxi == 1) CALL AxiPow(npow)

  IF (afrad == 1) CALL AsmFlux(1._DP)

  call print_tail()

  1792 format(I3, F9.2, F14.5, ES14.5, ES13.5, ES17.5)

  END SUBROUTINE cbsearcht

  !****************************************************************************!

  SUBROUTINE check_ppm(n, bcon)

  !
  ! Purpose:
  !    To check critical boron concentration search
  !

  USE io, ONLY: ounit

  IMPLICIT NONE

  integer, intent(in) :: n
  real(dp), intent(in) :: bcon

  IF (bcon > 3000.) THEN
      WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
      WRITE(ounit,*) '  KOMODO IS STOPPING'
      WRITE(*,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
      STOP
  END IF
  IF (bcon < 0.) THEN
      WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
      WRITE(ounit,*) '  KOMODO IS STOPPING'
      WRITE(*,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
      STOP
  END IF
  IF (n == 30) THEN
      WRITE(ounit,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
      WRITE(ounit,*) '  KOMODO IS STOPPING'
      WRITE(*,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
      STOP
  END IF

END SUBROUTINE check_ppm

!****************************************************************************!

subroutine get_power_dist (p)

  !
  ! Purpose:
  !    To calculate power for each nodes
  !
  
  
  USE sdata, ONLY: ng, nnod, sigf, f0, vdel, powtot, mode
  USE io,    ONLY: ounit
  
  implicit none
  
  real(dp), dimension(:), intent(out) :: p
  integer :: g, n
  real(dp) :: pow, vtot
  
  p = 0._dp
  do g= 1, ng
      do n= 1, nnod
        pow = f0(n,g) * sigf(n,g)
        if (pow < 0.) pow = 0.
        p(n) = p(n) + pow
      end do
  end do
  
  ! Normalize to 1._DP
  powtot = 0._DP
  vtot = 0.0
  do n = 1, nnod
      powtot = powtot + p(n) * vdel(n)
      vtot   = vtot + vdel(n)
  end do
  
  if (powtot <= 0 .AND. mode /= 'FIXEDSRC') THEN
     write(ounit, *) '   ERROR: TOTAL NODES POWER IS ZERO OR LESS'
     write(ounit, *) '   STOP IN subroutine get_power_dist'
     STOP
  end if
  
  if (powtot > 0.0) then
     do n = 1, nnod
         p(n) = p(n) * vtot / powtot
     end do
  end if
  
  end subroutine

  !****************************************************************************!

  SUBROUTINE print_head()

    !
    ! Purpose:
    !    To print header
    !

    use sdata, only: mode
    use io,    only: ounit, scr, bther

    IMPLICIT NONE
    if (mode == 'FORWARD' .and. bther == 1) then
      WRITE(ounit,*); WRITE(ounit,*)
      WRITE(ounit,3245); WRITE(ounit,3247) mode; WRITE(ounit,3245)
      WRITE(ounit,*)
      WRITE(ounit,3251); WRITE(ounit,1179)
      if (scr) then
        WRITE(*,*); WRITE(*,*)
        WRITE(*,3245); WRITE(*,3247) mode; WRITE(*,3245)
        WRITE(*,*)
        WRITE(*,3251); WRITE(*,1179)
      end if
    else if (mode == 'FORWARD' .or. mode == 'ADJOINT') then
      WRITE(ounit,*); WRITE(ounit,*)
      WRITE(ounit,3245); WRITE(ounit,3247) mode; WRITE(ounit,3245)
      WRITE(ounit,*)
      WRITE(ounit,3248); WRITE(ounit,3249)
      if (scr) then
        WRITE(*,*); WRITE(*,*)
        WRITE(*,3245); WRITE(*,3247) mode; WRITE(*,3245)
        WRITE(*,*)
        WRITE(*,3248); WRITE(*,3249)
      end if
    else if (mode == 'FIXEDSRC') then
      WRITE(ounit,*); WRITE(ounit,*)
      WRITE(ounit,3245); WRITE(ounit,3247) mode; WRITE(ounit,3245)
      WRITE(ounit,*)
      WRITE(ounit,3248); WRITE(ounit,3249)
      if (scr) then
        WRITE(*,*); WRITE(*,*)
        WRITE(*,3245); WRITE(*,3247) mode; WRITE(*,3245)
        WRITE(*,*)
        WRITE(*,3250); WRITE(*,3249)
      end if
    else
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
    end if



    3245 format (1X, ' ==============================================', &
    '================================')
    3247 format(23X, A8, ' CALCULATION RESULTS')
    3248 format(2X,'Itr     k-eff     Fis.Src Error     Flux error')
    3249 format(1X,'----------------------------------------------------')
    3250 format(2X,'Itr   Fis.Src Error     Flux error')
    2176 format(' ============================================================')
    2177 format(12X,'CRITICAL BORON CONCENTRATION SEARCH')
    2178 format('Itr  Boron Conc.   K-eff     Flux Error   Fiss. Src. Error')
    2179 format(' -----------------------------------------------------------')
    1176 format  &
    (' =========================================================================')
    1177 format(19X,'CRITICAL BORON CONCENTRATION SEARCH')
    1178 format &
    ('Itr  Boron Conc.   K-eff     Flux Err.    Fiss. Src. Err.  Fuel Temp. Error.')
    1179 format &
    (' -----------------------------------------------------------------------')
    3251 format(2X,'Itr     k-eff     Fis.Src Error     Flux error    Fuel Temp. Error')


  END SUBROUTINE print_head

  !****************************************************************************!

  SUBROUTINE print_tail()

    !
    ! Purpose:
    !    To print final th paramters
    !

    use sdata, only: ftem, mtem, cden, tfm
    use io, only: ounit, scr
    use th, only : par_ave, par_max, par_ave_out

    IMPLICIT NONE

    REAL(DP) :: tf, tm, mtm, mtf, otm, cd, ocd

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

    5001 FORMAT(2X, 'AVERAGE FUEL TEMPERATURE        : ', F7.1, ' K (', F7.1, ' C)')
    5002 FORMAT(2X, 'MAX FUEL CENTERLINE TEMPERATURE : ', F7.1, ' K (', F7.1, ' C)')
    5003 FORMAT(2X, 'AVERAGE MODERATOR TEMPERATURE   : ', F7.1, ' K (', F7.1, ' C)')
    5004 FORMAT(2X, 'MAXIMUM MODERATOR TEMPERATURE   : ', F7.1, ' K (', F7.1, ' C)')
    5005 FORMAT(2X, 'OUTLET MODERATOR TEMPERATURE    : ', F7.1, ' K (', F7.1, ' C)')
    5006 FORMAT(2X, 'AVERAGE MODERATOR DENSITY       : ', F7.1, ' kg/m3 (', F7.3, ' g/cc)')
    5007 FORMAT(2X, 'OUTLET MODERATOR DENSITY        : ', F7.1, ' kg/m3 (', F7.3, ' g/cc)')

  END SUBROUTINE print_tail



end module
