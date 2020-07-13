module xsec

  use sdata, only: dp
  implicit none
  save

contains

  !******************************************************************************!

  SUBROUTINE XS_updt (xbcon, xftem, xmtem, xcden, xbpos)
  !
  ! Purpose:
  !    To update XS for given TH paramaters and rod position
  !

  use io, only: bbcon, bftem, bmtem, bcden, bcrod, bcbcs
  use sdata, only: get_time, xs_time


  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: xbcon  ! Provided Boron Concentration
  REAL(DP), DIMENSION(:), INTENT(IN) :: xftem  ! Provided fuel temperature
  REAL(DP), DIMENSION(:), INTENT(IN) :: xmtem  ! Provided moderator temperature
  REAL(DP), DIMENSION(:), INTENT(IN) :: xcden  ! Provided coolant density
  REAL(DP), DIMENSION(:), INTENT(IN) :: xbpos  ! Provided control rod bank position

  REAL(DP) :: st, fn

  st = get_time()

  CALL base_updt()
  IF (bbcon == 1 .OR. bcbcs == 1) CALL bcon_updt(xbcon)
  IF (bftem == 1) CALL ftem_updt(xftem)
  IF (bmtem == 1) CALL mtem_updt(xmtem)
  IF (bcden == 1) CALL cden_updt(xcden)
  IF (bcrod == 1) CALL crod_updt(xbpos)
  CALL Dsigr_updt()

  call check_xs()

  fn = get_time()
  xs_time = xs_time + (fn-st)

  END SUBROUTINE XS_updt

  !******************************************************************************!

  SUBROUTINE XStab_updt (xbcon, xftem, xmtem, xcden, xbpos)
  !
  ! Purpose:
  !    To update XS for given TH paramaters and rod position (when xtab file exist)
  !

  USE sdata, ONLY : nnod, mat, sigtr, siga, nuf, sigf, sigs, dc
  use io,    only : bcrod
  use sdata, only: get_time, xs_time

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: xbcon  ! Provided Boron Concentration
  REAL(DP), DIMENSION(:), INTENT(IN) :: xftem  ! Provided fuel temperature
  REAL(DP), DIMENSION(:), INTENT(IN) :: xmtem  ! Provided moderator temperature
  REAL(DP), DIMENSION(:), INTENT(IN) :: xcden  ! Provided coolant density
  REAL(DP), DIMENSION(:), INTENT(IN) :: xbpos  ! Provided control rod bank position

  INTEGER :: n
  REAL(DP) :: st, fn

  st = get_time()

  DO n = 1, nnod
    CALL brInterp(0, mat(n), xcden(n),  xbcon, xftem(n), xmtem(n), sigtr(n,:), &
    siga(n,:), nuf(n,:), sigf(n,:), sigs(n,:,:), dc(n,:,:))
  END DO
  IF (bcrod == 1) CALL crod_tab_updt(xbpos)
  CALL Dsigr_updt()

  call check_xs()

  fn = get_time()
  xs_time = xs_time + (fn-st)


  END SUBROUTINE XStab_updt

  !******************************************************************************!

  SUBROUTINE check_xs ()
  !
  ! Purpose:
  !    To check errors in xsec
  !

  use sdata, only : D, sigr, nuf, chi, sigs, nnod, ng, mat, nmat
  use io,    only : ounit


  IMPLICIT NONE

  INTEGER :: n, g, h

  DO g = 1, ng
    DO n = 1, nnod
      IF (D(n,g) < 1.e-20) then
        WRITE(*,*)
        WRITE(*,3541) mat(n), g
        WRITE(*,*) ' DIFFUSION COEF. IS CLOSE TO ZERO OR NEGATIVE'
        WRITE(ounit,3541) mat(n), g
        WRITE(ounit,*) ' DIFFUSION COEF. IS CLOSE TO ZERO OR NEGATIVE'
        STOP
      END IF
      IF (sigr(n,g) < 0.) then
        WRITE(*,*)
        WRITE(*,3542) mat(n), g
        WRITE(*,*) ' REMOVAL XS IS  NEGATIVE'
        WRITE(ounit,3542) mat(n), g
        WRITE(ounit,*) ' REMOVAL XS IS NEGATIVE'
        STOP
      END IF
      IF (nuf(n,g) < 0.) then
        WRITE(*,*)
        WRITE(*,3544) mat(n), g
        WRITE(*,*) ' NU*FISSION XS IS  NEGATIVE'
        WRITE(ounit,3544) mat(n), g
        WRITE(ounit,*) ' NU*FISSION XS IS NEGATIVE'
        STOP
      END IF
    END DO
    DO n = 1, nmat
      IF (chi(n,g) < 0.) then
        WRITE(*,*)
        WRITE(*,3543) n, g
        WRITE(*,*) ' FISSION SPECTRUM IS  NEGATIVE'
        WRITE(ounit,3543) n, g
        WRITE(ounit,*) ' FISSION SPECTRUM IS NEGATIVE'
        STOP
      END IF
    END DO
    DO h = 1, ng
      DO n = 1, nnod
        IF (sigs(n,g,h) < 0.) then
          WRITE(*,*)
          WRITE(*,3545) mat(n), g, h
          WRITE(*,*) ' SCATTERING XS IS  NEGATIVE'
          WRITE(ounit,3545) mat(n), g, h
          WRITE(ounit,*) ' SCATTERING XS IS NEGATIVE'
          STOP
        END IF
      END DO
    END DO
  END DO

  3541 FORMAT (2X, 'ERROR IN THE DIFFUSION COEFFICIENT ', &
  'FOR MATERIAL NUMBER : ', I3, ' ENERGY GROUP : ', I2)
  3542 FORMAT (2X, 'ERROR IN THE REMOVAL CROSS SECTION ', &
  'FOR MATERIAL NUMBER : ', I3, ' ENERGY GROUP : ', I2)
  3543 FORMAT (2X, 'ERROR IN THE FISSION SPECTRUM ', &
  'FOR MATERIAL NUMBER : ', I3, ' ENERGY GROUP : ', I2)
  3544 FORMAT (2X, 'ERROR IN THE NU*FISSION XSEC ', &
  'FOR MATERIAL NUMBER : ', I3, ' ENERGY GROUP : ', I2)
  3545 FORMAT (2X, 'ERROR IN THE SCATTERING CROSS SECTION ', &
  'FOR MATERIAL NUMBER : ', I3, ' ENERGY GROUP : ', I2, &
   ' TO ENERGY GROUP : ', I2)


  END SUBROUTINE check_xs

  !******************************************************************************!

  SUBROUTINE base_updt ()
  !
  ! Purpose:
  !    To update current XS to base XS
  !


  USE sdata, ONLY: nnod, sigtr, siga, nuf, sigf, sigs, mat, &
                   xsigtr, xsiga, xnuf, xsigf, xsigs

  IMPLICIT NONE

  INTEGER :: n

  DO n = 1, nnod
      sigtr(n,1:)   = xsigtr(mat(n),1:)
      siga (n,1:)   = xsiga (mat(n),1:)
      nuf  (n,1:)   = xnuf  (mat(n),1:)
      sigf (n,1:)   = xsigf (mat(n),1:)
      sigs (n,1:,1:) = xsigs (mat(n),1:,1:)
  END DO


  END SUBROUTINE base_updt

  !******************************************************************************!

  SUBROUTINE Dsigr_updt ()
  !
  ! Purpose:
  !    To update diffusion constant and removal XS
  !


  USE sdata, ONLY: nnod, ng, sigtr, siga,  &
                   sigs, D, sigr

  IMPLICIT NONE

  INTEGER :: i, g, h
  REAL(DP) :: dum

  DO i = 1, nnod
    DO g = 1, ng
      IF (sigtr(i,g) < 1.e-5 ) STOP "Negative diffusion coefficient encountered"
      D(i,g) = 1._dp/(3._dp*sigtr(i,g))
      dum = 0.
      DO h= 1, ng
          IF (g /= h) dum = dum + sigs(i,g,h)
      END DO
      sigr(i,g) =  siga(i,g) + dum
    END DO
  END DO

  END SUBROUTINE Dsigr_updt

  !******************************************************************************!

  SUBROUTINE crod_updt (bpos)
  !
  ! Purpose: TO UPDATE AND CALCUALTE VOLUME WEIGHTED HOMOGENIZED CX FOR RODDED NODE
  !

  USE sdata, ONLY: ng, nxx, nyy, nzz, xyz, zdel, mat, &
                   sigtr, siga, nuf, sigf, sigs, &
                   dsigtr, dsiga, dnuf, dsigf, dsigs, &
                   coreh, fbmap, pos0, ssize

  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(IN) :: bpos

  INTEGER ::i, j, k, g, h
  REAL(DP) :: rodh, vfrac
  REAL(DP) :: dum

  ! For each node
  DO j = 1, nyy
    DO i = 1, nxx
       IF (fbmap(i,j) > 0) THEN
          !!!(rodh -> posistion the tip of the control rod the top of core)
           rodh = coreh - pos0  - bpos(fbmap(i,j))*ssize
           dum = 0._dp
           DO k = nzz, 1, -1
             ! For partially rodded node, get volume weighted homogenized CX (0 < vfrac < 1._DP)
             IF (rodh >= dum .AND. rodh <= dum+zdel(k)) THEN   ! If this node partially rodded
               vfrac = (rodh - dum) / zdel(k)
               sigtr(xyz(i,j,k),:) = sigtr(xyz(i,j,k),:) + &
                                  vfrac * dsigtr(mat(xyz(i,j,k)),:)
               siga(xyz(i,j,k),:)  = siga(xyz(i,j,k),:) + &
                                  vfrac * dsiga(mat(xyz(i,j,k)),:)
               nuf(xyz(i,j,k),:)   = nuf(xyz(i,j,k),:) + &
                                  vfrac * dnuf(mat(xyz(i,j,k)),:)
               sigf(xyz(i,j,k),:)  = sigf(xyz(i,j,k),:) + &
                                  vfrac * dsigf(mat(xyz(i,j,k)),:)
               sigs(xyz(i,j,k),:,:)  = sigs(xyz(i,j,k),:,:) + &
                                    vfrac * dsigs(mat(xyz(i,j,k)),:,:)
               EXIT
             END IF
             ! For fully rodded node, vfrac = 1.
             sigtr(xyz(i,j,k),:) = sigtr(xyz(i,j,k),:) + dsigtr(mat(xyz(i,j,k)),:)
             siga(xyz(i,j,k),:)  = siga(xyz(i,j,k),:) + dsiga(mat(xyz(i,j,k)),:)
             nuf(xyz(i,j,k),:)   = nuf(xyz(i,j,k),:) + dnuf(mat(xyz(i,j,k)),:)
             sigf(xyz(i,j,k),:)  = sigf(xyz(i,j,k),:) + dsigf(mat(xyz(i,j,k)),:)
             sigs(xyz(i,j,k),:,:)  = sigs(xyz(i,j,k),:,:) + dsigs(mat(xyz(i,j,k)),:,:)

             dum = dum + zdel(k)
           END DO
           ! if negative CX found, Surpress CX to zero and calculate D and sigr
           DO k = nzz, 1, -1
             DO g = 1, ng
                IF (siga(xyz(i,j,k),g) < 0.) siga(xyz(i,j,k),g) = 0.
                IF (nuf(xyz(i,j,k),g) < 0.) nuf(xyz(i,j,k),g) = 0.
                IF (sigf(xyz(i,j,k),g) < 0.) sigf(xyz(i,j,k),g) = 0.
                DO h = 1, ng
                  IF (sigs(xyz(i,j,k),g,h) < 0.) sigs(xyz(i,j,k),g,h) = 0.
                END DO
             END DO
           END DO
       END IF
    END DO
  END DO


  END SUBROUTINE crod_updt

  !******************************************************************************!

  SUBROUTINE crod_tab_updt (bpos)
  !
  ! Purpose: TO UPDATE AND CALCUALTE VOLUME WEIGHTED HOMOGENIZED CX FOR RODDED NODE
  !

  USE sdata, ONLY: ng, nxx, nyy, nzz, xyz, zdel, mat, &
                   sigtr, siga, nuf, sigf, sigs, dc, &
                   coreh, fbmap, pos0, ssize, m, &
                   nnod, cden, ftem, mtem, bcon

  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(IN) :: bpos

  INTEGER ::i, j, k, g, h
  REAL(DP) :: rodh, vfrac
  REAL(DP) :: dum

  INTEGER :: nd

  REAL(DP), DIMENSION(nnod,ng) :: rsigtr, rsiga, rnuf, rsigf
  REAL(DP), DIMENSION(nnod,ng,ng) :: rsigs
  REAL(DP), DIMENSION(nnod,ng,6) :: rdc

  DO j = 1, nyy
    DO i = 1, nxx
      IF (fbmap(i,j) > 0) THEN
         !(rodh -> posistion the tip of the control rod the top of core)
         rodh = coreh - pos0  - bpos(fbmap(i,j))*ssize
         dum = 0._DP
           DO k = nzz, 1, -1
             nd = xyz(i,j,k)                                  ! Node number
             IF (m(mat(nd))%trod == 1) THEN
               CALL brInterp(1, mat(nd), cden(nd), bcon, ftem(nd), &
               mtem(nd), rsigtr(nd,:), rsiga(nd,:), rnuf(nd,:), &
               rsigf(nd,:), rsigs(nd,:,:), rdc(nd,:,:))
             ELSE
               WRITE(101,1671) fbmap(i,j), mat(nd)
               WRITE(*,1671) fbmap(i,j), mat(nd)
               STOP
               1671 FORMAT (2X, 'CONTROL ROD BANK NUMBER ', I4, &
               ' COINCIDES WITH MATERIAL NUMBER ', I4, &
               ' THAT DOES NOT HAVE CONTROL ROD DATA IN XTAB FILE')
             END IF
             ! For partially rodded node, get volume weighted homogenized CX (0 < vfrac < 1._DP)
             IF (rodh >= dum .AND. rodh <= dum+zdel(k)) THEN   ! If this node partially rodded
               vfrac = (rodh - dum) / zdel(k)                    ! Rodded fraction
               sigtr(nd,:)  = (1._DP - vfrac) * sigtr(nd,:)  + &
               vfrac * rsigtr(nd,:)
               siga(nd,:)   = (1._DP - vfrac) * siga(nd,:)   + &
               vfrac * rsiga(nd,:)
               nuf(nd,:)    = (1._DP - vfrac) * nuf(nd,:)    + &
               vfrac * rnuf(nd,:)
               sigf(nd,:)   = (1._DP - vfrac) * sigf(nd,:)   + &
               vfrac * rsigf(nd,:)
               sigs(nd,:,:) = (1._DP - vfrac) * sigs(nd,:,:) + &
               vfrac * rsigs(nd,:,:)
               dc(nd,:,:)   = (1._DP - vfrac) * dc(nd,:,:)   + &
               vfrac * rdc(nd,:,:)
               EXIT
             END IF
             ! For fully rodded node, vfrac = 1.
             sigtr(nd,:)  = rsigtr(nd,:)
             siga(nd,:)   = rsiga(nd,:)
             nuf(nd,:)    = rnuf(nd,:)
             sigf(nd,:)   = rsigf(nd,:)
             sigs(nd,:,:) = rsigs(nd,:,:)
             dc(nd,:,:)   = rdc(nd,:,:)

             dum = dum + zdel(k)
           END DO
           ! if negative CX found, Surpress CX to zero  and calculate D and sigr
           DO k = nzz, 1, -1
             nd = xyz(i,j,k)
             DO g = 1, ng
               IF (siga(nd,g) < 0.) siga(nd,g)  = 0.
               IF (nuf(nd,g)  < 0.) nuf(nd,g)   = 0.
               IF (sigf(nd,g) < 0.) sigf(nd,g)  = 0.
               DO h = 1, ng
                 IF (sigs(nd,g,h) < 0.) sigs( nd,g,h) = 0.
               END DO
               DO h = 1, 6
                 IF (dc(nd,g,h) < 0.) dc(nd,g,h) = 0.
               END DO
             END DO
           END DO
       END IF
    END DO
  END DO

  END SUBROUTINE crod_tab_updt

  !******************************************************************************!

  SUBROUTINE bcon_updt (bcon)

  !
  ! Purpose:
  !    To update CX for given boron concentration

  USE sdata, ONLY: nnod, ng, sigtr, siga, nuf, sigf, sigs, mat, &
                   csigtr, csiga, cnuf, csigf, csigs, rbcon

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: bcon
  INTEGER :: n, g, h

  DO n = 1, nnod
      DO g = 1, ng
          sigtr(n,g) = sigtr(n,g) + csigtr(mat(n),g) * (bcon - rbcon)
          siga(n,g)  = siga(n,g)  + csiga(mat(n),g)  * (bcon - rbcon)
          nuf(n,g)   = nuf(n,g)   + cnuf(mat(n),g)   * (bcon - rbcon)
          sigf(n,g)  = sigf(n,g)  + csigf(mat(n),g)  * (bcon - rbcon)
          DO h = 1, ng
              sigs(n,g,h) = sigs(n,g,h) + csigs(mat(n),g,h) * (bcon - rbcon)
          END DO
      END DO
  END DO


  END SUBROUTINE bcon_updt

  !******************************************************************************!

  SUBROUTINE ftem_updt (ftem)

  !
  ! Purpose:
  !    To update CX for given fuel temp

  USE sdata, ONLY: nnod, ng, sigtr, siga, nuf, sigf, sigs, mat, &
                   fsigtr, fsiga, fnuf, fsigf, fsigs, rftem

  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(IN) :: ftem
  INTEGER :: n, g, h


  DO n = 1, nnod
      DO g = 1, ng
          sigtr(n,g) = sigtr(n,g) + fsigtr(mat(n),g) * (SQRT(ftem(n))- SQRT(rftem))
          siga(n,g)  = siga(n,g)  + fsiga(mat(n),g)  * (SQRT(ftem(n)) - SQRT(rftem))
          nuf(n,g)   = nuf(n,g)   + fnuf(mat(n),g)   * (SQRT(ftem(n)) - SQRT(rftem))
          sigf(n,g)  = sigf(n,g)  + fsigf(mat(n),g)  * (SQRT(ftem(n)) - SQRT(rftem))
          DO h = 1, ng
             sigs(n,g,h) = sigs(n,g,h) + fsigs(mat(n),g,h) * (SQRT(ftem(n)) - SQRT(rftem))
          END DO
        END DO
  END DO

  END SUBROUTINE ftem_updt

  !******************************************************************************!

  SUBROUTINE mtem_updt (mtem)

  !
  ! Purpose:
  !    To update CX for given moderator temperature

  USE sdata, ONLY: nnod, ng, sigtr, siga, nuf, sigf, sigs, mat, &
                   msigtr, msiga, mnuf, msigf, msigs, rmtem

  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(IN) :: mtem
  INTEGER :: n, g, h


  DO n = 1, nnod
      DO g = 1, ng
          sigtr(n,g) = sigtr(n,g) + msigtr(mat(n),g) * (mtem(n) - rmtem)
          siga(n,g)  = siga(n,g)  + msiga(mat(n),g)  * (mtem(n) - rmtem)
          nuf(n,g)   = nuf(n,g)   + mnuf(mat(n),g)   * (mtem(n) - rmtem)
          sigf(n,g)  = sigf(n,g)  + msigf(mat(n),g)  * (mtem(n) - rmtem)
          DO h = 1, ng
              sigs(n,g,h) = sigs(n,g,h) + msigs(mat(n),g,h) * (mtem(n) - rmtem)
          END DO
      END DO
  END DO


  END SUBROUTINE mtem_updt

  !******************************************************************************!

  SUBROUTINE cden_updt (cden)

  !
  ! Purpose:
  !    To update CX for given coolant density

  USE sdata, ONLY: nnod, ng, sigtr, siga, nuf, sigf, sigs, mat, &
                   lsigtr, lsiga, lnuf, lsigf, lsigs, rcden

  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(IN) :: cden
  INTEGER :: n, g, h


  DO n = 1, nnod
      DO g = 1, ng
          sigtr(n,g) = sigtr(n,g) + lsigtr(mat(n),g) * (cden(n) - rcden)
          siga(n,g)  = siga(n,g)  + lsiga(mat(n),g)  * (cden(n) - rcden)
          nuf(n,g)   = nuf(n,g)   + lnuf(mat(n),g)   * (cden(n) - rcden)
          sigf(n,g)  = sigf(n,g)  + lsigf(mat(n),g)  * (cden(n) - rcden)
          DO h = 1, ng
              sigs(n,g,h) = sigs(n,g,h) + lsigs(mat(n),g,h) * (cden(n) - rcden)
          END DO
      END DO
  END DO


  END SUBROUTINE cden_updt

  !******************************************************************************!

  SUBROUTINE brInterp (rod, mn, xcden, xbcon, xftem, xmtem, sigtr, siga, nuf, &
    sigf, sigs, dc)
  !Purpose: To interpolate the xsec data from xtab file for given bcon,
  ! ftem, mtem and cden

  USE sdata, ONLY: m, XBRANCH, ng
  use io,    only: ounit

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: rod, mn  ! CR indicator and material number
  REAL(DP), INTENT(IN) :: xbcon, xftem, xmtem, xcden  ! TH Parameters
  REAL(DP), DIMENSION(:), INTENT(OUT) :: sigtr, siga, nuf, sigf
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: dc, sigs

  INTEGER, PARAMETER :: nx = 8
  TYPE(XBRANCH), DIMENSION(nx) :: xs   !Temporary xsec for interpolation
  INTEGER :: s, t, u, v, mx
  INTEGER :: s1, s2, t1, t2, u1, u2, v1, v2  ! Dimnesion Position of the given parameters
  INTEGER :: i
  REAL(DP) :: radx

  !Define + and - operators for XBRANCH type addition and substitution respectively
  INTERFACE OPERATOR (+)
    MODULE PROCEDURE brAdd
  END INTERFACE
  INTERFACE OPERATOR (-)
    MODULE PROCEDURE brSubst
  END INTERFACE
  INTERFACE OPERATOR (*)
    MODULE PROCEDURE brRealMult
  END INTERFACE

  ! Set to default
  s1=1; s2=1; t1=1; t2=1; u1=1; u2=1; v1=1; v2=1

  ! Get 2 closest points for interpolation
  ! FOR COOLANT DENSITY
  IF (m(mn)%nd > 1) THEN
    mx = m(mn)%nd
    IF (xcden >= m(mn)%pd(1) .AND. xcden <= m(mn)%pd(mx)) THEN
      DO s = 2, mx
        IF (xcden >= m(mn)%pd(s-1) .AND. xcden <= m(mn)%pd(s)) THEN
          s1 = s-1
          s2 = s
          EXIT
        END IF
      END DO
    ELSE IF (xcden < m(mn)%pd(1) .AND. (m(mn)%pd(1) - xcden) / m(mn)%pd(1) < 0.2) THEN
      s1 = 1
      s2 = 2
    ELSE IF (xcden > m(mn)%pd(mx) .AND. (xcden - m(mn)%pd(mx)) / m(mn)%pd(mx) < 0.2) THEN
      s1 = mx - 1
      s2 = mx
    ELSE
      WRITE(ounit,1567) xcden
      WRITE(*,1567) xcden
      STOP
      1567 FORMAT(2X, '  ERROR: COOLANT DENSITY ', F7.3 ,' g/cm3 &
      & IS OUT OF THE RANGE OF THE BRANCH PARAMETER')
    END IF
  END IF

  ! FOR BORON CONCENTRATION
  IF (m(mn)%nb > 1) THEN
    mx = m(mn)%nb
    IF (xbcon >= m(mn)%pb(1) .AND. xbcon <= m(mn)%pb(mx)) THEN
      DO t = 2, mx
        IF (xbcon >= m(mn)%pb(t-1) .AND. xbcon <= m(mn)%pb(t)) THEN
          t1 = t-1
          t2 = t
          EXIT
        END IF
      END DO
    ELSE IF (xbcon < m(mn)%pb(1) .AND. (m(mn)%pb(1) - xbcon) < 100.) THEN
      t1 = 1
      t2 = 2
    ELSE IF (xbcon > m(mn)%pb(mx) .AND. (xbcon - m(mn)%pb(mx)) < 100.) THEN
      t1 = mx - 1
      t2 = mx
    ELSE
      WRITE(ounit,1568) xbcon
      WRITE(*,1568) xbcon
      STOP
      1568 FORMAT(2X, '  ERROR: BORON CONCENTRATION ', F8.1 ,' PPM &
      & IS OUT OF THE RANGE OF THE BRANCH PARAMETER')
    END IF
  END IF

  ! FOR FUEL TEMPERATURE
  IF (m(mn)%nf > 1) THEN
    mx = m(mn)%nf
    IF (xftem >= m(mn)%pf(1) .AND. xftem <= m(mn)%pf(mx)) THEN
      DO u = 2, mx
        IF (xftem >= m(mn)%pf(u-1) .AND. xftem <= m(mn)%pf(u)) THEN
          u1 = u-1
          u2 = u
          EXIT
        END IF
      END DO
    ELSE IF (xftem < m(mn)%pf(1) .AND. (m(mn)%pf(1) - xftem) / m(mn)%pf(1) < 0.2) THEN
      u1 = 1
      u2 = 2
    ELSE IF (xftem > m(mn)%pf(mx) .AND. (xftem - m(mn)%pf(mx)) / m(mn)%pf(mx) < 0.2) THEN
      u1 = mx - 1
      u2 = mx
    ELSE
      WRITE(ounit,1570) xftem
      WRITE(*,1570) xftem
      STOP
      1570 FORMAT(2X, '  ERROR: FUEL TEMPERATURE ', F7.1 ,' K &
      & IS OUT OF THE RANGE OF THE BRANCH PARAMETER')
    END IF
  END IF

  ! FOR MODERATOR TEMPERATURE
  IF (m(mn)%nm > 1) THEN
    mx = m(mn)%nm
    IF (xmtem >= m(mn)%pm(1) .AND. xmtem <= m(mn)%pm(mx)) THEN
      DO v = 2, mx
        IF (xmtem >= m(mn)%pm(v-1) .AND. xmtem <= m(mn)%pm(v)) THEN
          v1 = v-1
          v2 = v
          EXIT
        END IF
      END DO
    ELSE IF (xmtem < m(mn)%pm(1) .AND. (m(mn)%pm(1) - xmtem) / m(mn)%pm(1) < 0.2) THEN
      v1 = 1
      v2 = 2
    ELSE IF (xmtem > m(mn)%pm(mx) .AND. (xmtem - m(mn)%pm(mx)) / m(mn)%pm(mx) < 0.2) THEN
      v1 = mx - 1
      v2 = mx
    ELSE
      WRITE(ounit,1569) xmtem
      WRITE(*,1569) xmtem
      STOP
      1569 FORMAT(2X, '  ERROR: MODERATOR TEMPERATURE ', F7.1 ,' K &
      & IS OUT OF THE RANGE OF THE BRANCH PARAMETER')
    END IF
  END IF

  !Start doing interpolation
  DO i = 1, nx   !Allocate memory to xs
    ALLOCATE(xs(i)%sigtr(ng))
    ALLOCATE(xs(i)%siga(ng))
    ALLOCATE(xs(i)%nuf(ng))
    ALLOCATE(xs(i)%sigf(ng))
    ALLOCATE(xs(i)%dc(ng,6))
    ALLOCATE(xs(i)%sigs(ng,ng))
  END DO

  IF (rod == 0) THEN   ! For Unrodded XSEC
    !interpolation on Moderator Temperature
    IF (m(mn)%nm > 1) THEN
      radx = (xmtem - m(mn)%pm(v1)) / (m(mn)%pm(v2) - m(mn)%pm(v1))
      xs(1) = m(mn)%xsec(s1,t1,u1,v1) + &
      radx * (m(mn)%xsec(s1,t1,u1,v2) - m(mn)%xsec(s1,t1,u1,v1))
      xs(2) = m(mn)%xsec(s1,t1,u2,v1) + &
      radx * (m(mn)%xsec(s1,t1,u2,v2) - m(mn)%xsec(s1,t1,u2,v1))
      xs(3) = m(mn)%xsec(s1,t2,u1,v1) + &
      radx * (m(mn)%xsec(s1,t2,u1,v2) - m(mn)%xsec(s1,t2,u1,v1))
      xs(4) = m(mn)%xsec(s1,t2,u2,v1) + &
      radx * (m(mn)%xsec(s1,t2,u2,v2) - m(mn)%xsec(s1,t2,u2,v1))
      xs(5) = m(mn)%xsec(s2,t1,u1,v1) + &
      radx * (m(mn)%xsec(s2,t1,u1,v2) - m(mn)%xsec(s2,t1,u1,v1))
      xs(6) = m(mn)%xsec(s2,t1,u2,v1) + &
      radx * (m(mn)%xsec(s2,t1,u2,v2) - m(mn)%xsec(s2,t1,u2,v1))
      xs(7) = m(mn)%xsec(s2,t2,u1,v1) + &
      radx * (m(mn)%xsec(s2,t2,u1,v2) - m(mn)%xsec(s2,t2,u1,v1))
      xs(8) = m(mn)%xsec(s2,t2,u2,v1) + &
      radx * (m(mn)%xsec(s2,t2,u2,v2) - m(mn)%xsec(s2,t2,u2,v1))
    ELSE
      xs(1) = m(mn)%xsec(s1,t1,u1,v1)
      xs(2) = m(mn)%xsec(s1,t1,u2,v1)
      xs(3) = m(mn)%xsec(s1,t2,u1,v1)
      xs(4) = m(mn)%xsec(s1,t2,u2,v1)
      xs(5) = m(mn)%xsec(s2,t1,u1,v1)
      xs(6) = m(mn)%xsec(s2,t1,u2,v1)
      xs(7) = m(mn)%xsec(s2,t2,u1,v1)
      xs(8) = m(mn)%xsec(s2,t2,u2,v1)
    END IF
    !interpolation on Fuel Temperature
    IF (m(mn)%nf > 1) THEN
      radx = (xftem - m(mn)%pf(u1)) / (m(mn)%pf(u2) - m(mn)%pf(u1))
      xs(1) = xs(1) + radx * (xs(2) - xs(1))
      xs(3) = xs(3) + radx * (xs(4) - xs(3))
      xs(5) = xs(5) + radx * (xs(6) - xs(5))
      xs(7) = xs(7) + radx * (xs(8) - xs(7))
    END IF
    !interpolation on Boron concentration
    IF (m(mn)%nb > 1) THEN
      radx = (xbcon - m(mn)%pb(t1)) / (m(mn)%pb(t2) - m(mn)%pb(t1))
      xs(1) = xs(1) + radx * (xs(3) - xs(1))
      xs(5) = xs(5) + radx * (xs(7) - xs(5))
    END IF
    !interpolation on coolant density
    IF (m(mn)%nd > 1) THEN
      xs(1) = xs(1) + (xcden - m(mn)%pd(s1)) / (m(mn)%pd(s2) - m(mn)%pd(s1)) * &
      (xs(5) - xs(1))
    END IF
  ELSE   ! For Rodded XSEC
    !interpolation on Moderator Temperature
    IF (m(mn)%nm > 1) THEN
      radx = (xmtem - m(mn)%pm(v1)) / (m(mn)%pm(v2) - m(mn)%pm(v1))
      xs(1) = m(mn)%rxsec(s1,t1,u1,v1) + &
      radx * (m(mn)%rxsec(s1,t1,u1,v2) - m(mn)%rxsec(s1,t1,u1,v1))
      xs(2) = m(mn)%rxsec(s1,t1,u2,v1) + &
      radx * (m(mn)%rxsec(s1,t1,u2,v2) - m(mn)%rxsec(s1,t1,u2,v1))
      xs(3) = m(mn)%rxsec(s1,t2,u1,v1) + &
      radx * (m(mn)%rxsec(s1,t2,u1,v2) - m(mn)%rxsec(s1,t2,u1,v1))
      xs(4) = m(mn)%rxsec(s1,t2,u2,v1) + &
      radx * (m(mn)%rxsec(s1,t2,u2,v2) - m(mn)%rxsec(s1,t2,u2,v1))
      xs(5) = m(mn)%rxsec(s2,t1,u1,v1) + &
      radx * (m(mn)%rxsec(s2,t1,u1,v2) - m(mn)%rxsec(s2,t1,u1,v1))
      xs(6) = m(mn)%rxsec(s2,t1,u2,v1) + &
      radx * (m(mn)%rxsec(s2,t1,u2,v2) - m(mn)%rxsec(s2,t1,u2,v1))
      xs(7) = m(mn)%rxsec(s2,t2,u1,v1) + &
      radx * (m(mn)%rxsec(s2,t2,u1,v2) - m(mn)%rxsec(s2,t2,u1,v1))
      xs(8) = m(mn)%rxsec(s2,t2,u2,v1) + &
      radx * (m(mn)%rxsec(s2,t2,u2,v2) - m(mn)%rxsec(s2,t2,u2,v1))
    ELSE
      xs(1) = m(mn)%rxsec(s1,t1,u1,v1)
      xs(2) = m(mn)%rxsec(s1,t1,u2,v1)
      xs(3) = m(mn)%rxsec(s1,t2,u1,v1)
      xs(4) = m(mn)%rxsec(s1,t2,u2,v1)
      xs(5) = m(mn)%rxsec(s2,t1,u1,v1)
      xs(6) = m(mn)%rxsec(s2,t1,u2,v1)
      xs(7) = m(mn)%rxsec(s2,t2,u1,v1)
      xs(8) = m(mn)%rxsec(s2,t2,u2,v1)
    END IF
    !interpolation on Fuel Temperature
    IF (m(mn)%nf > 1) THEN
      radx = (xftem - m(mn)%pf(u1)) / (m(mn)%pf(u2) - m(mn)%pf(u1))
      xs(1) = xs(1) + radx * (xs(2) - xs(1))
      xs(3) = xs(3) + radx * (xs(4) - xs(3))
      xs(5) = xs(5) + radx * (xs(6) - xs(5))
      xs(7) = xs(7) + radx * (xs(8) - xs(7))
    END IF
    !interpolation on Boron concentration
    IF (m(mn)%nb > 1) THEN
      radx = (xbcon - m(mn)%pb(t1)) / (m(mn)%pb(t2) - m(mn)%pb(t1))
      xs(1) = xs(1) + radx * (xs(3) - xs(1))
      xs(5) = xs(5) + radx * (xs(7) - xs(5))
    END IF

    !interpolation on coolant density
    IF (m(mn)%nd > 1) THEN
      xs(1) = xs(1) + (xcden - m(mn)%pd(s1)) / (m(mn)%pd(s2) - m(mn)%pd(s1)) * &
      (xs(5) - xs(1))
    END IF
  END IF
  sigtr = xs(1)%sigtr
  siga = xs(1)%siga
  nuf = xs(1)%nuf
  sigf = xs(1)%sigf
  sigs = xs(1)%sigs
  dc   = xs(1)%dc

  DO i = 1, nx   !DeAllocate memory to xs
    DEALLOCATE(xs(i)%sigtr)
    DEALLOCATE(xs(i)%siga)
    DEALLOCATE(xs(i)%nuf)
    DEALLOCATE(xs(i)%sigf)
    DEALLOCATE(xs(i)%dc)
    DEALLOCATE(xs(i)%sigs)
  END DO


  END SUBROUTINE brInterp

  !******************************************************************************!

  FUNCTION brAdd(A, B) RESULT (C)

    ! To perform XBRANCH data type addition

  USE sdata, ONLY: XBRANCH

  IMPLICIT NONE

  TYPE(XBRANCH), INTENT(IN) :: A, B
  TYPE(XBRANCH) :: C

  C%sigtr = A%sigtr + B%sigtr
  C%siga = A%siga + B%siga
  C%nuf = A%nuf + B%nuf
  C%sigf = A%sigf + B%sigf
  C%sigs = A%sigs + B%sigs
  C%dc = A%dc + B%dc

  END FUNCTION brAdd

  !******************************************************************************!

  FUNCTION brSubst(A, B) RESULT (C)

      ! To perform XBRANCH data type substraction

  USE sdata, ONLY: XBRANCH

  IMPLICIT NONE

  TYPE(XBRANCH), INTENT(IN) :: A, B
  TYPE(XBRANCH) :: C

  C%sigtr = A%sigtr - B%sigtr
  C%siga = A%siga - B%siga
  C%nuf = A%nuf - B%nuf
  C%sigf = A%sigf - B%sigf
  C%sigs = A%sigs - B%sigs
  C%dc = A%dc - B%dc

  END FUNCTION brSubst

  !******************************************************************************!

  FUNCTION brRealMult(Re, A) RESULT (B)

      ! To perform XBRANCH data type substraction

  USE sdata, ONLY: XBRANCH, DP

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: Re
  TYPE(XBRANCH), INTENT(IN) :: A
  TYPE(XBRANCH) :: B

  B%sigtr = Re * A%sigtr
  B%siga = Re * A%siga
  B%nuf = Re * A%nuf
  B%sigf = Re * A%sigf
  B%sigs = Re * A%sigs
  B%dc = Re * A%dc

  END FUNCTION brRealMult

  !******************************************************************************!




end module
