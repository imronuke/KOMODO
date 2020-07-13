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

    use sdata, only: nnod, f0, aprad, apaxi, afrad, ftem, mtem, cden, bcon, bpos
    use io,    only: AsmPow, AxiPow, AsmFlux, inp_read
    use xsec,  only: XS_updt
    use cmfd,  only: outer, powdis

    IMPLICIT NONE

    REAL(DP), DIMENSION(:), ALLOCATABLE :: pow

    !Update xsec
    CALL XS_updt(bcon, ftem, mtem, cden, bpos)

    call print_head()

    !Outer iteration
    CALL outer(1)

    IF (aprad == 1 .OR. apaxi == 1) THEN
        ALLOCATE(pow(nnod))
        CALL PowDis(pow)
    END IF

    IF (aprad == 1) CALL AsmPow(pow)

    IF (apaxi == 1) CALL AxiPow(pow)

    IF (afrad == 1) CALL AsmFlux(f0, 1.e0_DP)


  END SUBROUTINE forward

  !******************************************************************************!

  SUBROUTINE adjoint()

    !
    ! Purpose:
    !    To solve adjoint problems
    !

    use sdata, only: nnod, f0, aprad, apaxi, afrad, ftem, mtem, cden, bcon, bpos
    use io,    only: AsmPow, AxiPow, AsmFlux, inp_read
    use xsec,  only: XS_updt
    use cmfd,  only: outer_ad, powdis

    IMPLICIT NONE

    REAL(DP), DIMENSION(:), ALLOCATABLE :: pow

    !Update xsec
    CALL XS_updt(bcon, ftem, mtem, cden, bpos)

    call print_head()

    !Outer iteration
    CALL outer_ad(1)

    IF (aprad == 1 .OR. apaxi == 1) THEN
        ALLOCATE(pow(nnod))
        CALL PowDis(pow)
    END IF

    IF (aprad == 1) CALL AsmPow(pow)

    IF (apaxi == 1) CALL AxiPow(pow)

    IF (afrad == 1) CALL AsmFlux(f0, 1.e0_DP)

  END SUBROUTINE adjoint

  !******************************************************************************!

  SUBROUTINE fixedsrc()

    !
    ! Purpose:
    !    To solve fixed source problems
    !

    use sdata, only: nnod, f0, aprad, apaxi, afrad, ftem, mtem, cden, bcon, bpos
    use io,    only: AsmPow, AxiPow, AsmFlux, inp_read
    use xsec,  only: XS_updt
    use cmfd,  only: outer_fs, powdis

    IMPLICIT NONE

    REAL(DP), DIMENSION(:), ALLOCATABLE :: pow

    !Update xsec
    CALL XS_updt(bcon, ftem, mtem, cden, bpos)

    call print_head()

    !Outer iteration
    CALL outer_fs(1)

    IF (aprad == 1 .OR. apaxi == 1) THEN
        ALLOCATE(pow(nnod))
        CALL PowDis(pow)
    END IF

    IF (aprad == 1) CALL AsmPow(pow)

    IF (apaxi == 1) CALL AxiPow(pow)

    IF (afrad == 1) CALL AsmFlux(f0)

  END SUBROUTINE fixedsrc

  !******************************************************************************!

  SUBROUTINE print_head()

    !
    ! Purpose:
    !    To print header
    !

    use sdata, only: mode
    use io,    only: ounit, scr

    IMPLICIT NONE

    if (mode == 'FORWARD' .or. mode == 'ADJOINT') then
      WRITE(ounit,*); WRITE(ounit,*)
      WRITE(ounit,3245); WRITE(ounit,3247); WRITE(ounit,3245)
      WRITE(ounit,*)
      WRITE(ounit,3248); WRITE(ounit,3249)
      if (scr) then
        WRITE(*,*); WRITE(*,*)
        WRITE(*,3245); WRITE(*,3247); WRITE(*,3245)
        WRITE(*,*)
        WRITE(*,3248); WRITE(*,3249)
      end if
    else
      WRITE(ounit,*); WRITE(ounit,*)
      WRITE(ounit,3245); WRITE(ounit,3247); WRITE(ounit,3245)
      WRITE(ounit,*)
      WRITE(ounit,3248); WRITE(ounit,3249)
      if (scr) then
        WRITE(*,*); WRITE(*,*)
        WRITE(*,3245); WRITE(*,3247); WRITE(*,3245)
        WRITE(*,*)
        WRITE(*,3250); WRITE(*,3249)
      end if
    end if



    3245 format (1X, ' ==============================================', &
    '================================')
    3247 format(27X,'CALCULATION RESULTS')
    3248 format(2X,'Itr     k-eff     Fis.Src Error   Inner Error')
    3249 format(1X,'----------------------------------------------------')
    3250 format(2X,'Itr   Fis.Src Error   Inner Error')


  END SUBROUTINE print_head



end module
