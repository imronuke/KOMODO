MODULE io

!=========================
! Input output module to read, process and echo input, as well as writing output
! =======================

USE sdata, ONLY: DP

IMPLICIT NONE

SAVE

! ind is sed to read x indicator in beginning of input buffer line.
! This to prevent reading next line
CHARACTER(LEN=1) :: ind

CHARACTER(LEN=200) :: message    ! error message
CHARACTER(LEN=200):: iline       ! Input line
!Ouput options
LOGICAL, PARAMETER :: ogeom = .TRUE.  ! Geometry output print option
LOGICAL, PARAMETER :: oxsec = .TRUE.  ! Macroscopic CXs output print option
LOGICAL, PARAMETER :: scr = .TRUE.    ! Terminal ouput print option

! Input, output and buffer input file unit number
INTEGER, PARAMETER :: iunit = 100   !input file unit number
INTEGER, PARAMETER :: ounit = 101   !output file unit number
INTEGER, PARAMETER :: buff  = 99    !input buffer file unit number (entire input)

! Input buffer file unit number for each card
INTEGER, PARAMETER :: umode = 111, uxsec = 112, ugeom = 113, ucase = 114
INTEGER, PARAMETER :: uesrc = 115, uiter = 118, uprnt = 119, uadf  = 120
INTEGER, PARAMETER :: ucrod = 121, ubcon = 122, uftem = 123, umtem = 124
INTEGER, PARAMETER :: ucden = 125, ucbcs = 126, uejct = 127, uther = 128
INTEGER, PARAMETER :: uxtab = 129, ukern = 130, uextr = 131, uthet = 132
INTEGER :: bunit

! Card active/inactive indicator (active = 1, inactive = 0)
INTEGER :: bmode = 0, bxsec = 0, bgeom = 0, bcase = 0, besrc = 0
INTEGER :: biter = 0, bprnt = 0, badf  = 0, bcrod = 0, bbcon = 0
INTEGER :: bftem = 0, bmtem = 0, bcden = 0, bcbcs = 0, bejct = 0
INTEGER :: bther = 0, bxtab = 0, bkern = 0, bextr = 0, bthet = 0

! This declaration is to notify that the error in separated card file (in case so)
INTEGER, PARAMETER :: ncard = 20                  ! Number of card
INTEGER, DIMENSION(ncard) :: uarr = &             ! Array of buffer unit number
(/umode, uxsec, ugeom, ucase, uesrc, uiter, uprnt, uadf,  ucrod, ubcon, &
  uftem, umtem, ucden, ucbcs, uejct, uther, uxtab, ukern, uextr, uthet /)
CHARACTER(LEN=4), DIMENSION(ncard) :: carr = &    ! Array of card name
(/'MODE', 'XSEC', 'GEOM', 'CASE', 'ESRC', 'ITER', 'PRNT', 'ADF ', 'CROD', 'BCON', &
  'FTEM', 'MTEM', 'CDEN', 'CBCS', 'EJCT', 'THER', 'XTAB', 'KERN', 'EXTR', 'THET' /)
CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: farr   ! Array of card file

! Geometry
INTEGER :: np                                           ! Number of planars
INTEGER, DIMENSION(:), ALLOCATABLE :: zpln              ! Planar assignment to z direction
REAL(DP), DIMENSION(:), ALLOCATABLE :: xsize, ysize, zsize  !Assembly size for each direction
TYPE :: MAT_ASGN                                        ! Material assignment
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: asm         ! Material assignment into assembly
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: node        ! Material assignment into nodes
END TYPE
TYPE(MAT_ASGN), DIMENSION(:), ALLOCATABLE :: plnr       ! planar
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: mnum


CONTAINS

!******************************************************************************!

SUBROUTINE inp_read()
!
! Purpose:
!    [Main subroutine in this module] To read input, echo the
!    input and gives the description
!    to the user about his/her input
!


USE sdata, ONLY: ng, nnod, mode, dc, exsrc

IMPLICIT NONE

INTEGER :: g, i, N
CHARACTER(LEN=100) :: iname, oname

!Got this trick from: http://web.utah.edu/thorne/computing/Handy_Fortran_Tricks.pdf
N = IARGC()
IF (N < 1) THEN
   WRITE(*,*) '  NOTE : You can also write the input directly after the command'
   WRITE(*,'(A,A100)',ADVANCE='NO') '  INPUT NAME : '
   READ(*,*) iname
ELSE
   CALL GETARG(1,iname) !Grab the first command line argument
ENDIF

iname = TRIM(iname)

CALL openFIle (iunit, iname, 'input', 'Input File Open Failed--status')

oname = TRIM(iname) // '.out'
oname = TRIM(oname)

OPEN (UNIT=ounit, FILE=oname, STATUS='REPLACE', ACTION='WRITE')

OPEN (UNIT=umode, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=uxsec, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=ugeom, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=ucase, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=uesrc, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=uiter, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=uprnt, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=uadf,  STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=ucrod, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=ubcon, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=uftem, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=umtem, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=ucden, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=ucbcs, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=uejct, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=uther, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=uxtab, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=ukern, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=uextr, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=uthet, STATUS='SCRATCH', ACTION='READWRITE')

! By default, card file names are the input file name
ALLOCATE (farr(ncard)); farr = ADJUSTL(iname)

! Echo the input to the output file
CALL inp_echo()

! Remove comments and write the entire input to the buffer file
CALL inp_comments (iunit, buff, '!')

! Break the buffer file and re-write into different input card buffer
CALL inp_rewrite(buff)

! Back to the first line for all input card buffer
REWIND(umode); REWIND(uxsec); REWIND(ugeom); REWIND(ucase); REWIND(uesrc)
REWIND(uiter); REWIND(uprnt); REWIND(uadf); REWIND(ucrod); REWIND(ubcon)
REWIND(uftem); REWIND(umtem); REWIND(ucden); REWIND(ucbcs); REWIND(uejct)
REWIND(uther); REWIND(uxtab); REWIND(ukern); REWIND(uextr); REWIND(uthet)

! Start reading buffer files for each input card buffer

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,1008)
WRITE(ounit,*) &
' ******************************************************************************'

! Card MODE
IF (bmode == 1) THEN
    CALL inp_mode(umode)
ELSE
    WRITE(ounit,1021) '%MODE'
    STOP
END IF

! Card CASE
IF (bcase == 1) CALL inp_case (ucase)

! Card %KERN
if (bkern==1) then
  call inp_kern(ukern)
else
  if (scr) then
    WRITE(*,*)
    WRITE(*,*) ' NODAL KERNEL  : SEMI-ANALYTIC NODAL METHOD'
  end if
end if

! Card XSEC
IF (bxsec == 1) THEN
    CALL inp_xsec(uxsec)
ELSE IF (bxtab == 1) THEN
    CALL inp_xtab(uxtab)
ELSE
    WRITE(ounit,1021) '%XSEC OR %XTAB'
    WRITE(*,1021) '%XSEC OR %XTAB'
    STOP
END IF

! Card GEOM
IF (bgeom == 1) THEN
    CALL inp_geom1(ugeom)
    CALL inp_geom2(ugeom)
ELSE
    WRITE(ounit,1021) '%GEOM'
    WRITE(*,1021) '%GEOM'
    STOP
END IF

! Card PRNT
IF (bprnt == 1) CALL inp_prnt (uprnt)

! Card CBCS
IF (mode == 'BCSEARCH' .AND. bcbcs == 1) THEN
    CALL inp_cbcs(ucbcs)
ELSE IF (mode == 'BCSEARCH' .AND. bcbcs == 0 .AND. bxtab == 0) THEN
    WRITE(ounit,*) '   ERROR: CALCULATION MODE IS CRITICAL BORON CONCENTRATION SEARCH'
    WRITE(ounit,1041) 'CBCS', 'CRITICAL BORON CONCENTRATION SEARCH'
    WRITE(*,*) '   ERROR: CALCULATION MODE IS CRITICAL BORON CONCENTRATION SEARCH'
    WRITE(*,1041) 'CBCS', 'CRITICAL BORON CONCENTRATION SEARCH'
    STOP
ELSE IF (mode /= 'BCSEARCH' .AND. bcbcs == 1) THEN
    WRITE(ounit,*) '   ERROR: CBCS CARD IS NOT NECESSARY FOR THIS CALCULATION MODE'
    WRITE(*,*) '   ERROR: CBCS CARD IS NOT NECESSARY FOR THIS CALCULATION MODE'
    STOP
ELSE IF (mode == 'BCSEARCH' .AND. bbcon == 1) THEN
    WRITE(ounit,*) '   ERROR: BCON CARD MUST NOT PRESENT FOR THIS CALCULATION MODE'
    WRITE(*,*) '   ERROR: BCON CARD MUST NOT PRESENT FOR THIS CALCULATION MODE'
    STOP
ELSE
    CONTINUE
END IF

IF (mode == 'BCSEARCH' .AND. bxtab == 1 .AND. bther == 0) THEN
  WRITE(ounit,*) '   ERROR: THER CARD MUST BE PRESENT IN THIS PROBLEM'
  WRITE(*,*) '   ERROR: THER CARD MUST BE PRESENT IN THIS PROBLEM'
  STOP
END IF


!CARD CROD
IF (bcrod == 1) CALL inp_crod (ucrod)

! Card EJCT (Rod Ejection)
IF (mode == 'RODEJECT' .AND. bejct == 1 .AND. bcrod == 1) THEN
    CALL inp_ejct(uejct)
ELSE IF (mode == 'RODEJECT' .AND. bejct /= 1) THEN
    WRITE(ounit,*) '   CALCULATION MODE ROD EJECTION'
    WRITE(ounit,1041) 'EJCT', 'ROD EJECTION - TRANSIENT'
    WRITE(*,*) '   CALCULATION MODE ROD EJECTION'
    WRITE(*,1041) 'EJCT', 'ROD EJECTION - TRANSIENT'
    STOP
ELSE IF (mode == 'RODEJECT' .AND. bcrod /= 1) THEN
    WRITE(ounit,*) '   CALCULATION MODE ROD EJECTION'
    WRITE(ounit,1041) 'CROD', 'CONTROL ROD'
    WRITE(*,*) '   CALCULATION MODE ROD EJECTION'
    WRITE(*,1041) 'CROD', 'CONTROL ROD'
    STOP
ELSE IF (mode /= 'RODEJECT' .AND. bejct == 1) THEN
    WRITE(ounit,*) '   EJCT CARD IS NOT NECESSARY FOR THIS CALCULATION MODE'
    WRITE(*,*) '   EJCT CARD IS NOT NECESSARY FOR THIS CALCULATION MODE'
    STOP
ELSE
    CONTINUE
END IF

! Miscellaneous things
CALL misc()

! Card ITER
IF (biter == 1) CALL inp_iter (uiter)

! Card THET
IF (bthet == 1) CALL inp_thet (uthet)

!!CARD THER
IF (bther == 1 .AND. mode == 'FORWARD') THEN
    WRITE(ounit,*)'   ERROR: %THER CARD NOT VALID FOR FORWARD CALCULATION MODE'
    WRITE(*,*)'   ERROR: %THER CARD NOT VALID FOR FORWARD CALCULATION MODE'
    STOP
ELSE IF (bther == 1 .AND. mode == 'FIXEDSRC') THEN
    WRITE(ounit,*)'   ERROR: %THER CARD NOT VALID FOR FIXED SOURCE CALCULATION MODE'
    WRITE(*,*)'   ERROR: %THER CARD NOT VALID FOR FIXED SOURCE CALCULATION MODE'
    STOP
ELSE IF (bther == 1 .AND. mode == 'ADJOINT') THEN
    WRITE(ounit,*)'   ERROR: %THER CARD NOT VALID FOR ADJOINT CALCULATION MODE'
    WRITE(*,*)'   ERROR: %THER CARD NOT VALID FOR ADJOINT CALCULATION MODE'
    STOP
ELSE IF (bther == 1 .AND. bftem == 1 .AND. (bmtem ==1 .OR. bcden == 1)) THEN
    CALL inp_ther (uther)
ELSE IF (bther == 1 .AND. bxtab == 1) THEN
      CALL inp_ther (uther)
ELSE IF (bther == 0) THEN
    CONTINUE
ELSE
    IF (bxtab /= 1) THEN
      WRITE(ounit,*)'   ERROR: WHEN %THER CARD PRESENT %FTEM AND,' // &
                    'AT LEAST ONE OF THE FOLLOWING CARDS MUST PRESENT'
      WRITE(ounit,*)'   1. %MTEM   2. %CDEN'
      WRITE(*,*)'   ERROR: WHEN %THER CARD PRESENT %FTEM AND,' // &
                    'AT LEAST ONE OF THE FOLLOWING CARDS MUST PRESENT'
      WRITE(*,*)'   1. %MTEM   2. %CDEN'
      STOP
    END IF
END IF

if (bextr == 1) call inp_extr()

!!CARD BCON
IF (bbcon == 1 .AND. bxtab == 0) CALL inp_bcon (ubcon)
IF (bbcon == 1 .AND. bxtab == 1 .AND. mode == 'RODEJECT') CALL inp_bcon (ubcon)

!!CARD FTEM
IF (bftem == 1 .AND. bxtab == 0) CALL inp_ftem (uftem)

!!CARD MTEM
IF (bmtem == 1 .AND. bxtab == 0) CALL inp_mtem (umtem)

!!CARD CDEN
IF (bcden == 1 .AND. bxtab == 0) CALL inp_cden (ucden)


!CARD ADF
allocate(dc(nnod,ng,6))
DO g = 1, ng
    DO n = 1, nnod
        dc(n,g,1:) = 1._dp     !by default, adf = 1
    END DO
END DO

IF (badf == 1 .AND. bxtab == 0) THEN
  CALL inp_adf (uadf)
ELSE IF (badf == 1 .AND. bxtab == 1) THEN
  WRITE(ounit,*) '  BOTH %ADF AND %XTAB CARDS CANNOT PRESENT TOGETHER'
  WRITE(*,*) '  BOTH %ADF AND %XTAB CARDS CANNOT PRESENT TOGETHER'
  STOP
ELSE
  CONTINUE
END IF

! Card ESRC
ALLOCATE(exsrc(nnod, ng))   ! For transient or rod ejection problem, this used to
exsrc = 0._DP               ! store transient terms
IF (mode == 'FIXEDSRC' .AND. besrc == 1) THEN
    CALL inp_esrc(uesrc)
ELSE IF (mode == 'FIXEDSRC' .AND. besrc /= 1) THEN
    WRITE(ounit,*) '   CALCULATION MODE IS FIXED SOURCE'
    WRITE(ounit,1041) 'ESRC', 'FIXED SOURCE'
    WRITE(*,*) '   CALCULATION MODE IS FIXED SOURCE'
    WRITE(*,1041) 'ESRC', 'FIXED SOURCE'
    STOP
ELSE IF (mode /= 'FIXEDSRC' .AND. besrc == 1) THEN
    WRITE(ounit,*) '   ESRC CARD IS NOT NECESSARY FOR THIS CALCULATION MODE'
    WRITE(*,*) '   ESRC CARD IS NOT NECESSARY FOR THIS CALCULATION MODE'
    STOP
ELSE
    CONTINUE
END IF



DEALLOCATE(mnum)
DO i= 1,np
    DEALLOCATE(plnr(i)%asm)
    DEALLOCATE(plnr(i)%node)
END DO
DEALLOCATE(plnr)
DEALLOCATE(zpln)

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) &
' ************************', '  STOP READING INPUT  ', '********************************'

1008 FORMAT (30X, 'START READING INPUT')
1021 FORMAT(2X, 'CARD ', A, ' DOES NOT PRESENT. THIS CARD IS MANDATORY')
1041 FORMAT(2X, 'CARD ', A, ' DOES NOT PRESENT. THIS CARD IS MANDATORY FOR ', A,' CALCULATION MODE')


CLOSE(UNIT=umode); CLOSE(UNIT=uxsec); CLOSE(UNIT=ugeom); CLOSE(UNIT=ucase)
CLOSE(UNIT=uesrc); CLOSE(UNIT=uiter); CLOSE(UNIT=uprnt); CLOSE(UNIT=uadf)
CLOSE(UNIT=ucrod); CLOSE(UNIT=ubcon); CLOSE(UNIT=uftem); CLOSE(UNIT=umtem)
CLOSE(UNIT=ucden); CLOSE(UNIT=ucbcs); CLOSE(UNIT=uejct); CLOSE(UNIT=uther)
CLOSE(UNIT=uxtab); CLOSE(UNIT=ukern); CLOSE(UNIT=uextr); CLOSE(UNIT=uthet)
CLOSE(UNIT=buff)


END SUBROUTINE inp_read

!******************************************************************************!

SUBROUTINE openFile(iunit, iname, file, message)

INTEGER :: iunit
CHARACTER(LEN=*) :: iname, file, message
INTEGER  :: iost

OPEN (UNIT=iunit, FILE=iname, STATUS='OLD', ACTION='READ', &
      IOSTAT = iost)

IF (iost /= 0) THEN
  WRITE(*,1020) message, iost
  WRITE(*,*) '  CANNOT FIND '// file //' FILE : ', iname
  1020 FORMAT    (2X, A, I6)
  STOP
END IF

END SUBROUTINE openFile

!******************************************************************************!

SUBROUTINE inp_echo()
!
! Purpose:
!    To rewrite the input
!

IMPLICIT NONE

INTEGER :: eof
INTEGER :: nline

WRITE(ounit, 2409)
WRITE(ounit, 2411)
WRITE(ounit, 2412)
WRITE(ounit, 2409)
WRITE(ounit, *)
WRITE(ounit, *)

if (scr) then
  WRITE(*, *)
  WRITE(*, *)
  WRITE(*, 2409)
  WRITE(*, 2411)
  WRITE(*, 2412)
  WRITE(*, 2409)
  WRITE(*, *)
  WRITE(*, *)
end if


WRITE(ounit,1002) 'STARTS'

nline = 0
DO
   READ(iunit, '(A200)', IOSTAT=eof) iline
   nline = nline + 1
   IF (eof < 0) THEN
       WRITE(ounit,1002) 'ENDS'
       WRITE(ounit,*)
       EXIT
   END IF
   WRITE(ounit, 1001) nline, iline
END DO

2409 FORMAT(11X, '###########################################################')
2411 FORMAT(11X, '#                     ADPRES 1.2                          #')
2412 FORMAT(11X, '#        ABU DHABI POLYTECHNIC REACTOR SIMULATOR          #')

1001 FORMAT (2X, I4, ': ', A200)
1002 FORMAT    (2X, '=============================INPUT DATA',A7, &
                    ' HERE===========================')
REWIND (iunit)

END SUBROUTINE inp_echo

!******************************************************************************!

SUBROUTINE inp_comments (inunit, buffer, mark)
!
! Purpose:
!    To remove the comments in input and rewrite the
!    input into input buffer. Comments marked by !.
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: inunit, buffer
CHARACTER(LEN=*), INTENT(IN) :: mark

INTEGER :: ln                  ! line number
INTEGER :: eof, comm

OPEN (UNIT=buffer, STATUS='SCRATCH', ACTION='READWRITE')

! Start removing comments and rewrite into one input buffer
ln = 0
DO
    ln = ln+1
    READ (inunit, '(A200)', IOSTAT=eof) iline
    IF (eof < 0) EXIT              !Check end of file
    iline = TRIM(ADJUSTL(iline))   ! Remove trailing blanks and adjust to left
    comm = INDEX(iline, mark)       ! Find position '!' if any
    ! If there is no '!' and no first 20 blank spaces (in case line is blank)
    IF (comm == 0 .AND. iline(1:20) /= '                    ')  THEN
        WRITE(buffer,1012)'x ',ln,iline
    END IF
    !If the first character is not '!'
    IF (comm > 1) THEN
        iline = iline(1:comm-1)       ! Take only part of input
        WRITE(buffer,1012)'x ',ln, iline
    END IF
END DO

REWIND(buffer)

1012 FORMAT(A2, I5,' ',A200)

END SUBROUTINE inp_comments

!******************************************************************************!

SUBROUTINE inp_rewrite (buffer)

! Purpose:
! Read previous input buffer and rewrite and break into different input buffer
! for particular cards. Cards identfied by %

IMPLICIT NONE

INTEGER, INTENT(IN) :: buffer

INTEGER :: ln                  ! line number
INTEGER :: eof, per, comm
CHARACTER(LEN=100) :: fname    ! Card File name
CHARACTER(LEN=4)   :: card    ! Card name
INTEGER, PARAMETER :: cunit = 996  !XTAB file unit number
INTEGER, PARAMETER :: xunit = 997  !XTAB Buffer unit number

DO
    per = 0
    READ (buffer, '(A2,I5,A200)', IOSTAT=eof) ind, ln, iline
    IF (eof < 0) EXIT                                                   !Check end of file

    ! If the card is placed in a separated file as indicated by keyword FILE
    IF (INDEX(iline,'FILE') > 0 .OR. INDEX(iline,'file') > 0) then
      BACKSPACE(buffer)
      READ (buffer, '(A2,I5,A200)') ind, ln, iline
      iline = ADJUSTL(iline)                                             ! Adjust to left
      comm = INDEX(iline, ' ')                                           ! Get space position
      fname = TRIM(ADJUSTL(iline(comm+1:200)))                           ! Get card file name
      farr(GETLOC(uarr, bunit)) = fname                                 ! Change default file name for error notification
      CALL openFile(xunit, fname, card, 'CARD File Open Failed--status') ! Open card file
      CALL inp_comments(xunit, cunit, '!')                               ! Remove comments in card file
      ! Begin read card in a separated file
      DO
        READ (cunit, '(A2,I5,A200)', IOSTAT=eof) ind, ln, iline
        IF (eof < 0) EXIT                                                !Check end of file
        WRITE(bunit, 1019) 'x ',ln, iline                                ! 'x' used to prevent reading next line
      END DO
      CLOSE(xunit); CLOSE(cunit)
    END IF

    per = INDEX(iline,'%')
    IF (per > 0) THEN              ! IF %card detected
      iline = iline(per+1:200)
      iline = TRIM(ADJUSTL(iline))
      card  = iline
      SELECT CASE (iline)
      CASE('MODE'); bunit = umode; bmode = 1
      CASE('XSEC'); bunit = uxsec; bxsec = 1
      CASE('GEOM'); bunit = ugeom; bgeom = 1
      CASE('CASE'); bunit = ucase; bcase = 1
      CASE('ESRC'); bunit = uesrc; besrc = 1
      CASE('ITER'); bunit = uiter; biter = 1
      CASE('PRNT'); bunit = uprnt; bprnt = 1
      CASE('ADF') ; bunit = uadf ; badf  = 1
      CASE('CROD'); bunit = ucrod; bcrod = 1
      CASE('EJCT'); bunit = uejct; bejct = 1
      CASE('CBCS'); bunit = ucbcs; bcbcs = 1
      CASE('FTEM'); bunit = uftem; bftem = 1
      CASE('MTEM'); bunit = umtem; bmtem = 1
      CASE('CDEN'); bunit = ucden; bcden = 1
      CASE('BCON'); bunit = ubcon; bbcon = 1
      CASE('THER'); bunit = uther; bther = 1
      CASE('XTAB'); bunit = uxtab; bxtab = 1
      CASE('KERN'); bunit = ukern; bkern = 1
      CASE('EXTR'); bunit = uextr; bextr = 1
      CASE('THET'); bunit = uthet; bthet = 1
      CASE DEFAULT
        WRITE(ounit,1014) ln, iline
        WRITE(*,1014) ln, iline
        STOP
      END SELECT
    END IF
    ! Write input buffer for each card
    IF (per == 0) WRITE(bunit, 1019) 'x ',ln, iline  ! 'x' used to prevent reading next line
END DO

1014 FORMAT(2X, 'AT LINE', I3, ' : THIS IS A WRONG INPUT CARD : ', A8)
1019 FORMAT(A2, I5,' ',A200)

END SUBROUTINE inp_rewrite

!******************************************************************************!

SUBROUTINE inp_mode (xbunit)
!
! Purpose:
!    To read case mode in input
!

USE sdata, ONLY: mode

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status
CHARACTER(LEN=60) :: mode_desc

READ(xbunit, '(A2, I5,A100)', IOSTAT=ios) ind, ln, mode
message = ' error in reading MODE'
CALL er_message(ounit, ios, ln, message, buf=xbunit)


mode = TRIM(ADJUSTL(mode))  ! ADJUSTL = MOVE PRECEDING BLANK TO TRAILING

SELECT CASE(mode)
CASE('FORWARD')
    mode_desc = TRIM(ADJUSTL('FORWARD CALCULATION'))
CASE('ADJOINT')
  mode_desc = TRIM(ADJUSTL('ADJOINT CALCULATION'))
CASE('FIXEDSRC')
  mode_desc = TRIM(ADJUSTL('FIXED SOURCE CALCULATION'))
CASE('RODEJECT')
  if (bther == 0) then
    mode_desc = TRIM(ADJUSTL('ROD EJECTION CALCULATION WITHOUT T-H FEEDBACK'))
  else
    mode_desc = TRIM(ADJUSTL('ROD EJECTION CALCULATION WITH T-H FEEDBACK'))
  end if
CASE('BCSEARCH')
  if (bther == 0) then
    mode_desc = TRIM(ADJUSTL('CRITICAL BORON CONCENTRATION SEARCH' &
    // ' WITHOUT T-H FEEDBACK'))
  else
    mode_desc = TRIM(ADJUSTL('CRITICAL BORON CONCENTRATION SEARCH' &
    // ' WITH T-H FEEDBACK'))
  end if
CASE DEFAULT
  WRITE(ounit,1032) mode
  WRITE(*,1032) mode
  STOP
END SELECT

WRITE(ounit,1031) mode_desc
WRITE(ounit,*)
if (scr) then
  WRITE(*,1031) mode_desc
  WRITE(*,*)
end if

1031 FORMAT(2X, 'CALCULATION MODE : ', A60)
1032 FORMAT(2X, 'MODE : ', A10, ' IS UNIDENTIFIED')

END SUBROUTINE inp_mode

!******************************************************************************!

SUBROUTINE inp_case (xbunit)
!
! Purpose:
!    To read case card in input
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

CHARACTER(LEN=100) :: case_id
CHARACTER(LEN=100) :: case_exp

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

READ(xbunit, '(A2, I5,A100)', IOSTAT=ios) ind, ln, case_id
message = ' error in CASE ID'
CALL er_message(ounit, ios, ln, message, buf=xbunit)
READ(xbunit, '(A2, I5,A100)', IOSTAT=ios) ind, ln, case_exp
message = ' error in case description'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

case_id = TRIM(ADJUSTL(case_id))
case_exp = TRIM(ADJUSTL(case_exp))

WRITE(ounit,1006) case_id
WRITE(ounit,1007) case_exp
if (scr) then
  WRITE(*,1006) case_id
  WRITE(*,1007) case_exp
end if
1006 FORMAT(2X, 'CASE ID : ', A100)
1007 FORMAT(2X, A100)

END SUBROUTINE inp_case

!******************************************************************************!

SUBROUTINE inp_xsec (xbunit)
!
! Purpose:
!    To read CROSS SECTIONS card in input
!

USE sdata, ONLY: nmat, ng, xsigtr, xsiga, xnuf, xsigf, &
                 xsigs, xD, xsigr, chi, mode, ccnuf, ccsigf

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: i, g, h
INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status
REAL(DP) :: dum
INTEGER, DIMENSION(:), ALLOCATABLE :: group
INTEGER :: comm1, comm2

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>> READING MACROSCOPIC CROSS SECTIONS <<<<'
WRITE(ounit,*) '           --------------------------------------------'

READ(xbunit, *, IOSTAT=ios) ind, ln, ng, nmat
message = ' error in material number'
CALL er_message(ounit, ios, ln, message, buf=xbunit)


ALLOCATE(group(ng))
DO g = 1,ng
   group(g) = g
END DO

ALLOCATE(xsigtr(nmat,ng))
ALLOCATE(xsiga (nmat,ng))
ALLOCATE(xnuf  (nmat,ng))
ALLOCATE(xsigf (nmat,ng))
ALLOCATE(xsigs (nmat,ng,ng))
ALLOCATE(xD    (nmat,ng))
ALLOCATE(xsigr (nmat,ng))
ALLOCATE(chi  (nmat,ng))

! To ancticipate users make mistake when they use %XTAB instead of %XSEC
READ(xbunit, '(A200)') iline
comm1 = INDEX(iline, '/')
comm2 = INDEX(iline, '\')
IF (comm1 > 0 .OR. comm2 > 0) THEN
  WRITE(ounit, *) '  ERROR: SLASH (/) OR BACKSLASH (\) NOT ACCEPTED IN %XSEC CARD '
  WRITE(*, *) '  ERROR: SLASH (/) OR BACKSLASH (\) NOT ACCEPTED IN %XSEC CARD '
  STOP
END IF
BACKSPACE(xbunit)

! Reading MACROSCOPIC CXs
DO i= 1, nmat
    DO g= 1, ng
        READ(xbunit, *, IOSTAT=ios) ind, ln, xsigtr(i,g), &
        xsiga(i,g), xnuf(i,g), xsigf(i,g), &
        chi(i,g), (xsigs(i,g,h), h = 1, ng)
        message = ' error in cross section data'
        CALL er_message(ounit, ios, ln, message, buf=xbunit)

        ! Check CXs values
        IF (xsigtr(i,g) <= 0.0) THEN
            WRITE(ounit,1020)i, g
            STOP
        END IF
        IF (xnuf(i,g) > 0.) ccnuf = .FALSE.
        IF (xsigf(i,g) > 0.) ccsigf = .FALSE.

        xD(i,g) = 1./(3.*xsigtr(i,g))
        dum = 0.0
        DO h= 1, ng
            IF (g /= h) dum = dum + xsigs(i,g,h)
        END DO
        xsigr(i,g) =  xsiga(i,g) + dum
    END DO
END DO


! Writing output
IF (oxsec) THEN
    DO i= 1, nmat
        WRITE(ounit,*)
        WRITE(ounit,1009) i
        WRITE(ounit,1011)'GROUP', 'TRANSPORT', 'DIFFUSION', 'ABSORPTION', &
        'REMOVAL', 'NU*FISS', 'KAP*FIS','FISS. SPCTR'
        DO g= 1, ng
            WRITE(ounit,1010) g, xsigtr(i,g), xD(i,g), xsiga(i,g), &
            xsigr(i,g), xnuf(i,g), xsigf(i,g), chi(i,g)
        END DO
        WRITE(ounit,*)'  --SCATTERING MATRIX--'
        WRITE(ounit,'(4X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
        DO g= 1, ng
            WRITE(ounit,1015)g, (xsigs(i,g,h), h=1,ng)
        END DO
    END DO
END IF

IF (ccnuf .AND. mode /= 'FIXEDSRC') THEN
    WRITE(ounit, *) "ERROR: The Problem has no fission material (nu*fission for all materials are zero)"
    STOP
END IF
IF (ccsigf .AND. mode /= 'FIXEDSRC') THEN
    WRITE(ounit, *) "ERROR: The Problem has no fission material (fission xsec for all materials are zero)"
    STOP
END IF

WRITE(ounit,*)
WRITE(ounit,*) ' ...Macroscopic CX Card is successfully read...'

1009 FORMAT(5X, 'MATERIAL', I3)
1011 FORMAT(2X, A7, A12, A11, A12, A9, A11, A11, A14)
1010 FORMAT(2X, I6, F13.6, F10.6, F11.6, F12.6, 2F11.6, F10.6)
1015 FORMAT(4X, I3, F16.6, 20F12.6)
1020 FORMAT(2X, 'ERROR: Transport cross section (sigtr)is zero or negative in material: ', I3, ' ;group: ', I3)

DEALLOCATE(group)
DEALLOCATE(xD, xsigr)

END SUBROUTINE inp_xsec

!******************************************************************************!

SUBROUTINE inp_geom1 (xbunit)
!
! Purpose:
!    To read geometry card in input (1st part)
!

USE sdata, ONLY: nx, ny, nz, nxx, nyy, nzz, xdel, ydel, zdel, &
                xdiv, ydiv, zdiv

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln, ios

INTEGER :: i, j, k, lx, ly, lz, xtot, ytot, ztot
REAL(DP) :: div

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>>READING CORE GEOMETRY<<<<<'
WRITE(ounit,*) '           -------------------------------'

! Read number of assemblies in x, y and z directions
READ(xbunit, *, IOSTAT=ios) ind, ln, nx, ny, nz
message = ' error in reading number assemblies'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

! Limit values of nx, ny and nz
IF (nx < 2) THEN
    WRITE(ounit,*) '  Error: nx shall be at least 2'
    STOP
END IF
IF (ny < 2) THEN
    WRITE(ounit,*) '  Error: ny shall be at least 2'
    STOP
END IF
IF (nz < 2) THEN
    WRITE(ounit,*) '  Error: nz shall be at least 2'
    STOP
END IF
IF (nx > 33) THEN
    WRITE(ounit,*) '  Error: nx shall be maximum 33'
    STOP
END IF
IF (ny > 33) THEN
    WRITE(ounit,*) '  Error: ny shall be maximum 33'
    STOP
END IF
IF (nz > 40) THEN
    WRITE(ounit,*) '  Error: nz shall be maximum 40'
    STOP
END IF

ALLOCATE(xsize(nx), ysize(ny), zsize(nz))
ALLOCATE(xdiv(nx), ydiv(ny), zdiv(nz))

! Read assembly sizes and assembly division in x, y and z directions
! x-direction
READ(xbunit, *, IOSTAT=ios) ind, ln, (xsize(i), i=1,nx)
message = ' error in reading x-direction assembly sizes'
CALL er_message(ounit, ios, ln, message, buf=xbunit)
READ(xbunit, *, IOSTAT=ios) ind, ln, (xdiv(i), i=1,nx)
message = ' error in reading x-direction assembly division'
CALL er_message(ounit, ios, ln, message, buf=xbunit)
! y-direction
READ(xbunit, *, IOSTAT=ios) ind, ln, (ysize(j), j=1,ny)
message = ' error in reading y-direction assembly sizes'
CALL er_message(ounit, ios, ln, message, buf=xbunit)
READ(xbunit, *, IOSTAT=ios) ind, ln, (ydiv(j), j=1,ny)
message = ' error in reading y-direction assembly division'
CALL er_message(ounit, ios, ln, message, buf=xbunit)
! z-direction
READ(xbunit, *, IOSTAT=ios) ind, ln, (zsize(k), k=1,nz)
message = ' error in reading z-direction assembly sizes'
CALL er_message(ounit, ios, ln, message, buf=xbunit)
READ(xbunit, *, IOSTAT=ios) ind, ln, (zdiv(k), k=1,nz)
message = ' error in reading z-direction assembly division'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

!Calculate number of nodes in x, y and z directions
!x-direction
nxx=0
DO i= 1,nx
   nxx = nxx+xdiv(i)
END DO
!y-direction
nyy=0
DO j= 1,ny
   nyy = nyy+ydiv(j)
END DO
!z-direction
nzz=0
DO k= 1,nz
   nzz = nzz+zdiv(k)
END DO

! Limit values of nxx, nyy and nzz
IF (nxx > 800) THEN
    WRITE(ounit,*) '  Error: nxx shall be maximum 800'
    WRITE(*,*) '  Error: nxx shall be maximum 800'
    STOP
END IF
IF (nyy > 800) THEN
    WRITE(ounit,*) '  Error: nyy shall be maximum 800'
    WRITE(*,*) '  Error: nyy shall be maximum 800'
    STOP
END IF
IF (nzz > 800) THEN
    WRITE(ounit,*) '  Error: nzz shall be maximum 800'
    WRITE(*,*) '  Error: nzz shall be maximum 800'
    STOP
END IF

ALLOCATE(xdel(nxx), ydel(nyy), zdel(nzz))

!Calculate delta x, y, and z (node sizes)
!Delta x
xtot=0
DO i= 1,nx
    div = xsize(i)/REAL(xdiv(i))
    DO lx= 1, xdiv(i)
    xtot = xtot+1
    xdel(xtot) = div
    END DO
END DO
!Delta y
ytot=0
DO j= 1,ny
    div = ysize(j)/REAL(ydiv(j))
    DO ly= 1, ydiv(j)
    ytot = ytot+1
    ydel(ytot) = div
    END DO
END DO
!Delta z
ztot=0
DO k= 1,nz
    div = zsize(k)/REAL(zdiv(k))
    DO lz= 1, zdiv(k)
    ztot = ztot+1
    zdel(ztot) = div
    END DO
END DO

END SUBROUTINE inp_geom1

!******************************************************************************!

SUBROUTINE inp_geom2 (xbunit)
!
! Purpose:
!    To read geometry card in input (2nd part)
!

USE sdata, ONLY: nx, ny, nz, nxx, nyy, nzz, xdel, ydel, zdel, &
                xwest, xeast, ynorth, ysouth, zbott, ztop, nnod, &
                xstag, ystag, xdiv, ydiv, zdiv, nmat, coreh

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln, ios

INTEGER :: i, j, k, lx, ly, lz, xtot, ytot, ztot
CHARACTER(LEN=2), DIMENSION(nxx, nyy) :: mmap
INTEGER, PARAMETER :: xm = 36
INTEGER :: ip, ipr, kp
INTEGER :: xs, xf

! Reading number of planar
READ(xbunit, *, IOSTAT=ios) ind, ln, np
message = ' error in reading number of planars'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

!Reading planar assignment into z-direction
ALLOCATE(zpln(nz))
READ(xbunit, *, IOSTAT=ios) ind, ln, (zpln(k), k=1,nz)
message = ' error in reading planar assignment'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

DO k = 1, nz
    IF (zpln(k) > np) THEN
        WRITE(ounit,'(2X,A15,I3,A35)') 'ERROR: PLANAR ', &
        zpln(k), ' IS GREATER THAN NUMBER OF PLANAR'
        STOP
    END IF
    IF (zpln(k) < 1) THEN
        WRITE(ounit,'(2X,A)') 'ERROR: PLANAR SHOULD BE AT LEAST 1'
        STOP
    END IF
END DO

! Reading material assignment for each planar
ALLOCATE(plnr(np))
DO k= 1,np
    ALLOCATE(plnr(k)%asm (nx,ny))
    ALLOCATE(plnr(k)%node(nxx,nyy))
    DO j= ny, 1, -1
        READ(xbunit, *, IOSTAT=ios) ind, ln, (plnr(k)%asm(i,j), i=1,nx)
        message = ' error in reading planar'
        CALL er_message(ounit, ios, ln, message, buf=xbunit)
    END DO
END DO

! Material assignment into nodes (not part for geom output)
ALLOCATE (mnum(nxx, nyy, nzz))

ztot = 0
DO k= 1, nz
    DO lz= 1, zdiv(k)
        ztot = ztot+1
        ytot = 0
        DO j= 1, ny
            DO ly= 1, ydiv(j)
                ytot = ytot+1
                xtot = 0
                DO i= 1, nx
                    DO lx= 1, xdiv(i)
                        xtot = xtot+1
                        mnum(xtot, ytot, ztot) = plnr(zpln(k))%asm(i,j)
                        IF (mnum(xtot, ytot, ztot) > nmat) THEN
                            WRITE(ounit,'(2X,A17,I3,A37)') 'ERROR: MATERIAL ', &
                            mnum(xtot, ytot, ztot), ' IS GREATER THAN NUMBER OF MATERIAL'
                            WRITE(*,'(2X,A17,I3,A37)') 'ERROR: MATERIAL ', &
                            mnum(xtot, ytot, ztot), ' IS GREATER THAN NUMBER OF MATERIAL'
                            STOP
                        END IF
                        IF (mnum(xtot, ytot, ztot) < 0) THEN
                            WRITE(ounit,'(2X,A)') 'ERROR: NEGATIVE MATERIAL FOUND'
                            WRITE(*,'(2X,A)') 'ERROR: NEGATIVE MATERIAL FOUND'
                            STOP
                        END IF
                    END DO
                END DO
            END DO
        END DO
    END DO
END DO

! Assign assembly wise planar to node wise planar for output
ztot = 0
DO k= 1, np
  ytot = 0
  DO j= 1, ny
      DO ly= 1, ydiv(j)
          ytot = ytot+1
          xtot = 0
          DO i= 1, nx
              DO lx= 1, xdiv(i)
                  xtot = xtot+1
                  plnr(k)%node(xtot,ytot) = plnr(k)%asm(i,j)
              END DO
          END DO
      END DO
  END DO
END DO

! -Indexing non zero material for staggered mesh-
ALLOCATE(ystag(nyy), xstag(nxx))
!Indexing non zero material for staggered mesh along y direction
DO j= 1, nyy
    ystag(j)%smin = nxx
    DO i = 1, nxx
        IF (mnum(i,j,1) /= 0) THEN
            ystag(j)%smin = i
            EXIT
        END IF
    END DO
END DO

DO j= 1, nyy
    ystag(j)%smax = 0
    DO i = nxx, 1, -1
        IF (mnum(i,j,1) /= 0) THEN
            ystag(j)%smax = i
            EXIT
        END IF
    END DO
END DO

!Indexing non zero material for staggered mesh along x direction
DO i= 1, nxx
    xstag(i)%smin = nyy
    DO j = 1, nyy
        IF (mnum(i,j,1) /= 0) THEN
            xstag(i)%smin = j
            EXIT
        END IF
    END DO
END DO

DO i= 1, nxx
    xstag(i)%smax = 0
    DO j = nyy, 1, -1
        IF (mnum(i,j,1) /= 0) THEN
            xstag(i)%smax = j
            EXIT
        END IF
    END DO
END DO

! Checking zero material between non-zero material
DO k = 1, nzz
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax
            IF (mnum(i,j,k) == 0) THEN
                WRITE(ounit,*) 'Zero material found inside core. Check material assignment'
                STOP
            END IF
        END DO
    END DO
END DO
DO k = 1, nzz
    DO i = 1, nxx
        DO j = xstag(i)%smin, xstag(i)%smax
            IF (mnum(i,j,k) == 0) THEN
                WRITE(ounit,*) 'Zero material found inside core. Check material assignment'
                STOP
            END IF
        END DO
    END DO
END DO

!Reading Boundary Conditions
READ(xbunit, *, IOSTAT=ios) ind, ln, xeast, xwest, ynorth, ysouth, zbott, ztop
message = ' error in reading boundary conditions'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

IF (xeast > 2 .OR. xwest > 2 .OR. ynorth > 2 .OR. ysouth > 2 &
.OR. zbott > 2 .OR. ztop > 2) then
  write(*,1019) ln
  stop
END IF


! Wrting core geometry output
IF (ogeom) THEN
    WRITE(ounit,*)' Number of assembly in x, y and z directions respectively :'
    WRITE(ounit,*) nx, ny, nz
    WRITE(ounit,*)' Number of nodes in x, y and z directions respectively :'
    WRITE(ounit,*) nxx, nyy, nzz
    WRITE(ounit,*)
    WRITE(ounit,1016) 'x','x'
    WRITE(ounit,'(2X,10F7.2)')(xdel(i), i=1,nxx)
    WRITE(ounit,1016) 'y','y'
    WRITE(ounit,'(2X,10F7.2)')(ydel(j), j=1,nyy)
    WRITE(ounit,*)

    IF (nxx < 100) THEN
      ip = nxx/xm
      ipr = MOD(nxx,xm) - 1
      DO k= 1,np

          DO j = 1, nyy
              DO i = 1, nxx
                  IF (plnr(k)%node(i,j) == 0) THEN
                      mmap(i,j) = '  '
                  ELSE
                      WRITE (mmap(i,j),'(I2)') plnr(k)%node(i,j)
                      mmap(i,j) = TRIM(ADJUSTL(mmap(i,j)))
                  END IF
              END DO
          END DO

          WRITE(ounit,1017) k
          xs = 1; xf = xm
          DO kp = 1, ip
              WRITE(ounit,'(6X,100I3)') (i, i = xs, xf)
              DO j= nyy, 1, -1
                  WRITE(ounit,'(2X,I4,1X,100A3)') j, (mmap(i,j), i=xs, xf)
              END DO
              xs = xs + xm
              xf = xf + xm
          END DO

          WRITE(ounit,'(6X,100I3)') (i, i = xs, xs+ipr)
          IF (xs+ipr > xs) THEN
              DO j= nyy, 1, -1
                  WRITE(ounit,'(2X,I4,1X,100A3)') j, (mmap(i,j), i=xs, xs+ipr)
              END DO
          END IF
      END DO
    END IF

    WRITE(ounit,*)
    WRITE(ounit,1018)
    WRITE(ounit,*) '--------------------------------------'
    WRITE(ounit,*) '  Plane Number     Planar Region    delta-z'
    ztot = nzz
    DO k= nz, 1, -1
        DO lz= 1, zdiv(k)
            IF (ztot == nzz) THEN
                WRITE(ounit,'(I9, A6, I13, F15.2)') ztot, ' (TOP)', zpln(k), zdel(ztot)
            ELSE IF (ztot == 1) THEN
                 WRITE(ounit,'(I9, A9, I10, F15.2)') ztot, ' (BOTTOM)', zpln(k), zdel(ztot)
            ELSE
                WRITE(ounit,'(I9, I19, F15.2)') ztot, zpln(k), zdel(ztot)
            END IF
            ztot = ztot - 1
        END DO
    END DO
    WRITE(ounit,*)
    WRITE(ounit,*) '  Boundary conditions'

    IF (xwest == 0) THEN
        WRITE(ounit,*)' X-directed West   : ZERO FLUX'
    ELSE IF (xwest == 1) THEN
        WRITE(ounit,*)' X-directed West   : ZERO INCOMING CURRENT'
    ELSE
        WRITE(ounit,*)' X-directed West   : REFLECTIVE'
    END IF

    IF (xeast == 0) THEN
        WRITE(ounit,*)' X-directed East   : ZERO FLUX'
    ELSE IF (xeast == 1) THEN
        WRITE(ounit,*)' X-directed East   : ZERO INCOMING CURRENT'
    ELSE
        WRITE(ounit,*)' X-directed East   : REFLECTIVE'
    END IF

    IF (ynorth == 0) THEN
        WRITE(ounit,*)' Y-directed North  : ZERO FLUX'
    ELSE IF (ynorth == 1) THEN
        WRITE(ounit,*)' Y-directed North  : ZERO INCOMING CURRENT'
    ELSE
        WRITE(ounit,*)' Y-directed North  : REFLECTIVE'
    END IF

    IF (ysouth == 0) THEN
        WRITE(ounit,*)' Y-directed South  : ZERO FLUX'
    ELSE IF (ysouth == 1) THEN
        WRITE(ounit,*)' Y-directed South  : ZERO INCOMING CURRENT'
    ELSE
        WRITE(ounit,*)' Y-directed South  : REFLECTIVE'
    END IF

    IF (zbott == 0) THEN
        WRITE(ounit,*)' Z-directed Bottom : ZERO FLUX'
    ELSE IF (zbott == 1) THEN
        WRITE(ounit,*)' Z-directed Bottom : ZERO INCOMING CURRENT'
    ELSE
        WRITE(ounit,*)' Z-directed Bottom : REFLECTIVE'
    END IF

    IF (ztop == 0) THEN
        WRITE(ounit,*)' Z-directed Top    : ZERO FLUX'
    ELSE IF (ztop == 1) THEN
        WRITE(ounit,*)' Z-directed Top    : ZERO INCOMING CURRENT'
    ELSE
        WRITE(ounit,*)' Z-directed Top    : REFLECTIVE'
    END IF
END IF


1016 FORMAT(2X,A,'-directed nodes division (delta-',A,')')
1017 FORMAT(3X, 'Material Map for Planar Region : ', I2)
1018 FORMAT(2X, 'Planar Region Assignment to planes.')
1019 FORMAT(2X, 'ERROR: Wrong boundary conditions in line : ', I4)

WRITE(ounit,*)
WRITE(ounit,*) ' ...Core geometry is successfully read...'

! Calculate core height
coreh = 0._DP
DO k = 1, nzz
    coreh = coreh + zdel(k)
END DO

! set Number of nodes (nnod)
nnod = 0
DO k = 1, nzz
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax
            nnod = nnod + 1
        END DO
    END DO
END DO

END SUBROUTINE inp_geom2

!******************************************************************************!

SUBROUTINE misc ()
!
! Purpose:
!    To assign material xsec to nodes
!    To arranges the nodes into 1-D array instead of 3-D array
!    To allocate CMFD matrix index
!    etc


USE sdata, ONLY: nxx, nyy, nzz, ix, iy, iz, xyz, &
                 nnod, sigtr, siga, nuf, sigf, &
                 sigs, D, sigr, ng, ystag, xstag, &
                 xdel, ydel, zdel, vdel, nupd, &
                 mat, ind

IMPLICIT NONE

INTEGER :: i, j, k, n, noz

ALLOCATE(ix(nnod), iy(nnod), iz(nnod))
ALLOCATE(xyz(nxx, nyy, nzz))
ALLOCATE(mat(nnod))

! Set ix, iy, iz and xyz
n = 0
xyz = 0
DO k = 1, nzz      ! VERY IMPORTANT: NEVER CHANGE THE ORDERS OF THESE NESTED LOOPS
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax
             n = n + 1
             ix(n) = i
             iy(n) = j
             iz(n) = k
             xyz(i,j,k) = n
             mat(n) = mnum(i,j,k)
        END DO
    END DO
END DO

ALLOCATE(sigtr(nnod,ng))
ALLOCATE(siga (nnod,ng))
ALLOCATE(nuf  (nnod,ng))
ALLOCATE(sigf (nnod,ng))
ALLOCATE(sigs (nnod,ng,ng))
ALLOCATE(D    (nnod,ng))
ALLOCATE(sigr (nnod,ng))

! Calculate nodes' volume
ALLOCATE(vdel(nnod))
DO i = 1, nnod
    vdel(i) = xdel(ix(i)) * ydel(iy(i)) * zdel(iz(i))
END DO

! Calculate number of non-zero element in the CMFD Matrix
ALLOCATE(ind(nnod))
do n = 1, nnod
  noz = 0
  i = ix(n); j = iy(n); k = iz(n)         ! Set i, j, k
  if (k /= 1) noz = noz + 1
  if (j /= xstag(i)%smin) noz = noz + 1
  if (i /= ystag(j)%smin) noz = noz + 1
  noz = noz + 1
  if (i /= ystag(j)%smax) noz = noz + 1
  if (j /= xstag(i)%smax) noz = noz + 1
  if (k /= nzz) noz = noz + 1
  ind(n)%ncol = noz
  ALLOCATE(ind(n)%col(noz))
end do

! calaculate default nodal update interval
nupd = ceiling((nxx + nyy + nzz) / 2.5)

END SUBROUTINE misc

!******************************************************************************!

SUBROUTINE inp_esrc (xbunit)
!
! Purpose:
!    To read extra sources if any
!

USE sdata, ONLY: exsrc, ng, nx, ny, nz, &
                 xdiv, ydiv, zdiv, xyz

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln, ios  ! line number and IOSTAT status
INTEGER :: g, i, j, k, n
INTEGER :: xt, yt, zt
INTEGER :: it, jt, kt
INTEGER :: nsrc
REAL(DP)    :: sden                                       ! Source density
REAL(DP), DIMENSION(:), ALLOCATABLE :: spec               ! Source Spectrum
CHARACTER(LEN=1), DIMENSION(:,:), ALLOCATABLE :: spos ! Source position
INTEGER :: xpos, ypos, zpos
REAL(DP) :: summ

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>>READING EXTRA SOURCE<<<<<'
WRITE(ounit,*) '           -------------------------------'

! Read number of source
READ(xbunit, *, IOSTAT=ios) ind, ln, nsrc
message = ' error in reading number of extra source'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

ALLOCATE(spec(ng), spos(nx,ny))

DO n = 1, nsrc
    ! Read source density
    READ(xbunit, *, IOSTAT=ios) ind, ln, sden
    message = ' error in reading source density'
    CALL er_message(ounit, ios, ln, message, buf=xbunit)

    IF (sden <= 0.0) THEN
        WRITE(ounit,*) '  ERROR: SOURCE DENSITY SHALL BE GREATER THAN ZERO'
        STOP
    END IF

    ! Read source spectrum
    READ(xbunit, *, IOSTAT=ios) ind, ln, (spec(g), g = 1, ng)
    message = ' error in reading source spectrum'
    CALL er_message(ounit, ios, ln, message, buf=xbunit)

    ! Is total spectrum = 1._DP?
    summ = 0._DP
    DO g = 1, ng
        summ = summ + spec(g)
    END DO
    ! Check total spectrum
    IF (ABS(summ - 1._DP) > 1.e-5_DP) THEN
        WRITE(ounit,*) 'TOTAL SOURCE SPECTRUM AT LINE', ln, ' IS NOT EQUAL TO 1._DP'
        STOP
    END IF

    ! WRITE OUTPUT
    WRITE(ounit,'(A12,I3)') '     SOURCE ', n
    WRITE(ounit,*)         '-----------------'
    WRITE(ounit,'(A20,ES10.3, A11)') '  Source Density  : ', sden, '  n/(m^3*s)'
    WRITE(ounit,'(A19,100F6.2)') '  Source Spectrum : ', (spec(g), g = 1, ng)
    WRITE(ounit,*) ' Source Position '

    ! Read source position
    DO
        READ(xbunit, *, IOSTAT=ios) ind, ln, zpos
        message = ' error in reading axial position (zpos) of extra source'
        CALL er_message(ounit, ios, ln, message, buf=xbunit)
        IF (zpos < 1) EXIT
        IF (zpos > nz) THEN
            WRITE(ounit,* ) '  ERROR: WRONG EXTRA SOURCES POSITION (ZPOS)'
            WRITE(ounit, 2033) ln, zpos
            STOP
        END IF
        spos = '0'
        DO
            READ(xbunit, *, IOSTAT=ios) ind, ln, xpos, ypos
            message = ' error in reading radial position (xpos and ypos) of extra source'
            CALL er_message(ounit, ios, ln, message, buf=xbunit)
            IF (xpos < 1 .OR. ypos < 1) EXIT

            IF (xpos > nx) THEN
                WRITE(ounit,* ) '  ERROR: WRONG EXTRA SOURCES POSITION (XPOS)'
                WRITE(ounit, 2033) ln, xpos, ypos
                STOP
            END IF

            IF (ypos > ny) THEN
                WRITE(ounit,* ) '  ERROR: WRONG EXTRA SOURCES POSITION (YPOS)'
                WRITE(ounit, 2033) ln, xpos, ypos
                STOP
            END IF

            spos(xpos,ypos) = 'X'

            zt = 0
            kt = 1
            DO k = 1, zpos
                IF (k > 1) kt = zt + 1
                zt = zt + zdiv(k)
            END DO

            yt = 0
            jt = 1
            DO j = 1, ypos
                IF (j > 1) jt = yt + 1
                yt = yt + ydiv(j)
            END DO

            xt = 0
            it = 1
            DO i = 1, xpos
                IF (i > 1) it = xt + 1
                xt = xt + xdiv(i)
            END DO

            DO k = kt, zt
                DO j = jt, yt
                    DO i = it, xt
                        DO g = 1, ng
                            exsrc(xyz(i,j,k), g) = exsrc(xyz(i,j,k), g) + &
                            sden * spec(g)
                        END DO
                    END DO
                END DO
            END DO

        END DO

        WRITE(ounit,'(A18,I3)') '   Plane number : ', zpos
        WRITE(ounit,'(7X,100I3)') (i, i = 1, nx)
        DO j = ny, 1, -1
            WRITE(ounit,'(4X,I3, 100A3 )') j, (spos(i,j), i=1, nx)
        END DO
        WRITE(ounit,*)
    END DO
END DO

2033 FORMAT(2X,'LINE', I4, ' : ', I3, I3)

DEALLOCATE(spec, spos)

END SUBROUTINE inp_esrc

!******************************************************************************!

SUBROUTINE inp_iter (xbunit)

!
! Purpose:
!    To read iteration control if any

USE sdata, ONLY: nout, nin, serc, ferc, nac, nupd, th_niter, nth, kern

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

READ(xbunit, *, IOSTAT=ios) ind, ln, nout, nin, serc, ferc, nac, nupd, &
th_niter, nth
message = ' error in reading iteration control'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>>READING ITERATION CONTROL<<<<<'
WRITE(ounit,*) '           ------------------------------------'

WRITE(ounit,'(A,I5)') '  MAXIMUM NUMBER OF OUTER ITERATION                     : ', nout
WRITE(ounit,'(A,I5)') '  MAXIMUM NUMBER OF INNER ITERATION                     : ', nin
WRITE(ounit,'(A,ES12.3)') '  FISSION SOURCE ERROR CRITERIA                     : ', serc
WRITE(ounit,'(A,ES12.3)') '  FLUX ERROR CRITERIA                               : ', ferc
WRITE(ounit,'(A,I5)') '  OUTER ITERATION FISSION SOURCE EXTRAPOLATION INTERVAL : ', nac
WRITE(ounit,'(A,I5)') '  NODAL UPDATE INTERVAL                                 : ', nupd
WRITE(ounit,'(A,I5)') '  MAX. NUMBER OF T-H ITERATION                          : ', th_niter
WRITE(ounit,'(A,I5)') '  MAX. NUMBER OF OUTER ITERATION PER T-H ITERATION      : ', nth

if (nout < nupd .and. kern /= ' FDM') then
  write(*,*) "ERROR: MAX. NUMBER OF OUTER ITERATION SHOULD BE BIGGER THAN NODAL UPDATE INTERVAL"
  write(ounit,*) "ERROR: MAX. NUMBER OF OUTER ITERATION SHOULD BE BIGGER THAN NODAL UPDATE INTERVAL"
  stop
end if

if (nth < nupd .and. kern /= ' FDM' .and. bther == 1) then
  write(*,*) "ERROR: MAX. NUMBER OF OUTER ITERATION PER T-H ITERATION SHOULD BE BIGGER THAN NODAL UPDATE INTERVAL"
  write(ounit,*) "ERROR: MAX. NUMBER OF OUTER ITERATION PER T-H ITERATION SHOULD BE BIGGER THAN NODAL UPDATE INTERVAL"
  stop
end if


END SUBROUTINE inp_iter

!******************************************************************************!

SUBROUTINE inp_extr ()

!
! Purpose:
!    To tell user the flux transformation is active

IMPLICIT NONE

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>>READING FLUX TRANSFORMATION<<<<<'
WRITE(ounit,*) '           -------------------------------------'
WRITE(ounit,1571)

1571 format (2X, 'EXPONENTIAL FLUX TRANSFORMATION IS ACTIVE')

END SUBROUTINE inp_extr

!******************************************************************************!

SUBROUTINE inp_kern (xbunit)

!
! Purpose:
!    To determine nodal kernel

USE sdata, ONLY:  kern

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status
CHARACTER (LEN=30) :: kern_desc


READ(xbunit, *, IOSTAT=ios) ind, ln, kern
message = ' error in reading iteration control'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

kern = TRIM(ADJUSTR(kern))
if (kern /= ' FDM' .and. kern /= ' PNM' .and. kern /= 'SANM') then
  write(*,1646) kern
  write(*,*) ' NODAL KERNEL OPTIONS: FDM, PNM OR SANM'
  write(ounit,1646) kern
  write(ounit,*) ' NODAL KERNEL OPTIONS: FDM, PNM OR SANM'
  stop
end if

if (kern == ' FDM') then
  kern_desc = "FINITE DIFFERENCE METHOD"
else if (kern == ' PNM') then
  kern_desc = "POLYNOMIAL NODAL METHOD"
else
  kern_desc = "SEMI-ANALYTIC NODAL METHOD"
end if
kern_desc = TRIM(ADJUSTL(kern_desc))

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '             >>>>>READING NODAL KERNEL<<<<<'
WRITE(ounit,*) '           ------------------------------------'

WRITE(ounit,1648) kern_desc
if (scr) then
  WRITE(*,*)
  WRITE(*,1648) kern_desc
end if

1646 format (2X, 'ERROR: COULD NOT RECOGNIZE NODAL KERNEL: ', A4)
1648 format (2X, 'NODAL KERNEL  : ', A30)

END SUBROUTINE inp_kern

!******************************************************************************!

SUBROUTINE inp_thet (xbunit)

!
! Purpose:
!    To read iteration theta value for transient problem

USE sdata, ONLY: sth, bth

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

READ(xbunit, *, IOSTAT=ios) ind, ln, sth
message = ' error in theta in %THETA card'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

if (sth < 0.001) then
  write(*,*)
  write(*,*) "ERROR: THETA VALUE IS TOO SMALL"
  stop
end if
if (sth > 1.0) then
  write(*,*)
  write(*,*) "ERROR: THETA VALUE SHALL NOT GREATER THAT 1.0"
  stop
end if
bth = (1._dp - sth) / sth

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '               >>>>>READING THETA CARD<<<<<'
WRITE(ounit,*) '           ------------------------------------'

WRITE(ounit,'(A,F6.2)') '  THETA IS  : ', sth


END SUBROUTINE inp_thet

!******************************************************************************!

SUBROUTINE inp_prnt (xbunit)

!
! Purpose:
!    To read output print option if any

USE sdata, ONLY: aprad, apaxi, afrad

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

CHARACTER(LEN=3) :: caprad='YES', capaxi='YES', cafrad='YES'

READ(xbunit, *, IOSTAT=ios) ind, ln, aprad, apaxi, afrad
message = ' error in reading output print option'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

IF (aprad == 0) caprad='NO'
IF (apaxi == 0) capaxi='NO'
IF (afrad == 0) cafrad='NO'

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>>READING OUTPUT PRINT OPTION<<<<<'
WRITE(ounit,*) '           -------------------------------------'

WRITE(ounit,'(A,A)') '  RADIAL ASSEMBLY POWER DISTRIBUTION : ', caprad
WRITE(ounit,'(A,A)') '  AXIAL ASSEMBLY POWER DISTRIBUTION  : ', capaxi
WRITE(ounit,'(A,A)') '  RADIAL FLUX POWER DISTRIBUTION     : ', cafrad


END SUBROUTINE inp_prnt

!******************************************************************************!

SUBROUTINE inp_adf (xbunit)

!
! Purpose:
!    To read ADF values if any

USE sdata, ONLY: ng, nmat, nx, ny, nz, nxx, nyy, nzz, &
                 xdiv, ydiv, zdiv, xyz, ystag, dc

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

TYPE :: ADF_TYPE
    REAL(DP), DIMENSION(6) :: dc
END TYPE
TYPE(ADF_TYPE), DIMENSION(nmat,ng) :: mdc
TYPE(ADF_TYPE), DIMENSION(nx,ny,nz,ng) :: xdc
TYPE(ADF_TYPE), DIMENSION(nxx,nyy,nzz,ng) :: xxdc

INTEGER :: g, i, j, k, u
INTEGER :: rot, x1, x2, y1, y2, z1, z2
INTEGER :: xtot, ytot, ztot
INTEGER :: lz, ly, lx
INTEGER, DIMENSION(nx+1) :: tx
INTEGER, DIMENSION(ny+1) :: ty
INTEGER, DIMENSION(nz+1) :: tz
INTEGER :: zp
CHARACTER(LEN=6), DIMENSION(nx, ny) :: cadf
INTEGER, PARAMETER :: xm = 12
INTEGER :: ip, ipr
INTEGER :: xs, xf

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '         >>>>>READING ASSEMBLY DISCONTINUITY FACTOR<<<<<'
WRITE(ounit,*) '         -----------------------------------------------'

DO i = 1, nmat
    DO g = 1, ng
        READ(xbunit, *, IOSTAT=ios) ind, ln, (mdc(i,g)%dc(j), j = 1, 6)
        message = ' error in reading Assembly Discontinuity Factor (ADF)'
        CALL er_message(ounit, ios, ln, message, buf=xbunit)
    END DO
END DO

DO g = 1, ng
    DO k = 1, nz
        DO j = 1, ny
            DO i = 1, nx
                IF (plnr(zpln(k))%asm(i,j) /= 0) xdc(i,j,k,g) = mdc(plnr(zpln(k))%asm(i,j),g)
            END DO
        END DO
    END DO
END DO

!!! ADF ROTATION
DO
    READ(xbunit, *, IOSTAT=ios) ind, ln, rot
    message = ' error in reading ADF Rotation'
    CALL er_message(ounit, ios, ln, message, buf=xbunit)
    IF (rot < 1) EXIT
    IF (rot > 3) THEN
        WRITE(ounit,*) '  ERROR: MAXIMUM ADF ROTATION IS 3 TIMES'
        WRITE(ounit,2030) ln, rot
    END IF

    DO
        READ(xbunit, *, IOSTAT=ios) ind, ln, x1, x2, y1, y2, z1, z2
        message = ' error in reading ADF Rotation'
        CALL er_message(ounit, ios, ln, message, buf=xbunit)
        IF (x1 < 1 .OR. x2 < 1 .OR. y1 < 1 .OR. y2 < 1 .OR. z1 < 1 .OR. z2 < 1) EXIT
        IF (x1 > nx .OR. x2 > nx .OR. y1 > ny .OR. y2 > ny .OR. z1 > nz .OR. z2 > nz) THEN
            WRITE(ounit,*) '  ERROR: WRONG POSITION FOR ADF ROTATION. ' // &
                           'OUT OF DIMENSION OF THE CORE'
            WRITE(ounit,2032) ln, x1, x2, y1, y2, z1, z2
            STOP
        END IF
        IF (x2 < x1 .OR. y2 < y1 .OR. z2 < z1) THEN
            WRITE(ounit,*) '  ERROR: WRONG POSITION FOR ADF ROTATION. ' // &
                           'INITIAL POSITION IS SMALLER'
            WRITE(ounit,2032) ln, x1, x2, y1, y2, z1, z2
        END IF

        DO g = 1, ng
            DO k = z1, z2
                DO j = y1, y2
                    DO i = x1, x2
                        CALL rotate(rot, xdc(i,j,k,g)%dc(1), xdc(i,j,k,g)%dc(2), &
                                         xdc(i,j,k,g)%dc(3), xdc(i,j,k,g)%dc(4))
                    END DO
                END DO
            END DO
        END DO

    END DO
END DO

! ADF PRINT OPTION
READ(xbunit, *, IOSTAT=ios) ind, ln, zp
IF (ios == 0 .AND. zp >=1) THEN
    WRITE(ounit,*)
    WRITE(ounit,'(A,I3)') '  ADF VALUES ON PLANAR NUMBER : ', zp

    ip = nx/xm
    ipr = MOD(nx,xm) - 1
    DO g = 1, ng
        WRITE(ounit,*)
        WRITE(ounit, 1999) g
        WRITE(ounit,*)

        WRITE(ounit,*) '  EAST ADF'
        DO j = 1, ny
            DO i = 1, nx
                !!! If ADF > 0, Convert ADF to character
                IF ((xdc(i,j,zp,g)%dc(1) - 0.) < 1.e-5_DP) THEN
                    cadf(i,j) = '      '
                ELSE
                    WRITE (cadf(i,j),'(F6.4)') xdc(i,j,zp,g)%dc(1)
                    cadf(i,j) = TRIM(ADJUSTL(cadf(i,j)))
                END IF
            END DO
        END DO

        xs = 1; xf = xm
        DO k = 1, ip
            WRITE(ounit,'(4X,100I8)') (i, i = xs, xf)
            DO j= ny, 1, -1
                WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xf)
            END DO
            WRITE(ounit,*)
            xs = xs + xm
            xf = xf + xm
        END DO

        WRITE(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
        IF (xs+ipr > xs) THEN
            DO j= ny, 1, -1
                WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xs+ipr)
            END DO
        END IF


        WRITE(ounit,*) '  WEST ADF'
        DO j = 1, ny
            DO i = 1, nx
                !!! If ADF > 0, Convert ADF to character
                IF ((xdc(i,j,zp,g)%dc(2) - 0.) < 1.e-5_DP)  THEN
                    cadf(i,j) = '      '
                ELSE
                    WRITE (cadf(i,j),'(F6.4)') xdc(i,j,zp,g)%dc(2)
                    cadf(i,j) = TRIM(ADJUSTL(cadf(i,j)))
                END IF
            END DO
        END DO

        xs = 1; xf = xm
        DO k = 1, ip
            WRITE(ounit,'(4X,100I8)') (i, i = xs, xf)
            DO j= ny, 1, -1
                WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xf)
            END DO
            WRITE(ounit,*)
            xs = xs + xm
            xf = xf + xm
        END DO

        WRITE(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
        IF (xs+ipr > xs) THEN
            DO j= ny, 1, -1
                WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xs+ipr)
            END DO
        END IF


        WRITE(ounit,*) '  NORTH ADF'
        DO j = 1, ny
            DO i = 1, nx
                !!! If ADF > 0, Convert ADF to character
                IF ((xdc(i,j,zp,g)%dc(3) - 0.) < 1.e-5_DP) THEN
                    cadf(i,j) = '      '
                ELSE
                    WRITE (cadf(i,j),'(F6.4)') xdc(i,j,zp,g)%dc(3)
                    cadf(i,j) = TRIM(ADJUSTL(cadf(i,j)))
                END IF
            END DO
        END DO

        xs = 1; xf = xm
        DO k = 1, ip
            WRITE(ounit,'(4X,100I8)') (i, i = xs, xf)
            DO j= ny, 1, -1
                WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xf)
            END DO
            WRITE(ounit,*)
            xs = xs + xm
            xf = xf + xm
        END DO

        WRITE(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
        IF (xs+ipr > xs) THEN
            DO j= ny, 1, -1
                WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xs+ipr)
            END DO
        END IF


        WRITE(ounit,*) '  SOUTH ADF'
        DO j = 1, ny
            DO i = 1, nx
                !!! If ADF > 0, Convert ADF to character
                IF ((xdc(i,j,zp,g)%dc(4) - 0.) < 1.e-5_DP) THEN
                    cadf(i,j) = '      '
                ELSE
                    WRITE (cadf(i,j),'(F6.4)') xdc(i,j,zp,g)%dc(4)
                    cadf(i,j) = TRIM(ADJUSTL(cadf(i,j)))
                END IF
            END DO
        END DO

        xs = 1; xf = xm
        DO k = 1, ip
            WRITE(ounit,'(4X,100I8)') (i, i = xs, xf)
            DO j= ny, 1, -1
                WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xf)
            END DO
            WRITE(ounit,*)
            xs = xs + xm
            xf = xf + xm
        END DO

        WRITE(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
        IF (xs+ipr > xs) THEN
            DO j= ny, 1, -1
                WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xs+ipr)
            END DO
        END IF


        WRITE(ounit,*) '  TOP ADF'
        DO j = 1, ny
            DO i = 1, nx
                !!! If ADF > 0, Convert ADF to character
                IF ((xdc(i,j,zp,g)%dc(5) - 0.) < 1.e-5_DP) THEN
                    cadf(i,j) = '      '
                ELSE
                    WRITE (cadf(i,j),'(F6.4)') xdc(i,j,zp,g)%dc(5)
                    cadf(i,j) = TRIM(ADJUSTL(cadf(i,j)))
                END IF
            END DO
        END DO

        xs = 1; xf = xm
        DO k = 1, ip
            WRITE(ounit,'(4X,100I8)') (i, i = xs, xf)
            DO j= ny, 1, -1
                WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xf)
            END DO
            WRITE(ounit,*)
            xs = xs + xm
            xf = xf + xm
        END DO

        WRITE(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
        IF (xs+ipr > xs) THEN
            DO j= ny, 1, -1
                WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xs+ipr)
            END DO
        END IF


        WRITE(ounit,*) '  BOTTOM ADF'
        DO j = 1, ny
            DO i = 1, nx
                !!! If ADF > 0, Convert ADF to character
                IF ((xdc(i,j,zp,g)%dc(6) - 0.) < 1.e-5_DP) THEN
                    cadf(i,j) = '      '
                ELSE
                    WRITE (cadf(i,j),'(F6.4)') xdc(i,j,zp,g)%dc(6)
                    cadf(i,j) = TRIM(ADJUSTL(cadf(i,j)))
                END IF
            END DO
        END DO

        xs = 1; xf = xm
        DO k = 1, ip
            WRITE(ounit,'(4X,100I8)') (i, i = xs, xf)
            DO j= ny, 1, -1
                WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xf)
            END DO
            WRITE(ounit,*)
            xs = xs + xm
            xf = xf + xm
        END DO

        WRITE(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
        IF (xs+ipr > xs) THEN
            DO j= ny, 1, -1
                WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xs+ipr)
            END DO
        END IF
        WRITE(ounit,*)
    END DO
END IF

1999 FORMAT (4X, 'GROUP : ', I3)

WRITE(ounit,*) '  ...Assembly Discontinuity Factors are successfully read...'


tx(1) = 1
DO i = 2, nx+1
    tx(i) = tx(i-1) + xdiv(i-1)
END DO
ty(1) = 1
DO j = 2, ny+1
    ty(j) = ty(j-1) + ydiv(j-1)
END DO
tz(1) = 1
DO k = 2, nz+1
    tz(k) = tz(k-1) + zdiv(k-1)
END DO

DO g = 1, ng
    ztot = 0
    DO k= 1, nz
        DO lz= 1, zdiv(k)
            ztot = ztot+1
            ytot = 0
            DO j= 1, ny
                DO ly= 1, ydiv(j)
                    ytot = ytot+1
                    xtot = 0
                    DO i= 1, nx
                        DO lx= 1, xdiv(i)
                            xtot = xtot+1
                            xxdc(xtot,ytot,ztot,g)%dc = 0._DP
                            IF (mnum(xtot, ytot, ztot) /= 0) xxdc(xtot,ytot,ztot,g)%dc = 1._DP
                            IF (xtot == tx(i))     xxdc(xtot,ytot,ztot,g)%dc(2) = xdc(i,j,k,g)%dc(2)
                            IF (xtot == tx(i+1)-1) xxdc(xtot,ytot,ztot,g)%dc(1) = xdc(i,j,k,g)%dc(1)
                            IF (ytot == ty(j))     xxdc(xtot,ytot,ztot,g)%dc(4) = xdc(i,j,k,g)%dc(4)
                            IF (ytot == ty(j+1)-1) xxdc(xtot,ytot,ztot,g)%dc(3) = xdc(i,j,k,g)%dc(3)
                            IF (ztot == tz(k))     xxdc(xtot,ytot,ztot,g)%dc(5) = xdc(i,j,k,g)%dc(5)
                            IF (ztot == tz(k+1)-1) xxdc(xtot,ytot,ztot,g)%dc(6) = xdc(i,j,k,g)%dc(6)
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO
END DO

DO g = 1, ng
    DO k = 1, nzz
        DO j = 1, nyy
            DO i = ystag(j)%smin, ystag(j)%smax
               DO u = 1, 6
                 dc(xyz(i,j,k),g,u) = xxdc(i,j,k,g)%dc(u)
               END DO
            END DO
        END DO
    END DO
END DO


2030 FORMAT(3X,'LINE', I4, ' : ', I3)
2032 FORMAT(3X,'LINE', I4, ' : ', I3, I3, I3, I3, I3, I3)


END SUBROUTINE inp_adf

!******************************************************************************!

SUBROUTINE rotate(rot, a1, a2, a3, a4)

! Purpose:
!           To rotate ADF values (necessary for BWR assemblies)

INTEGER, INTENT(IN) :: rot
REAL(DP), INTENT(INOUT) :: a1, a2, a3, a4
REAL(DP) :: x1, x2, x3, x4


x1 = a1
x2 = a2
x3 = a3
x4 = a4

IF (rot == 1) THEN
    a1 = x4
    a2 = x3
    a3 = x1
    a4 = x2
END IF

IF (rot == 2) THEN
    a1 = x2
    a2 = x1
    a3 = x4
    a4 = x3
END IF

IF (rot == 3) THEN
    a1 = x3
    a2 = x4
    a3 = x2
    a4 = x1
END IF

END SUBROUTINE rotate

!******************************************************************************!

SUBROUTINE inp_crod (xbunit)

!
! Purpose:
!    To read control rod position

USE sdata, ONLY: nx, ny, nmat, ng, xdiv, ydiv, &
                 nxx, nyy, bpos, nb, nstep, pos0, ssize, &
                 dsigtr, dsiga, dnuf, dsigf, dsigs, ddc, nnod, coreh, fbmap

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios

INTEGER :: i, j, g, h
INTEGER, DIMENSION(nx, ny) :: bmap       ! Radial control rod bank map (assembly wise)
INTEGER :: popt
INTEGER :: xtot, ytot, ly, lx
INTEGER, DIMENSION(ng) :: group
INTEGER, DIMENSION(:), ALLOCATABLE :: bank

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>> READING CONTROL RODS INSERTION <<<<'
WRITE(ounit,*) '           --------------------------------------------'

READ(xbunit, *, IOSTAT=ios) ind, ln, nb, nstep
message = ' error in reading number of control rod bank and max. number of steps'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

READ(xbunit, *, IOSTAT=ios) ind, ln, pos0, ssize
message = ' error in reading zeroth step rod position and step size'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

ALLOCATE(bpos(nb))

!!! READ CONTROL ROD BANK POSITIONS
READ(xbunit, *, IOSTAT=ios) ind, ln, (bpos(i), i = 1, nb)
message = ' error in reading bank position'
CALL er_message(ounit, ios, ln, message, buf=xbunit)


!!! Check Control Rod Bank POSITION
DO i = 1, nb
    IF (bpos(i) > REAL(nstep)) THEN
        WRITE(ounit,1999) 'ERROR: POSITION OF CONTROL ROD BANK ', i, ' IS ', bpos(i), ' WHICH IS HIGHER THAN NUMBER OF STEPS.'
        STOP
    END IF
    IF (bpos(i) < 0.) THEN
        WRITE(ounit,1999) 'ERROR: POSITION OF CONTROL ROD BANK ', i, ' IS ', bpos(i), ' WHICH IS LOWER THAN 0.'
        STOP
    END IF
    IF (coreh < bpos(i)*ssize) THEN
        WRITE(ounit,1998) 'ERROR: CORE HEIGHT ', coreh, ' IS LOWER THAN CONTROL ROD POSITION ', bpos(i)*ssize+pos0
        WRITE(ounit,*) ' BANK NUMBER ', i
        STOP
    END IF
END DO
1999 FORMAT (2X, A, I2, A, F5.1, A)
1998 FORMAT (2X, A, F6.2, A, F6.2)

!!! READ CONTROL ROD BANK MAP
DO j = ny, 1, -1
    READ(xbunit, *, IOSTAT=ios) ind, ln, (bmap(i,j), i = 1, nx)
    message = ' error in reading control rod bank map'
    CALL er_message(ounit, ios, ln, message, buf=xbunit)
    DO i = 1, nx
        IF (bmap(i,j) > nb) THEN
            WRITE(ounit,*) '  ERROR: BANK NUMBER ON CR BANK MAP IS GREATER THAN NUMBER OF BANK'
            STOP
        END IF
    END DO
END DO

IF (bxtab == 1) THEN  !IF XTAB FILE PRESENT
  ALLOCATE(dsigtr(nnod,ng))
  ALLOCATE(dsiga (nnod,ng))
  ALLOCATE(dnuf  (nnod,ng))
  ALLOCATE(dsigf (nnod,ng))
  ALLOCATE(dsigs (nnod,ng,ng))
  ALLOCATE(ddc (nnod,6,ng))
ELSE
  ALLOCATE(dsigtr(nmat,ng))
  ALLOCATE(dsiga (nmat,ng))
  ALLOCATE(dnuf  (nmat,ng))
  ALLOCATE(dsigf (nmat,ng))
  ALLOCATE(dsigs (nmat,ng,ng))

  ! Reac CX changes due to control rod increment or dcrement
  DO i = 1, nmat
      DO g= 1, ng
          READ(xbunit, *, IOSTAT=ios) ind, ln, dsigtr(i,g), &
          dsiga(i,g), dnuf(i,g), dsigf(i,g), (dsigs(i,g,h), h = 1, ng)
          message = ' error in reading macro xs changes due to control rod insertion'
          CALL er_message(ounit, ios, ln, message, buf=xbunit)
      END DO
  END DO
END IF

!! CROD PRINT OPTION
READ(xbunit, *, IOSTAT=ios) ind, ln, popt
IF (ios == 0 .AND. popt > 0) THEN

    WRITE(ounit,1201) nb
    WRITE(ounit,1216) NINT(nstep)
    WRITE(ounit,1202) pos0
    WRITE(ounit,1203) ssize

    ALLOCATE(bank(nb))
    DO i = 1, nb
        bank(i) = i
    END DO
    WRITE(ounit,*) ' INITIAL CONTROL ROD BANK POSITION (STEPS) : '
    WRITE(ounit,*) ' (0 means fully inserted) '
    WRITE(ounit, 1204)(bank(i), i = 1, nb)
    WRITE(ounit, 1205)(bpos(i), i = 1, nb)

    WRITE(ounit,*)
    WRITE(ounit,*) ' CONTROL ROD BANK MAP : '
    DO j = ny, 1, -1
        WRITE(ounit,'(100I3)' ) (bmap(i,j), i = 1, nx)
    END DO

    IF (bxtab == 0) THEN  ! If xtab file present
      WRITE(ounit,*)
      WRITE(ounit,*) ' MATERIAL CX INCREMENT OR DECREMENT DUE TO CR INSERTION : '
      DO i= 1, nmat
         WRITE(ounit,1209) i
          WRITE(ounit,1211)'GROUP', 'TRANSPORT', 'ABSORPTION', &
          'NU*FISS', 'FISSION'
          DO g= 1, ng
              WRITE(ounit,1210) g, dsigtr(i,g), dsiga(i,g), &
              dnuf(i,g), dsigf(i,g)
              group(g) = g
          END DO
          WRITE(ounit,*)'  --SCATTERING MATRIX--'
          WRITE(ounit,'(4X, A5, 20I9)') "G/G'", (group(g), g=1,ng)
          DO g= 1, ng
              WRITE(ounit,1215)g, (dsigs(i,g,h), h=1,ng)
          END DO
      END DO
    END IF
    DEALLOCATE(bank)
END IF

1201 FORMAT(2X, 'NUMBER OF CONTROL ROD BANK  :', I3)
1216 FORMAT(2X, 'MAX. NUMBER OF STEPS        :', I4)
1202 FORMAT(2X, 'FULLY INSERTED POSITION (cm): ', F4.1, ' (FROM BOTTOM OF THE CORE)')
1203 FORMAT(2X, 'STEP SIZE (cm)              : ', F8.4)
1204 FORMAT(2X, 10(:, 2X, 'Bank ', I2))
1205 FORMAT(10(:, 2X, F7.1), /)
1209 FORMAT(4X, 'MATERIAL', I3)
1211 FORMAT(2X, A7, A12, A12, 2A13)
1210 FORMAT(2X, I6, F13.6, F12.6, 2F13.6)
1215 FORMAT(4X, I3, F14.6, 20F10.6)


!!! Convert assembly wise CR bank map to node wise CR bank map
ALLOCATE(fbmap(nxx,nyy))
ytot = 0
DO j= 1, ny
    DO ly= 1, ydiv(j)
        ytot = ytot+1
        xtot = 0
        DO i= 1, nx
            DO lx= 1, xdiv(i)
                 xtot = xtot+1
                 fbmap(xtot, ytot) = bmap(i,j)
            END DO
        END DO
    END DO
END DO


WRITE(ounit,*)
WRITE(ounit,*) ' ...Control Rods Insertion card is successfully read...'

END SUBROUTINE inp_crod

!******************************************************************************!

SUBROUTINE inp_ejct (xbunit)

!
! Purpose:
!    To read rod ejection input

USE sdata, ONLY: nf, ng, lamb, iBeta, velo, nb, tbeta, nmat, &
                 ttot, tstep1, tdiv, tstep2, bpos, fbpos, tmove, &
                 bspeed, mdir, nstep

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

INTEGER :: i, g
INTEGER :: popt
INTEGER, DIMENSION(nb) :: bank
CHARACTER(LEN=4) :: cnb         ! number of bank (character type)

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>     READING ROD EJECTION DATA      <<<<'
WRITE(ounit,*) '           --------------------------------------------'

ALLOCATE(tmove(nb), bspeed(nb), mdir(nb), fbpos(nb))
ALLOCATE(tbeta(nmat))

! Read Final CR bank position, time to start, and speed
DO i = 1, nb
    READ(xbunit, *, IOSTAT=ios) ind, ln, fbpos(i), tmove(i), bspeed(i)
    WRITE (cnb,'(I4)') nb
    cnb = TRIM(ADJUSTL(cnb))
    message = ' error in reading Final CR Bank Position, time to move and speed for bank : ' // cnb
    CALL er_message(ounit, ios, ln, message, buf=xbunit)
    IF (fbpos(i) > nstep) THEN
      WRITE(ounit, 1889) ln
      WRITE(*, 1889) ln
      STOP
      1889 FORMAT(2X, "ERROR AT LINE ", I4, ": WRONG FINAL CONTROL ROD POSITION")
    END IF
    IF (ABS(fbpos(i)-bpos(i)) < 1.e-5_DP) THEN
        mdir(i) = 0
    ELSE IF (fbpos(i)-bpos(i) > 1.e-5_DP) THEN
        mdir(i) = 2
    ELSE
        mdir(i) = 1
    END IF
END DO

! Read time for CR to be ejected
READ(xbunit, *, IOSTAT=ios) ind, ln, ttot, tstep1, tdiv, tstep2
message = ' error in time parameters'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

IF (bxtab == 0) THEN  ! IF XTAB File does not present
  ! Read beta (delayed neutron fraction)
  READ(xbunit, *, IOSTAT=ios) ind, ln, (iBeta(i), i = 1, nf)
  message = ' error in reading delayed netron fraction (beta)'
  CALL er_message(ounit, ios, ln, message, buf=xbunit)

  ! Read precusor decay constant
  READ(xbunit, *, IOSTAT=ios) ind, ln, (lamb(i), i = 1, nf)
  message = ' error in reading precusor decay constant'
  CALL er_message(ounit, ios, ln, message, buf=xbunit)

  ! Read neutron velocity
  ALLOCATE(velo(ng))
  READ(xbunit, *, IOSTAT=ios) ind, ln, (velo(g), g = 1, ng)
  message = ' error in reading neutron velocity'
  CALL er_message(ounit, ios, ln, message, buf=xbunit)
END IF


!! EJCT PRINT OPTION
READ(xbunit, *, IOSTAT=ios) ind, ln, popt
IF (ios == 0 .AND. popt > 0) THEN

    DO i = 1, nb
        bank(i) = i
    END DO
    WRITE(ounit, 1294)(bank(i), i = 1, nb)
    WRITE(ounit, 1295)(fbpos(i), i = 1, nb)
    WRITE(ounit, 1281)(tmove(i), i = 1, nb)
    WRITE(ounit, 1282)(bspeed(i), i = 1, nb)

    WRITE(ounit,*)
    WRITE(ounit,*) ' TIME PARAMETERS IN SECONDS : '
    WRITE(ounit,1297) ttot
    WRITE(ounit,1298) tstep1
    WRITE(ounit,1299) tstep2
    WRITE(ounit,1300) tdiv

    IF (bxtab == 0) THEN  ! IF XTAB File does not present
      WRITE(ounit,*)
      WRITE(ounit,*) ' DELAYED NEUTRON FRACTION : '
      WRITE(ounit,'(100F11.5)') (iBeta(i), i = 1, nf)

      WRITE(ounit,*)
      WRITE(ounit,*) ' PRECUSOR DECAY CONSTANT (1/s): '
      WRITE(ounit,'(100F11.5)') (lamb(i), i = 1, nf)

      WRITE(ounit,*)
      WRITE(ounit,*) ' NEUTRON VELOCITY (cm/s) : '
      WRITE(ounit,'(100ES15.5)') (velo(g), g = 1, ng)
    END IF
END IF

! ttot must be bigger than tstep1 and tstep2
IF ((ttot < tstep1) .OR. (ttot < tstep2)) THEN
    WRITE(ounit,*) 'ERROR: TOTAL SIMULATION TIME SHALL BE GREATER THAN TIME STEPS'
    WRITE(*,*) 'ERROR: TOTAL SIMULATION TIME SHALL BE GREATER THAN TIME STEPS'
    STOP
END IF

! tdiv must be bigger than tstep1
IF (tdiv < tstep1) THEN
    WRITE(ounit,*) 'ERROR: THE TIME WHEN SECOND TIME STEP STARTS SHALL BE GREATER THAN FIRST TIME STEP'
    WRITE(*,*) 'ERROR: THE TIME WHEN SECOND TIME STEP STARTS SHALL BE GREATER THAN FIRST TIME STEP'
    STOP
END IF

! tdiv must be less than ttot
IF (tdiv > ttot) THEN
    WRITE(ounit,*) 'ERROR: THE TIME WHEN SECOND TIME STEP STARTS SHALL BE LESS THAN TOTAL TIME'
    WRITE(*,*) 'ERROR: THE TIME WHEN SECOND TIME STEP STARTS SHALL BE LESS THAN TOTAL TIME'
    STOP
END IF

! number of steps shall be less than 10,000
IF (NINT(tdiv/tstep1)+NINT((ttot-tdiv)/tstep2) > 10000) THEN
    WRITE(ounit,*) 'ERROR: NUMBER OF TOTAL TIME STEPS ARE MORE THAN 10,000'
    WRITE(*,*) 'ERROR: NUMBER OF TOTAL TIME STEPS ARE MORE THAN 10,000'
    STOP
END IF

WRITE(ounit,*)
WRITE(ounit,*) ' ...Rod Ejection Card is successfully read...'

1294 FORMAT(25X, 99(:, 2X, 'Bank ', I2))
1295 FORMAT(2X, 'FINAL BANK POS. (STEP)', 99(:, 2X, F7.1), /)
1281 FORMAT(2X, 'STARTS MOVE (SECOND)  ', 99(:, 2X, F7.1), /)
1282 FORMAT(2X, 'SPEED (STEP/SECOND)   ', 99(:, 2X, F7.1), /)
1297 FORMAT(4X, 'TOTAL SIMULATION TIME         : ', F6.2)
1298 FORMAT(4X, 'FIRST TIME STEP               : ', F6.4)
1299 FORMAT(4X, 'SECOND TIME STEP              : ', F6.4)
1300 FORMAT(4X, 'WHEN SECOND TIME STEP APPLY?  : ', F6.2)

END SUBROUTINE inp_ejct

!******************************************************************************!

SUBROUTINE inp_cbcs (xbunit)

!
! Purpose:
!    To read boron concentration for critical boron search

USE sdata, ONLY: nmat, ng, rbcon, &
                 csigtr, csiga, cnuf, csigf, csigs


IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

INTEGER :: i, g, h
INTEGER :: popt
INTEGER, DIMENSION(ng) :: group

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>> READING BORON CONCENTRATION FOR BC SEARCH <<<<'
WRITE(ounit,*) '           --------------------------------------------------'

! Read Boron Concentration
READ(xbunit, *, IOSTAT=ios) ind, ln, rbcon
message = ' error in reading bc guess and bc reference'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

IF (bxtab == 0) THEN  ! IF XTAB File does not present
  ALLOCATE(csigtr(nmat,ng))
  ALLOCATE(csiga (nmat,ng))
  ALLOCATE(cnuf  (nmat,ng))
  ALLOCATE(csigf (nmat,ng))
  ALLOCATE(csigs (nmat,ng,ng))

  ! Read CX changes per ppm born change
  DO i = 1, nmat
      DO g= 1, ng
          READ(xbunit, *, IOSTAT=ios) ind, ln, csigtr(i,g), &
          csiga(i,g), cnuf(i,g), csigf(i,g), (csigs(i,g,h), h = 1, ng)
          message = ' error in reading macro xs changes per ppm boron changes'
          CALL er_message(ounit, ios, ln, message, buf=xbunit)
      END DO
  END DO
END IF

!! BCON PRINT OPTION
READ(xbunit, *, IOSTAT=ios) ind, ln, popt
IF (ios == 0 .AND. popt > 0) THEN

    WRITE(ounit,1422) rbcon
    IF (bxtab == 0) THEN  ! IF XTAB File does not present
      WRITE(ounit,*)
      WRITE(ounit,*) ' MATERIAL CX CHANGES PER PPM BORON CHANGES : '
      DO i= 1, nmat
         WRITE(ounit,1429) i
          WRITE(ounit,1431)'GROUP', 'TRANSPORT', 'ABSORPTION', &
          'NU*FISS', 'FISSION'
          DO g= 1, ng
              WRITE(ounit,1430) g, csigtr(i,g), csiga(i,g), &
              cnuf(i,g), csigf(i,g)
              group(g) = g
          END DO
          WRITE(ounit,*)'  --SCATTERING MATRIX--'
          WRITE(ounit,'(4X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
          DO g= 1, ng
              WRITE(ounit,1435)g, (csigs(i,g,h), h=1,ng)
          END DO
      END DO
    END IF
END IF

1422 FORMAT(2X, 'BORON CONCENTRATION REFERENCE :', F8.2)
1429 FORMAT(4X, 'MATERIAL', I3)
1431 FORMAT(2X, A9, A12, A14, A13, A14)
1430 FORMAT(2X, I6, E16.5, 3E14.5)
1435 FORMAT(4X, I3, E17.5, 20E13.5)


WRITE(ounit,*)
WRITE(ounit,*) ' ...Critical Boron Search card is successfully read...'

END SUBROUTINE inp_cbcs

!******************************************************************************!

SUBROUTINE inp_bcon (xbunit)

!
! Purpose:
!    To read boron concentration

USE sdata, ONLY: nmat, ng, bcon, rbcon, &
                 csigtr, csiga, cnuf, csigf, csigs

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

INTEGER :: i, g, h
INTEGER :: popt
INTEGER, DIMENSION(ng) :: group

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>       READING BORON CONCENTRATION        <<<<'
WRITE(ounit,*) '           --------------------------------------------------'

! Read Boron Concentration
READ(xbunit, *, IOSTAT=ios) ind, ln, bcon, rbcon
message = ' error in reading boron concentration and boron concentration reference'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

IF (bxtab == 0) THEN  ! If xtab file present
  ! Read CX changes per ppm born change
  ALLOCATE(csigtr(nmat,ng))
  ALLOCATE(csiga (nmat,ng))
  ALLOCATE(cnuf  (nmat,ng))
  ALLOCATE(csigf (nmat,ng))
  ALLOCATE(csigs (nmat,ng,ng))
  DO i = 1, nmat
      DO g= 1, ng
          READ(xbunit, *, IOSTAT=ios) ind, ln, csigtr(i,g), &
          csiga(i,g), cnuf(i,g), csigf(i,g), (csigs(i,g,h), h = 1, ng)
          message = ' error in reading macro xs changes per ppm boron changes'
          CALL er_message(ounit, ios, ln, message, buf=xbunit)
      END DO
  END DO
END IF

!! BCON PRINT OPTION
READ(xbunit, *, IOSTAT=ios) ind, ln, popt
IF (ios == 0 .AND. popt > 0) THEN

    WRITE(ounit,1221) bcon
    WRITE(ounit,1222) rbcon

    IF (bxtab == 0) THEN  ! If xtab file present
      WRITE(ounit,*)
      WRITE(ounit,*) ' MATERIAL CX CHANGES PER PPM BORON CHANGES : '
      DO i= 1, nmat
         WRITE(ounit,1229) i
          WRITE(ounit,1231)'GROUP', 'TRANSPORT', 'ABSORPTION', &
          'NU*FISS', 'FISSION'
          DO g= 1, ng
              WRITE(ounit,1230) g, csigtr(i,g), csiga(i,g), &
              cnuf(i,g), csigf(i,g)
              group(g) = g
          END DO
          WRITE(ounit,*)'  --SCATTERING MATRIX--'
          WRITE(ounit,'(4X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
          DO g= 1, ng
              WRITE(ounit,1235)g, (csigs(i,g,h), h=1,ng)
          END DO
      END DO
    END IF
END IF

1221 FORMAT(2X, 'BORON CONCENTRATION SET       :', F8.2)
1222 FORMAT(2X, 'BORON CONCENTRATION REFERENCE :', F8.2)
1229 FORMAT(4X, 'MATERIAL', I3)
1231 FORMAT(2X, A9, A12, A14, A13, A14)
1230 FORMAT(2X, I6, E16.5, 3E14.5)
1235 FORMAT(4X, I3, E17.5, 20E13.5)



WRITE(ounit,*)
WRITE(ounit,*) ' ...Boron Concentration card is successfully read...'

END SUBROUTINE inp_bcon

!******************************************************************************!

SUBROUTINE inp_ftem (xbunit)

!
! Purpose:
!    To read fuel temperature

USE sdata, ONLY: nmat, ng, nnod, ftem, rftem, &
                 fsigtr, fsiga, fnuf, fsigf, fsigs

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

REAL(DP) :: cftem
INTEGER :: i, g, h
INTEGER :: popt
INTEGER, DIMENSION(ng) :: group

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>      READING FUEL TEMPERATURE      <<<<'
WRITE(ounit,*) '           --------------------------------------------'

! Read Fuel Temperature
READ(xbunit, *, IOSTAT=ios) ind, ln, cftem, rftem
message = ' error in reading average fuel temperature and fuel temperature reference'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

! ASSIGN CFTEM to FTEM
IF (bther == 0) ALLOCATE (ftem(nnod))
ftem = cftem   !Initial guess for fuel temperature

IF (bxtab == 0) THEN  ! IF XTAB File does not present
  ! Read CX changes fuel temperature change
  ALLOCATE(fsigtr(nmat,ng), fsiga(nmat,ng), fnuf(nmat,ng), fsigf(nmat,ng), fsigs(nmat,ng,ng))
  DO i = 1, nmat
      DO g= 1, ng
          READ(xbunit, *, IOSTAT=ios) ind, ln, fsigtr(i,g), &
          fsiga(i,g), fnuf(i,g), fsigf(i,g), (fsigs(i,g,h), h = 1, ng)
          message = ' error in reading macro xs changes per fuel temperature changes'
          CALL er_message(ounit, ios, ln, message, buf=xbunit)
      END DO
  END DO
END IF

!! FTEM PRINT OPTION
READ(xbunit, *, IOSTAT=ios) ind, ln, popt
IF (ios == 0 .AND. popt > 0) THEN

    IF (bther == 0) THEN
        WRITE(ounit,1241) cftem
    ELSE
        WRITE(ounit,1256) cftem
    END IF
    WRITE(ounit,1242) rftem

    IF (bxtab == 0) THEN  ! IF XTAB File does not present
      WRITE(ounit,*)
      WRITE(ounit,*) ' MATERIAL CX CHANGES PER FUEL TEMPERATURE CHANGES : '
      DO i= 1, nmat
         WRITE(ounit,1249) i
          WRITE(ounit,1251)'GROUP', 'TRANSPORT', 'ABSORPTION', &
          'NU*FISS', 'FISSION'
          DO g= 1, ng
              WRITE(ounit,1250) g, fsigtr(i,g), fsiga(i,g), &
              fnuf(i,g), fsigf(i,g)
              group(g) = g
          END DO
          WRITE(ounit,*)'  --SCATTERING MATRIX--'
          WRITE(ounit,'(4X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
          DO g= 1, ng
              WRITE(ounit,1255)g, (fsigs(i,g,h), h=1,ng)
          END DO
      END DO
    END IF
END IF

1241 FORMAT(2X, 'AVERAGE FUEL TEMPERATURE   :', F6.2)
1242 FORMAT(2X, 'FUEL TEMPERATURE REFERENCE :', F6.2)
1249 FORMAT(4X, 'MATERIAL', I3)
1251 FORMAT(2X, A9, A12, A14, A13, A14)
1250 FORMAT(2X, I6, E16.5, 3E14.5)
1255 FORMAT(4X, I3, E17.5, 20E13.5)
1256 FORMAT(2X, 'AVERAGE FUEL TEMPERATURE   :', F6.2, '  (NOT USED)')



WRITE(ounit,*)
WRITE(ounit,*) ' ...Fuel Temperature is card successfully read...'

END SUBROUTINE inp_ftem

!******************************************************************************!

SUBROUTINE inp_mtem (xbunit)

!
! Purpose:
!    To read moderator temperature

USE sdata, ONLY: nmat, ng, nnod, mtem, rmtem, &
                 msigtr, msiga, mnuf, msigf, msigs

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

REAL(DP) :: cmtem
INTEGER :: i, g, h
INTEGER :: popt
INTEGER, DIMENSION(ng) :: group

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>   READING MODERATOR TEMPERATURE    <<<<'
WRITE(ounit,*) '           --------------------------------------------'

! Read Moderator Temperature
READ(xbunit, *, IOSTAT=ios) ind, ln, cmtem, rmtem
message = ' error in reading Moderator temperature and Moderator temperature reference'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

! ASSIGN CMTEM to MTEM
IF (bther == 0) ALLOCATE (mtem(nnod))
mtem = cmtem

IF (bxtab == 0) THEN  ! IF XTAB File does not present
  ! Read CX changes per moderator temperature change
  ALLOCATE(msigtr(nmat,ng), msiga(nmat,ng), mnuf(nmat,ng), msigf(nmat,ng), msigs(nmat,ng,ng))
  DO i = 1, nmat
      DO g= 1, ng
          READ(xbunit, *, IOSTAT=ios) ind, ln, msigtr(i,g), &
          msiga(i,g), mnuf(i,g), msigf(i,g), (msigs(i,g,h), h = 1, ng)
          message = ' error in reading macro xs changes per Moderator temperature changes'
          CALL er_message(ounit, ios, ln, message, buf=xbunit)
      END DO
  END DO
END IF

!! MTEM PRINT OPTION
READ(xbunit, *, IOSTAT=ios) ind, ln, popt
IF (ios == 0 .AND. popt > 0) THEN

    IF (bther == 0) THEN
        WRITE(ounit,1261) cmtem
    ELSE
        WRITE(ounit,1276) cmtem
    END IF
    WRITE(ounit,1262) rmtem
    IF (bxtab == 0) THEN  ! IF XTAB File does not present
      WRITE(ounit,*)
      WRITE(ounit,*) ' MATERIAL CX CHANGES PER MODERATOR TEMPERATURE CHANGES : '
      DO i= 1, nmat
         WRITE(ounit,1269) i
          WRITE(ounit,1271)'GROUP', 'TRANSPORT', 'ABSORPTION', &
          'NU*FISS', 'FISSION'
          DO g= 1, ng
              WRITE(ounit,1270) g, msigtr(i,g), msiga(i,g), &
              mnuf(i,g), msigf(i,g)
              group(g) = g
          END DO
          WRITE(ounit,*)'  --SCATTERING MATRIX--'
          WRITE(ounit,'(4X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
          DO g= 1, ng
              WRITE(ounit,1275)g, (msigs(i,g,h), h=1,ng)
          END DO
      END DO
    END IF
END IF

1261 FORMAT(2X, 'AVERAGE MODERATOR TEMPERATURE   :', F6.2)
1262 FORMAT(2X, 'MODERATOR TEMPERATURE REFERENCE :', F6.2)
1269 FORMAT(4X, 'MATERIAL', I3)
1271 FORMAT(2X, A9, A12, A14, A13, A14)
1270 FORMAT(2X, I6, E16.5, 3E14.5)
1275 FORMAT(4X, I3, E17.5, 20E13.5)
1276 FORMAT(2X, 'AVERAGE MODERATOR TEMPERATURE   :', F6.2, '  (NOT USED)')



WRITE(ounit,*)
WRITE(ounit,*) ' ...Moderator Temperature Card is successfully read...'


END SUBROUTINE inp_mtem

!******************************************************************************!

SUBROUTINE inp_cden (xbunit)

!
! Purpose:
!    To read Coolant Density

USE sdata, ONLY: nmat, ng, nnod, cden, rcden, &
                 lsigtr, lsiga, lnuf, lsigf, lsigs

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

REAL(DP) :: ccden
INTEGER :: i, g, h
INTEGER :: popt
INTEGER, DIMENSION(ng) :: group

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>       READING COOLANT DENSITY      <<<<'
WRITE(ounit,*) '           --------------------------------------------'

! Read Coolant Density
READ(xbunit, *, IOSTAT=ios) ind, ln, ccden, rcden
message = ' error in reading Coolant Density and Coolant Density reference'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

!ASSIGN CCDEN TO CDEN
IF (bther == 0) ALLOCATE (cden(nnod))
cden = ccden

IF (bxtab == 0) THEN  ! IF XTAB File does not present
  ! Read CX changes per Coolant Density change
  ALLOCATE(lsigtr(nmat,ng), lsiga(nmat,ng), lnuf(nmat,ng), lsigf(nmat,ng), lsigs(nmat,ng,ng))
  DO i = 1, nmat
      DO g= 1, ng
          READ(xbunit, *, IOSTAT=ios) ind, ln, lsigtr(i,g), &
          lsiga(i,g), lnuf(i,g), lsigf(i,g), (lsigs(i,g,h), h = 1, ng)
          message = ' error in reading macro xs changes per Coolant Density changes'
          CALL er_message(ounit, ios, ln, message, buf=xbunit)
      END DO
  END DO
END IF

!! CDEN PRINT OPTION
READ(xbunit, *, IOSTAT=ios) ind, ln, popt
IF (ios == 0 .AND. popt > 0) THEN

    IF (bther == 0) THEN
        WRITE(ounit,1361) ccden
    ELSE
        WRITE(ounit,1376) ccden
    END IF
    WRITE(ounit,1362) rcden
    IF (bxtab == 0) THEN  ! IF XTAB File does not present
      WRITE(ounit,*)
      WRITE(ounit,*) ' MATERIAL CX CHANGES PER COOLANT DENSITY CHANGES : '
      DO i= 1, nmat
         WRITE(ounit,1369) i
          WRITE(ounit,1371)'GROUP', 'TRANSPORT', 'ABSORPTION', &
          'NU*FISS', 'FISSION'
          DO g= 1, ng
              WRITE(ounit,1370) g, lsigtr(i,g), lsiga(i,g), &
              lnuf(i,g), lsigf(i,g)
              group(g) = g
          END DO
          WRITE(ounit,*)'  --SCATTERING MATRIX--'
          WRITE(ounit,'(4X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
          DO g= 1, ng
              WRITE(ounit,1375)g, (lsigs(i,g,h), h=1,ng)
          END DO
      END DO
    END IF
END IF


1361 FORMAT(2X, 'AVERAGE COOLANT DENSITY   :', F8.4)
1362 FORMAT(2X, 'COOLANT DENSITY REFERENCE :', F8.4)
1369 FORMAT(4X, 'MATERIAL', I3)
1371 FORMAT(2X, A9, A12, A14, A13, A14)
1370 FORMAT(2X, I6, E16.5, 3E14.5)
1375 FORMAT(4X, I3, E17.5, 20E13.5)
1376 FORMAT(2X, 'AVERAGE COOLANT DENSITY   :', F8.4, '  (USED AS GUESS)')



WRITE(ounit,*)
WRITE(ounit,*) ' ...Coolant Density Card is successfully read...'


END SUBROUTINE inp_cden

!******************************************************************************!

SUBROUTINE inp_ther (xbunit)

!
! Purpose:
!    To read thermalhydraulics parameters input

USE sdata, ONLY: pow, tin, nx, ny, nxx, nyy, ystag, &
                 rf, tg, tc, rg, rc, ppitch, cf, dia, cflow, dh, pi, &
                 farea, xdiv, ydiv, ystag, node_nf, nm, nt, rdel, rpos, &
                 nnod, tfm, ppow, ent, heatf, ntem, stab, ftem, mtem, cden, &
                 thunit

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

INTEGER :: i, j
INTEGER :: nfpin, ngt                              ! Number of fuel pin and guide tubes

INTEGER :: ly, lx, ytot, xtot
REAL(DP) :: cmflow
REAL(DP), DIMENSION(nx,ny) :: area
REAL(DP) :: barea, div

REAL(DP) :: dum

INTEGER :: popt

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>   READING THERMAL-HYDRAULIC DATA   <<<<'
WRITE(ounit,*) '           --------------------------------------------'

ALLOCATE (ftem(nnod), mtem(nnod), cden(nnod))

! Read Percent Power
READ(xbunit, *, IOSTAT=ios) ind, ln, ppow
message = ' error in reading percent power'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

! Read reactor Power
READ(xbunit, *, IOSTAT=ios) ind, ln, pow
message = ' error in reading reactor full thermal power'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

! Read inlet coolant temp. (Kelvin) and  FA flow rate (kg/s)
READ(xbunit, *, IOSTAT=ios) ind, ln, tin, cmflow
message = ' error in reading coolant inlet temp. and Fuel Assembly mass flow rate'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

! Read fuel pin geometry in meter
READ(xbunit, *, IOSTAT=ios) ind, ln, rf, tg, tc, ppitch
message = ' error in reading fuel meat rad., gap thickness, clad thickness and pin pitch'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

! Check gap and clad thickness
IF (tg > 0.25 * rf) THEN
    WRITE(ounit,*) '  ERROR: GAP THICKNESS IS TO LARGE (> 0.25*rf)'
    STOP
END IF
IF (tc > 0.25 * rf) THEN
    WRITE(ounit,*) '  ERROR: CLADDING THICKNESS IS TO LARGE (> 0.25*rf)'
    STOP
END IF

! Read Number of fuel pins and guide tubes
READ(xbunit, *, IOSTAT=ios) ind, ln, nfpin, ngt
message = ' error in reading number of fuel pins and guide tubes'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

! Read Fraction of heat deposited in the coolant
READ(xbunit, *, IOSTAT=ios) ind, ln, cf
message = ' error in reading fraction of heat deposited in the coolant'
CALL er_message(ounit, ios, ln, message, buf=xbunit)

if (cf < 0. .or. cf > 1.0) stop "The value of the fraction of heat " &
// "deposited in the coolant is incorrect"

! Calculate outer radius of gap and cladding
rg = rf + tg
rc = rg + tc

! Calculate fuel pin diameter
dia = 2. * rc

! Calculate hydraulic diameter
dh = dia * ((4./pi) * (ppitch/dia)**2 - 1.)

! Calculate sub-channel area
farea = ppitch**2 - 0.25*pi*dia**2

! Calculate sub-channel mass flow rate
cflow = cmflow / REAL(nfpin)

! Calculate total coolant mass flow rate and number of fuel pins per node
barea = 0.
DO j = 1, ny
    DO i = 1, nx
        area(i,j) = xsize(i)*ysize(j)             ! assembly area
        IF (area(i,j) > barea) barea = area(i,j)  ! barea => largest assembly area for ref.
    END DO
END DO

ALLOCATE(node_nf(nxx, nyy))
node_nf = 0.

ytot = 0
DO j= 1, ny
    DO ly= 1, ydiv(j)
      ytot = ytot+1
        xtot = 0
        DO i= 1, nx
            DO lx= 1, xdiv(i)
                xtot = xtot+1
                IF ((xtot >= ystag(ytot)%smin) .AND. (xtot <= ystag(ytot)%smax )) THEN
                    div = REAL (ydiv(j) * xdiv(i))            ! Number of nodes in current assembly
                    node_nf(xtot,ytot) = area(i,j) * REAL(nfpin) / (barea * div)   ! Number of fuel pin for this node
                END IF
            END DO
        END DO
    END DO
END DO

! Calculate fuel pin mesh delta and position
dum = rf / REAL(nm)

ALLOCATE(rdel(nt), rpos(nt))
DO i = 1, nm
    rdel(i) = dum
END DO

! Fuel pin mesh size
rdel(nm+1) = tg
rdel(nm+2) = tc

! Fuel pin mesh position
rpos(1) = 0.5 * rdel(1)
DO i = 2, nt
    rpos(i) = rpos(i-1) + 0.5 * (rdel(i-1) + rdel(i))
END DO

! Guess fuel and moderator temperature
ALLOCATE(tfm(nnod, nt+1)) ! Allocate fuel pin mesh temperature
tfm = 900.                !Initial guess for radial fuel pin temperature distribution
IF (bxtab == 1) THEN
  ftem = 900.; cden = 0.711; mtem = 500.
END IF

ALLOCATE(ent(nnod))
ALLOCATE(heatf(nnod))

! Initial heat-flux rate
heatf = 0.

! Open steam table file
OPEN (UNIT=thunit, FILE='st155bar', STATUS='OLD', ACTION='READ', IOSTAT = ios)
IF (ios == 0) THEN  !If steam table present in the current directory
  ! Save steam table data to stab
  DO i = 1, ntem
     READ(thunit,*) stab(i,1), stab(i,2), stab(i,3), &
                    stab(i,4), stab(i,5), stab(i,6)
  END DO
  CLOSE(UNIT=thunit)
ELSE  ! if steam table not present, use steam table data at 15.5 MPa
  stab(1,1)=543.15; stab(1,2)=0.78106745; stab(1,3)=1182595.0;
  stab(1,4)=0.820773; stab(1,5)=0.128988; stab(1,6)=0.60720
  stab(2,1)=553.15; stab(2,2)=0.76428125; stab(2,3)=1232620.0;
  stab(2,4)=0.829727; stab(2,5)=0.126380; stab(2,6)=0.59285
  stab(3,1)=563.15; stab(3,2)=0.74619716; stab(3,3)=1284170.0;
  stab(3,4)=0.846662; stab(3,5)=0.124118; stab(3,6)=0.57710
  stab(4,1)=573.15; stab(4,2)=0.72650785; stab(4,3)=1337630.0;
  stab(4,4)=0.863597; stab(4,5)=0.121856; stab(4,6)=0.55970
  stab(5,1)=583.15; stab(5,2)=0.70475081; stab(5,3)=1393570.0;
  stab(5,4)=0.915035; stab(5,5)=0.120105; stab(5,6)=0.54045
  stab(6,1)=593.15; stab(6,2)=0.68018488; stab(6,3)=1452895.0;
  stab(6,4)=0.966472; stab(6,5)=0.118354; stab(6,6)=0.51880
  stab(7,1)=603.15; stab(7,2)=0.65150307; stab(7,3)=1517175.0;
  stab(7,4)=1.166745; stab(7,5)=0.143630; stab(7,6)=0.49420
  stab(8,1)=613.15; stab(8,2)=0.61590149; stab(8,3)=1589770.0;
  stab(8,4)=1.515852; stab(8,5)=0.195931; stab(8,6)=0.46550
  stab(9,1)=617.91; stab(9,2)=0.59896404; stab(9,3)=1624307.1;
  stab(9,4)=1.681940; stab(9,5)=0.220813; stab(9,6)=0.45185
END IF

!! THER PRINT OPTION
READ(xbunit, *, IOSTAT=ios) ind, ln, popt
IF (ios == 0 .AND. popt > 0) THEN
    WRITE(ounit,1309) ppow
    WRITE(ounit,1301) pow
    WRITE(ounit,1302) tin
    WRITE(ounit,1303) cmflow
    WRITE(ounit,1304) rf
    WRITE(ounit,1305) tg
    WRITE(ounit,1306) tc
    WRITE(ounit,1310) ppitch
    WRITE(ounit,1307) cf
END IF

WRITE(ounit,*)
WRITE(ounit,*) ' ...Thermal-hydraulic Card is successfully read...'

1309 FORMAT(2X, 'REACTOR PERCENT POWER (%)                : ', F12.5)
1301 FORMAT(2X, 'REACTOR POWER (Watt)                     : ', ES12.4)
1302 FORMAT(2X, 'COOLANT INLET TEMPERATURE (Kelvin)       : ', ES12.4)
1303 FORMAT(2X, 'FUEL ASSEMBLY MASS FLOW RATE (Kg/s)      : ', ES12.4)
1304 FORMAT(2X, 'FUEL MEAT RADIUS (m)                     : ', ES12.4)
1305 FORMAT(2X, 'GAP THICKNESS (m)                        : ', ES12.4)
1306 FORMAT(2X, 'CLAD THICKNESS (m)                       : ', ES12.4)
1310 FORMAT(2X, 'PIN PITCH(m)                             : ', ES12.4)
1307 FORMAT(2X, 'FRACTION OF HEAT DEPOSITED IN COOL.      : ', ES12.4)

DEALLOCATE(xsize, ysize, zsize)

END SUBROUTINE inp_ther

!******************************************************************************!

SUBROUTINE er_message (funit, iost, ln, mess, xtab, buf)
!
! Purpose:
!    To provide error message
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: funit, iost, ln
CHARACTER(LEN=*), INTENT(IN) :: mess
INTEGER, OPTIONAL, INTENT(IN) :: xtab, buf

IF (iost < 0) THEN
  WRITE(funit,*)
  WRITE(*,*)
  WRITE(*,*)''//achar(27)//'[31m OOPS! WE FOUND AN ERROR.'//achar(27)//'[0m'

  IF (PRESENT(xtab)) THEN
    WRITE(funit, 1014) ln, xtab
  ELSE
    WRITE(funit, 1006) carr(GETLOC(uarr,buf))
    WRITE(funit, 1013) ln, farr(GETLOC(uarr,buf))
  END IF
  WRITE(funit,*) mess
  IF (PRESENT(xtab)) THEN
    WRITE(*, 1014) ln, xtab
  ELSE
    WRITE(*, 1006) carr(GETLOC(uarr,buf))
    WRITE(*, 1013) ln, farr(GETLOC(uarr,buf))
  END IF
  WRITE(*,*) mess
  1013 FORMAT(2x, 'THIS LINE NEEDS MORE INPUT DATA. LINE', I4, &
  ' IN FILE : ', A100)
  1014 FORMAT(2x, 'ERROR: LINE', I4, &
  'IN XTAB FILE FOR MATERIAL NUMBER' , I4, '. IT NEEDS MORE DATA')
  STOP
END IF

IF (iost > 0) THEN
  WRITE(funit,*)
  WRITE(*,*)
  WRITE(*,*)''//achar(27)//'[31m OOPS! WE FOUND AN ERROR.'//achar(27)//'[0m'

  IF (PRESENT(xtab)) THEN
    WRITE(funit, 1005) ln, xtab
  ELSE
    WRITE(funit, 1006) carr(GETLOC(uarr,buf))
    WRITE(funit, 1004) ln, farr(GETLOC(uarr,buf))
  END IF
  WRITE(funit,*) mess
  IF (PRESENT(xtab)) THEN
    WRITE(*, 1005) ln, xtab
  ELSE
    WRITE(*, 1006) carr(GETLOC(uarr,buf))
    WRITE(*, 1004) ln, farr(GETLOC(uarr,buf))
  END IF
  WRITE(*,*) mess
  1006 FORMAT(2X, 'ERROR: THERE IS AN ERROR IN CARD %', A4)
  1004 FORMAT(2X, 'PLEASE CHECK LINE NUMBER', I4, ' IN FILE : ', A100)
  1005 FORMAT(2X, 'ERROR: PLEASE CHECK LINE NUMBER', I4, &
  ' IN XTAB FILE FOR MATERIAL NUMBER ', I4)
  STOP
END IF

END SUBROUTINE er_message

!****************************************************************************!

FUNCTION GETLOC(arr, elm) RESULT (loc)

  !Purpose: To  get location of an element in an array

  INTEGER, DIMENSION(:), INTENT(IN)  :: arr    ! input array
  INTEGER, INTENT(IN)  :: elm                  ! element whose location is wanted to find
  INTEGER :: loc

  integer :: i

  DO i = 1, SIZE(arr)
    IF (arr(i) == elm) THEN
      loc = i
      EXIT
    END IF
  END DO

END FUNCTION GETLOC

!******************************************************************************!

SUBROUTINE AsmPow(fn)

!
! Purpose:
!    To print axially averaged assembly-wise power distribution
!

USE sdata, ONLY: nx, ny, nxx, nyy, nzz, zdel, &
                xdel, ydel, ystag, nnod, ix, iy, iz, &
                xdiv, ydiv

IMPLICIT NONE

REAL(DP), DIMENSION(:), INTENT(IN) :: fn

REAL(DP), DIMENSION(nxx, nyy, nzz) :: fx
INTEGER :: i, j, k, n
INTEGER :: ly, lx, ys, xs, yf, xf
REAL(DP) :: summ, vsumm
REAL(DP), DIMENSION(nxx, nyy) :: fnode
REAL(DP), DIMENSION(nx, ny) :: fasm
REAL(DP) :: totp
INTEGER :: nfuel
REAL(DP) :: fmax
INTEGER :: xmax, ymax
CHARACTER(LEN=6), DIMENSION(nx, ny) :: cpow

INTEGER, PARAMETER :: xm = 12
INTEGER :: ip, ipr

fx = 0._DP
DO n = 1, nnod
    fx(ix(n), iy(n), iz(n)) = fn(n)
END DO

!Calculate axially averaged node-wise distribution
fnode = 0._DP
DO j = 1, nyy
    DO i = ystag(j)%smin, ystag(j)%smax
        summ = 0._DP
        vsumm = 0._DP
        DO k = 1, nzz
            summ = summ + fx(i,j,k)*zdel(k)
            vsumm = vsumm + zdel(k)
        END DO
        fnode(i,j)= summ/vsumm
    END DO
END DO

!Calculate assembly power
nfuel = 0
totp  = 0._DP
ys = 1
yf = 0
DO j= 1, ny
    yf = yf + ydiv(j)
    xf = 0
    xs = 1
    DO i= 1, nx
        xf = xf + xdiv(i)
        summ = 0._DP
        vsumm = 0._DP
        DO ly= ys, yf
            DO lx= xs, xf
                summ = summ + fnode(lx,ly)*xdel(lx)*ydel(ly)
                vsumm = vsumm + xdel(lx)*ydel(ly)
            END DO
        END DO
        fasm(i,j) = summ / vsumm
        xs = xs + xdiv(i)
        IF (fasm(i,j) > 0._DP) nfuel = nfuel + 1
        IF (fasm(i,j) > 0._DP) totp  = totp + fasm(i,j)
    END DO
    ys = ys + ydiv(j)
END DO


! Normalize assembly power to 1._DP
xmax = 1; ymax = 1
fmax = 0._DP
DO j = 1, ny
    DO i = 1, nx
        IF (totp > 0.) fasm(i,j) = REAL(nfuel) / totp * fasm(i,j)
        IF (fasm(i,j) > fmax) THEN     ! Get max position
            xmax = i
            ymax = j
            fmax = fasm(i,j)
        END IF
        ! Convert power to character (If power == 0 convert to blank spaces)
        IF ((fasm(i,j) - 0.) < 1.e-5_DP) THEN
            cpow(i,j) = '     '
        ELSE
            WRITE (cpow(i,j),'(F6.3)') fasm(i,j)
            cpow(i,j) = TRIM(cpow(i,j))
        END IF
    END DO
END DO


! Print assembly power distribution
WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '    Radial Power Distribution'
WRITE(ounit,*) '  =============================='

ip = nx/xm
ipr = MOD(nx,xm) - 1
xs = 1; xf = xm
DO k = 1, ip
    WRITE(ounit,'(4X,100I8)') (i, i = xs, xf)
    DO j= ny, 1, -1
        WRITE(ounit,'(2X,I4,100A8)') j, (cpow(i,j), i=xs, xf)
    END DO
    WRITE(ounit,*)
    xs = xs + xm
    xf = xf + xm
END DO


WRITE(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
IF (xs+ipr > xs) THEN
    DO j= ny, 1, -1
        WRITE(ounit,'(2X,I4,100A8)') j, (cpow(i,j), i=xs, xs+ipr)
    END DO
END IF



WRITE(ounit,*)

WRITE(ounit,*) '  MAX POS.       Maximum Value'
WRITE(ounit,1101) ymax, xmax, fasm(xmax, ymax)

1101 FORMAT(2X, '(' , I3, ',', I3,')', F15.3)


END SUBROUTINE AsmPow

!******************************************************************************!

SUBROUTINE  AxiPow(fn)

!
! Purpose:
!    To print radially averaged  power distribution
!

USE sdata, ONLY: nxx, nyy, nzz, nz, zdiv, &
                vdel, ystag, nnod, ix, iy, iz, xyz, zdel, ystag, coreh

IMPLICIT NONE

REAL(DP), DIMENSION(:), INTENT(IN) :: fn

REAL(DP), DIMENSION(nxx, nyy, nzz) :: fx
INTEGER :: i, j, k, n, ztot
INTEGER :: lz
REAL(DP) :: summ, vsumm
REAL(DP), DIMENSION(nz) :: faxi
REAL(DP) :: totp
INTEGER :: nfuel
REAL(DP) :: fmax
INTEGER :: amax

fx = 0._DP
DO n = 1, nnod
    fx(ix(n), iy(n), iz(n)) = fn(n)
END DO

! Calculate Axial Power
nfuel = 0
totp  = 0._DP
ztot = 0
DO k= 1, nz
    summ = 0._DP
    vsumm = 0._DP
    DO lz= 1, zdiv(k)
        ztot = ztot + 1
        DO j = 1, nyy
            DO i = ystag(j)%smin, ystag(j)%smax
                summ = summ + fx(i,j,ztot)
                vsumm = vsumm + vdel(xyz(i,j,ztot))
            END DO
        END DO
    END DO
    faxi(k) = summ/vsumm
    IF (faxi(k) > 0._DP) nfuel = nfuel + 1
    IF (faxi(k) > 0._DP) totp  = totp + faxi(k)
END DO

! Normalize Axial power to 1._DP
fmax = 0._DP
amax = 1
DO k = 1, nz
    faxi(k) = REAL(nfuel) / totp * faxi(k)
    IF (faxi(k) > fmax) THEN
        amax = k   ! Get max position
        fmax = faxi(k)
    END IF
END DO

! Print Axial power distribution
WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '    Axial Power Density Distribution'
WRITE(ounit,*) '  ===================================='
WRITE(ounit,*)
WRITE(ounit,*) '    Plane Number        Power      Height'
WRITE(ounit,*) '   -----------------------------------------'
summ = 0.
ztot = nzz
DO k= nz, 1, -1
    IF (k == nz) THEN
        WRITE(ounit,'(2X,I8,A7,F13.3, F12.2)') k, ' (TOP)', faxi(k), coreh-summ
    ELSE IF (k == 1) THEN
        WRITE(ounit,'(2X,I8,A10,F10.3, F12.2)') k, ' (BOTTOM)', faxi(k), coreh-summ
    ELSE
        WRITE(ounit,'(2X,I8,F20.3, F12.2)') k, faxi(k), coreh-summ
    END IF
    DO lz = 1, zdiv(k)
        summ = summ + zdel(ztot)
        ztot = ztot - 1
    END DO
END DO
WRITE(ounit,*)
WRITE(ounit,*) '  MAX POS.       Maximum Value'
WRITE(ounit,1102)  amax, faxi(amax)

1102 FORMAT(4X, '(' , I3, ')', F18.3)


END SUBROUTINE AxiPow

!******************************************************************************!

SUBROUTINE  AsmFlux(fn, norm)

!
! Purpose:
!    To print axially averaged assembly-wise flux distribution
!

USE sdata, ONLY: ng, nx, ny, nxx, nyy, nzz, zdel, &
                xdel, ydel, ystag, nnod, ix, iy, iz, &
                xdiv, ydiv

IMPLICIT NONE

REAL(DP), DIMENSION(:,:), INTENT(IN) :: fn
REAL(DP), OPTIONAL, INTENT(IN) :: norm

REAL(DP), DIMENSION(nxx, nyy, nzz, ng) :: fx
INTEGER :: g, i, j, k, n
INTEGER :: ly, lx, ys, xs, yf, xf
REAL(DP) :: summ, vsumm
REAL(DP), DIMENSION(nxx, nyy, ng) :: fnode
REAL(DP), DIMENSION(nx, ny, ng) :: fasm
REAL(DP), DIMENSION(ng) :: totp
CHARACTER(LEN=10), DIMENSION(nx, ny) :: cflx

INTEGER, PARAMETER :: xm = 12
INTEGER :: ip, ipr
INTEGER :: negf

fx = 0._DP
DO g = 1, ng
    DO n = 1, nnod
        fx(ix(n), iy(n), iz(n), g) = fn(n,g)
    END DO
END DO

!Calculate axially averaged node-wise distribution
fnode = 0._DP
DO g = 1, ng
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax
            summ = 0._DP
            vsumm = 0._DP
            DO k = 1, nzz
                summ = summ + fx(i,j,k,g)*xdel(i)*ydel(j)*zdel(k)
                vsumm = vsumm + xdel(i)*ydel(j)*zdel(k)
            END DO
            fnode(i,j,g)= summ/vsumm
        END DO
    END DO
END DO

!Calculate Radial Flux (assembly wise)
negf = 0
DO g = 1, ng
    totp(g)  = 0._DP
    ys = 1
    yf = 0
    DO j= 1, ny
        yf = yf + ydiv(j)
        xf = 0
        xs = 1
        DO i= 1, nx
            xf = xf + xdiv(i)
            summ = 0._DP
            vsumm = 0._DP
            DO ly= ys, yf
                DO lx= xs, xf
                    summ = summ + fnode(lx,ly,g)*xdel(lx)*ydel(ly)
                    vsumm = vsumm + xdel(lx)*ydel(ly)
                END DO
            END DO
            fasm(i,j,g) = summ / vsumm
            xs = xs + xdiv(i)
            IF (fasm(i,j,g) > 0._DP) THEN
                totp(g)  = totp(g) + fasm(i,j,g)
            END IF
            ! Check if there is negative flux
            IF (fasm(i,j,g) < 0._DP) negf = 1
        END DO
        ys = ys + ydiv(j)
    END DO
END DO

! Normalize Flux to norm
IF (PRESENT(norm)) THEN
    DO g = 1, ng
        DO j = 1, ny
            DO i = 1, nx
                fasm(i,j,g) = norm / totp(g) * fasm(i,j,g) * norm
            END DO
        END DO
    END DO
END IF


! Print assembly Flux distribution
WRITE(ounit,*)
IF (negf > 0) WRITE(ounit,*) '    ....WARNING: NEGATIVE FLUX ENCOUNTERED....'
WRITE(ounit,*)
WRITE(ounit,*) '    Radial Flux Distribution'
WRITE(ounit,*) '  =============================='

ip = nx/xm
ipr = MOD(nx,xm) - 1

DO g = 1, ng
    WRITE(ounit,'(A,I3)') '    Group : ', g
    !!! Convert to character (zero flux convert to blank spaces)
    DO j = 1, ny
        DO i = 1, nx
            IF ((fasm(i,j,g) - 0.) < 1.e-5_DP) THEN
                cflx(i,j) = '         '
            ELSE
                WRITE (cflx(i,j),'(ES10.3)') fasm(i,j,g)
                cflx(i,j) = TRIM(ADJUSTL(cflx(i,j)))
            END IF
        END DO
    END DO

    xs = 1; xf = xm
    DO k = 1, ip
        WRITE(ounit,'(3X,100I11)') (i, i = xs, xf)
        DO j= ny, 1, -1
            WRITE(ounit,'(2X,I4,2X,100A11)') j, (cflx(i,j), i=xs, xf)
        END DO
        WRITE(ounit,*)
        xs = xs + xm
        xf = xf + xm
    END DO

    WRITE(ounit,'(3X,100I11)') (i, i = xs, xs+ipr)
    IF (xs+ipr > xs) THEN
        DO j= ny, 1, -1
            WRITE(ounit,'(2X,I4,2X,100A11)') j, (cflx(i,j), i=xs, xs+ipr)
        END DO
    END IF
   WRITE(ounit,*)

END DO


END SUBROUTINE AsmFlux

!******************************************************************************!

SUBROUTINE inp_xtab(xbunit)

!
! Purpose:
!    To read tabular xsec file
!

USE sdata, ONLY: ng, nmat, nf, m, chi

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln !Line number
INTEGER :: iost  ! IOSTAT status

INTEGER, DIMENSION(nf) :: precf
INTEGER :: i, j, g, h, s,t,u,v

! XTAB TYPE
TYPE :: XFILE
    CHARACTER(LEN=100) :: fname             ! XTAB File name
    INTEGER :: cnum                        ! Composition number in the XTAB Files
END TYPE
TYPE(XFILE), DIMENSION(:), ALLOCATABLE :: xtab

LOGICAL, DIMENSION(:), ALLOCATABLE :: noty   !to check if this buffer was read or not?
INTEGER, PARAMETER :: xunit = 998  !XTAB file unit number
INTEGER, PARAMETER :: tunit = 999  !XTAB Buffer unit number
INTEGER :: comm  ! Position of comment mark
INTEGER :: popt
INTEGER, DIMENSION(:), ALLOCATABLE :: group
INTEGER :: nskip  !Number of lines to skip

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>> READING TABULAR XSEC FILE <<<<'
WRITE(ounit,*) '           --------------------------------------------'

READ(xbunit, *, IOSTAT=iost) ind, ln, ng, nmat  !Read numbef of group and material
message = ' error in material number'
CALL er_message(ounit, iost, ln, message)

ALLOCATE(chi(nmat,ng))
ALLOCATE(xtab(nmat), noty(nmat), m(nmat))
noty = .TRUE.

! Reading input file for XTAB file names and composition number
DO i= 1, nmat
    READ(xbunit, '(A2,I5,A200)', IOSTAT=iost) ind, ln, iline
    iline = ADJUSTL(iline)                !Adjust to left
    comm = INDEX(iline, ' ')              ! Get space position
    xtab(i)%fname = iline(1:comm-1)       !Get xtab file name
    READ(iline(comm:100),*) xtab(i)%cnum  ! Get composition number (convert to integer)
    message = ' error in reading XTAB files'
    CALL er_message(ounit, iost, ln, message)
END DO

! Starting to read XTAB File, remove comments and write to buffer file
DO i = 1, nmat
  IF (noty(i)) THEN  !If this composition was not written in buffer
    ! Open XTAB File
    CALL openFile(xunit, xtab(i)%fname, 'XTAB', 'XTAB File Open Failed--status')

    ! Start removing comments and rewrite into one input XTAB buffer
    CALL inp_comments(xunit, tunit, '*')

    !This loop to read another composition in the same XTAB file
    DO j = i, nmat
      ! If next material has the same file name
      IF (TRIM(ADJUSTL(xtab(i)%fname)) == TRIM(ADJUSTL(xtab(j)%fname))) THEN

        !Read buffer file and saved the xsec and transient data
        REWIND(tunit)
        READ(tunit, *,IOSTAT=iost) ind, ln, m(j)%tadf, m(j)%trod     ! Read input control
        message = ' ERROR IN XTAB FILE ' // TRIM(ADJUSTL(xtab(j)%fname)) &
        // ': CANNOT READ CONTROL PARAMETERS'
        CALL er_message(ounit, iost, ln, message, XTAB=j)

        READ(tunit, *, IOSTAT=iost) ind, ln, m(j)%nd, m(j)%nb, m(j)%nf, m(j)%nm    ! Read number of branch
        message = ' ERROR IN XTAB FILE '// TRIM(ADJUSTL(xtab(j)%fname))&
         // ': CANNOT READ BRANCH DIMENSION'
        CALL er_message(ounit, iost, ln, message)

        ! Check branch dimension
        IF (m(j)%nd < 1 .OR. m(j)%nb < 1 .OR. m(j)%nf < 1 .OR. m(j)%nm < 1) THEN
          WRITE(ounit, *) ' ERROR: MINIMUM NUMBER OF BRANCH IS 1'
          WRITE(*, *) ' ERROR: MINIMUM NUMBER OF BRANCH IS 1'
          STOP
        END IF

        ! Allocate and read branch paramaters
        CALL branchPar(tunit, m(j)%nd, j, m(j)%pd, xtab(j)%fname, 'COOLANT DENSITY')  ! Allocate and read coolant dens. branc paramaters
        CALL branchPar(tunit, m(j)%nb, j, m(j)%pb, xtab(j)%fname, 'BORON CONCENTRATION')  ! Allocate and read Boron conc. branc paramaters
        CALL branchPar(tunit, m(j)%nf, j, m(j)%pf, xtab(j)%fname, 'FUEL TEMPERATURE')  ! Allocate and read fule temp. branc paramaters
        CALL branchPar(tunit, m(j)%nm, j, m(j)%pm, xtab(j)%fname, 'MODERATOR TEMPERATURE')  ! Allocate and read moderator temp. branc paramaters

        ! ALLOCATE XSEC DATA
        ALLOCATE(m(j)%velo(ng))
        ALLOCATE(m(j)%xsec(m(j)%nd, m(j)%nb, m(j)%nf, m(j)%nm))
        IF (m(j)%trod == 1) ALLOCATE(m(j)%rxsec(m(j)%nd, m(j)%nb, m(j)%nf, m(j)%nm))
        DO s = 1, m(j)%nd
          DO t = 1, m(j)%nb
            DO u = 1, m(j)%nf
              DO v = 1, m(j)%nm
                ALLOCATE(m(j)%xsec(s,t,u,v)%sigtr(ng))
                ALLOCATE(m(j)%xsec(s,t,u,v)%siga(ng))
                ALLOCATE(m(j)%xsec(s,t,u,v)%sigf(ng))
                ALLOCATE(m(j)%xsec(s,t,u,v)%nuf(ng))
                ALLOCATE(m(j)%xsec(s,t,u,v)%dc(ng,6))
                ALLOCATE(m(j)%xsec(s,t,u,v)%sigs(ng,ng))
                IF (m(j)%trod == 1) THEN
                  ALLOCATE(m(j)%rxsec(s,t,u,v)%sigtr(ng))
                  ALLOCATE(m(j)%rxsec(s,t,u,v)%siga(ng))
                  ALLOCATE(m(j)%rxsec(s,t,u,v)%sigf(ng))
                  ALLOCATE(m(j)%rxsec(s,t,u,v)%nuf(ng))
                  ALLOCATE(m(j)%rxsec(s,t,u,v)%dc(ng,6))
                  ALLOCATE(m(j)%rxsec(s,t,u,v)%sigs(ng,ng))
                END IF
              END DO
            END DO
          END DO
        END DO

        ! Skip lines to read desired composition in the xtab file
        nskip = ng*m(j)%nb*m(j)%nf*m(j)%nm
        IF (m(j)%tadf  == 1) THEN  ! IF dc present
          IF (m(j)%trod  == 1) THEN
            CALL skipRead(tunit, j, (xtab(j)%cnum-1)*(10*nskip+2*ng*nskip+4))
          ELSE
            CALL skipRead(tunit, j, (xtab(j)%cnum-1)*(5*nskip+ng*nskip+4))
          END IF
        ELSE IF (m(j)%tadf  == 2) THEN
          IF (m(j)%trod  == 1) THEN
            CALL skipRead(tunit, j, (xtab(j)%cnum-1)*(20*nskip+2*ng*nskip+4))
          ELSE
            CALL skipRead(tunit, j, (xtab(j)%cnum-1)*(10*nskip+ng*nskip+4))
          END IF
        ELSE
          IF (m(j)%trod  == 1) THEN
            CALL skipRead(tunit, j, (xtab(j)%cnum-1)*(8*nskip+2*ng*nskip+4))
          ELSE
            CALL skipRead(tunit, j, (xtab(j)%cnum-1)*(4*nskip+ng*nskip+4))
          END IF
        END IF

        ! Read unrodded XSEC
        CALL readXS (tunit, j, 0, m(j)%xsec)
        ! Read rodded XSEC
        IF (m(j)%trod == 1) CALL readXS (tunit, j, 1, m(j)%rxsec)
        !Read fission spectrum
        READ(tunit, *, IOSTAT=iost) ind, ln, (chi(j,g), g = 1, ng)
        message = ' ERROR IN XTAB FILE: CANNOT READ FISSION SPECTRUM'
        CALL er_message(ounit, iost, ln, message, XTAB=j)
        !Read neutron Inverse velocity
        READ(tunit, *, IOSTAT=iost) ind, ln, (m(j)%velo(g), g = 1, ng)
        message = ' ERROR IN XTAB FILE: CANNOT READ Inverse Velocity'
        CALL er_message(ounit, iost, ln, message, XTAB=j)
        DO g = 1, ng
          m(j)%velo(g) = 1._DP/m(j)%velo(g)  !COnvert to velocity
        END DO
        ! Read decay constant
        READ(tunit, *, IOSTAT=iost) ind, ln, (m(j)%lamb(t), t = 1, nf)
        message = ' ERROR IN XTAB FILE: CANNOT READ DECAY CONSTANT'
        CALL er_message(ounit, iost, ln, message, XTAB=j)
        ! Read beta
        READ(tunit, *, IOSTAT=iost) ind, ln, (m(j)%iBeta(t), t = 1, nf)
        message = ' ERROR IN XTAB FILE: CANNOT READ DELAYED NEUTRON FRACTION'
        CALL er_message(ounit, iost, ln, message, XTAB=j)

        ! If read, indicate that it has been read
        noty(j) = .FALSE.
      END IF
    END DO

    ! Close XTAB File and buffer file
    CLOSE(UNIT=xunit); CLOSE(UNIT=tunit)
  END IF
END DO

! XTAB PRINT OPTION
READ(xbunit, *, IOSTAT=iost) ind, ln, popt
IF (iost == 0 .AND. popt > 0) THEN
  s = 1; t = 1; u = 1; v = 1
  ALLOCATE(group(ng))
  DO g = 1, ng
    group(g) = g
  END DO
  DO j = 1, nf
    precf(j) = j
  END DO
  DO i= 1, nmat
      WRITE(ounit,*)
      WRITE(ounit,1709) i
      WRITE(ounit,'(A,I3)') '     XTAB FILE '// TRIM(ADJUSTL(xtab(i)%fname)) &
      // '. COMPOSITION NUMBER', xtab(i)%cnum
      WRITE(ounit,1707)'GROUP', 'TRANSPORT', 'DIFFUSION', 'ABSORPTION', &
      'NU*FISS', 'KAP*FIS','FISS. SPCTR', 'NEUTRON VELOCITY'
      DO g= 1, ng
          WRITE(ounit,1706) g, m(i)%xsec(s,t,u,v)%sigtr(g), &
          1./(3.*m(i)%xsec(s,t,u,v)%sigtr(g)), m(i)%xsec(s,t,u,v)%siga(g), &
           m(i)%xsec(s,t,u,v)%nuf(g), m(i)%xsec(s,t,u,v)%sigf(g), chi(i,g), &
           m(i)%velo(g), m(i)%xsec(s,t,u,v)%dc(g,1)
      END DO
      WRITE(ounit,*)'  --SCATTERING MATRIX--'
      WRITE(ounit,'(4X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
      DO g= 1, ng
          WRITE(ounit,1015)g, (m(i)%xsec(s,t,u,v)%sigs(g,h), h=1,ng)
      END DO
      WRITE(ounit,*)'  --BETA AND LAMBDA--'
      WRITE(ounit,'(6X, A5, I9, 20I13)') "GROUP", (precf(j), j=1,nf)
      WRITE(ounit,1016)'BETA ', (m(i)%ibeta(j), j = 1, nf)
      WRITE(ounit,1016)'LAMBDA ', (m(i)%lamb(j), j = 1, nf)
  END DO
END IF

WRITE(ounit,*)
WRITE(ounit,*) ' ...XTAB FILE Card is successfully read...'

1709 FORMAT(5X, 'MATERIAL', I3)
1707 FORMAT(2X, A7, A12, A13, A12, A11, A13, A15, A18)
1706 FORMAT(2X, I6, F13.6, 2F12.6, F13.6, ES14.5, F12.6, 2ES16.5)
1015 FORMAT(4X, I3, F16.6, 20F12.6)
1016 FORMAT(4X, A9, 20ES13.5)

DEALLOCATE(xtab, noty)

END SUBROUTINE inp_xtab

!******************************************************************************!

SUBROUTINE branchPar (tunit, dim, matnum, par, fname, messPar)

!Purpose: To allocate and read branch paramaters

IMPLICIT NONE

INTEGER, INTENT(IN) :: tunit, dim, matnum
CHARACTER(LEN=*), INTENT(IN) :: messPar, fname
REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: par

INTEGER :: k
INTEGER :: ln, iost

IF (dim > 1) THEN                         ! If branch DIMENSION > 1
  ALLOCATE(par(dim))
  READ(tunit, *, IOSTAT=iost) ind, ln, par(1:dim)
  message = ' ERROR IN XTAB FILE '// TRIM(ADJUSTL(fname))&
   // ': CANNOT READ BRANCH PARAMETERS ' // messPar
  CALL er_message(ounit, iost, ln, message, XTAB=matnum)
  DO k = 2, dim
    IF (par(k-1) > par(k)) THEN
      WRITE(ounit,*) "  ERROR IN XTAB FILE  ", fname
      WRITE(ounit,*) "  ", messPar, " PARAMETER SHALL BE IN ORDER, SMALL to BIG"
      WRITE(*,*) "  ERROR IN XTAB FILE  ", fname
      WRITE(*,*) "  ", messPar, " PARAMETER SHALL BE IN ORDER, SMALL to BIG"
      STOP
    END IF
  END DO
ELSE
  ALLOCATE(par(1))
  par(1) = 0.0    !Arbitrary
END IF

END SUBROUTINE branchPar

!******************************************************************************!

SUBROUTINE readXS (tunit, mnum, rod, xsec)
!Purpose: To read xsec in XTAB file

USE sdata, ONLY: m, ng, XBRANCH
IMPLICIT NONE

INTEGER, INTENT(IN) :: tunit, mnum, rod  ! file unit number, material number, and rod indicator
TYPE(XBRANCH), DIMENSION(:,:,:,:), INTENT(INOUT) :: xsec  !Set INOUT, see: http://www.cs.rpi.edu/~szymansk/OOF90/bugs.html#2
INTEGER :: iost, ln
INTEGER :: g, h, s, t, u, v, k

!Read sigtr
DO g = 1, ng
  DO v = 1, m(mnum)%nm
    DO u = 1, m(mnum)%nf
      DO t = 1, m(mnum)%nb
        READ(tunit, *, IOSTAT=iost) ind, ln, &
        (xsec(s,t,u,v)%sigtr(g), s = 1, m(mnum)%nd)
        IF (rod == 0) THEN
          message = ' ERROR IN XTAB FILE: CANNOT READ TRANSPORT XSEC'
        ELSE
          message = ' ERROR IN XTAB FILE: CANNOT READ RODDED TRANSPORT XSEC'
        END IF
        CALL er_message(ounit, iost, ln, message, XTAB=mnum)
      END DO
    END DO
  END DO
END DO
!Read siga
DO g = 1, ng
  DO v = 1, m(mnum)%nm
    DO u = 1, m(mnum)%nf
      DO t = 1, m(mnum)%nb
        READ(tunit, *, IOSTAT=iost) ind, ln, &
        (xsec(s,t,u,v)%siga(g), s = 1, m(mnum)%nd)
        IF (rod == 0) THEN
          message = ' ERROR IN XTAB FILE: CANNOT READ ABSORPTION XSEC'
        ELSE
          message = ' ERROR IN XTAB FILE: CANNOT READ RODDED ABSORPTION XSEC'
        END IF
        CALL er_message(ounit, iost, ln, message, XTAB=mnum)
      END DO
    END DO
  END DO
END DO
!Read nu*sigf
DO g = 1, ng
  DO v = 1, m(mnum)%nm
    DO u = 1, m(mnum)%nf
      DO t = 1, m(mnum)%nb
        READ(tunit, *, IOSTAT=iost) ind, ln, &
        (xsec(s,t,u,v)%nuf(g), s = 1, m(mnum)%nd)
        IF (rod == 0) THEN
          message = ' ERROR IN XTAB FILE: CANNOT READ NU*SIGF XSEC'
        ELSE
          message = ' ERROR IN XTAB FILE: CANNOT READ RODDED NU*SIGF XSEC'
        END IF
        CALL er_message(ounit, iost, ln, message, XTAB=mnum)
      END DO
    END DO
  END DO
END DO
!Read kappa*sigf
DO g = 1, ng
  DO v = 1, m(mnum)%nm
    DO u = 1, m(mnum)%nf
      DO t = 1, m(mnum)%nb
        READ(tunit, *, IOSTAT=iost) ind, ln, &
        (xsec(s,t,u,v)%sigf(g), s = 1, m(mnum)%nd)
        IF (rod == 0) THEN
          message = ' ERROR IN XTAB FILE: CANNOT READ KAPPA*SIGF XSEC'
        ELSE
          message = ' ERROR IN XTAB FILE: CANNOT READ RODDED KAPPA*SIGF XSEC'
        END IF
        CALL er_message(ounit, iost, ln, message, XTAB=mnum)
      END DO
    END DO
  END DO
END DO
!Read sigs
DO g = 1, ng
  DO h = 1, ng
    DO v = 1, m(mnum)%nm
      DO u = 1, m(mnum)%nf
        DO t = 1, m(mnum)%nb
          READ(tunit, *, IOSTAT=iost) ind, ln, &
          (xsec(s,t,u,v)%sigs(g,h), s = 1, m(mnum)%nd)
          IF (rod == 0) THEN
            message = ' ERROR IN XTAB FILE: CANNOT READ SCATTERING XSEC'
          ELSE
            message = ' ERROR IN XTAB FILE: CANNOT READ RODDED SCATTERING XSEC'
          END IF
          CALL er_message(ounit, iost, ln, message, XTAB=mnum)
        END DO
      END DO
    END DO
  END DO
END DO
!Read dc
IF (m(mnum)%tadf  == 1) THEN  ! IF dc present
  DO g = 1, ng
    DO v = 1, m(mnum)%nm
      DO u = 1, m(mnum)%nf
        DO t = 1, m(mnum)%nb
          READ(tunit, *, IOSTAT=iost) ind, ln, &
          (xsec(s,t,u,v)%dc(g,1), s = 1, m(mnum)%nd)
          DO k = 1, 6
            DO s = 1, m(mnum)%nd
              xsec(s,t,u,v)%dc(g,k) = xsec(s,t,u,v)%dc(g,1)
            END DO
          END DO
          IF (rod == 0) THEN
            message = ' ERROR IN XTAB FILE: CANNOT READ ADFs'
          ELSE
            message = ' ERROR IN XTAB FILE: CANNOT READ RODDED ADFs'
          END IF
          CALL er_message(ounit, iost, ln, message, XTAB=mnum)
        END DO
      END DO
    END DO
  END DO
ELSE IF (m(mnum)%tadf  == 2) THEN
  DO g = 1, ng
    DO k = 1, 6
      DO v = 1, m(mnum)%nm
        DO u = 1, m(mnum)%nf
          DO t = 1, m(mnum)%nb
            READ(tunit, *, IOSTAT=iost) ind, ln, &
            (xsec(s,t,u,v)%dc(g,k), s = 1, m(mnum)%nd)
            IF (rod == 0) THEN
              message = ' ERROR IN XTAB FILE: CANNOT READ ADFs'
            ELSE
              message = ' ERROR IN XTAB FILE: CANNOT READ RODDED TRANSPORT ADFs'
            END IF
            CALL er_message(ounit, iost, ln, message, XTAB=mnum)
          END DO
        END DO
      END DO
    END DO
  END DO
ELSE
  CONTINUE
END IF


END SUBROUTINE readXS

!******************************************************************************!

SUBROUTINE skipRead (iunit,matnum, nskip)
!Purpose: To allocate and read branch paramaters

IMPLICIT NONE

INTEGER, INTENT(IN) :: iunit, matnum, nskip
INTEGER :: i, eof

DO i = 1, nskip
  READ (iunit, *, IOSTAT=eof)
  IF (eof < 0) THEN              !Check end of file
    WRITE(ounit,1131) matnum
    WRITE(ounit,1132)
    WRITE(*,1131) matnum
    WRITE(*,1132)
    STOP
  END IF
END DO

1131 FORMAT(2X, 'ERROR: END OF FILE REACHED FOR XTAB FILE IN MATERIAL NUMBER ', I3)
1132 FORMAT(2X, 'ADPRES IS STOPPING')

END SUBROUTINE skipRead

END MODULE
