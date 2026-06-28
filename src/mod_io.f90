module io

   !=========================
   ! Input output module to read, process and echo input, as well as writing output
   ! =======================

   use sdata, only: dp

   implicit none

   save

   character(len = 100) :: iname, oname

   ! ind is sed to read x indicator in beginning of input buffer line.
   ! This to prevent reading next line
   character(len = 1) :: ind

   character(len = 200) :: message    ! error message
   character(len = 200) :: iline       ! Input line
   !Ouput options
   logical, parameter :: ogeom = .true.  ! Geometry output print option
   logical, parameter :: oxsec = .true.  ! Macroscopic CXs output print option
   logical, parameter :: scr = .true.    ! Terminal ouput print option

   ! Input, output and buffer input file unit number
   integer, parameter :: iunit = 100   !input file unit number
   integer, parameter :: ounit = 101   !output file unit number
   integer, parameter :: vunit = 103   !vtk output file unit number
   integer, parameter :: buff = 99    !input buffer file unit number (entire input)

   ! Input buffer file unit number for each card
   integer, parameter :: umode = 111, uxsec = 112, ugeom = 113, ucase = 114
   integer, parameter :: uesrc = 115, uiter = 116, uprnt = 117, uadf = 118
   integer, parameter :: ucrod = 119, ubcon = 120, uftem = 121, umtem = 122
   integer, parameter :: ucden = 123, ucbcs = 124, uejct = 125, uther = 126
   integer, parameter :: uxtab = 127, ukern = 128, uextr = 129, uthet = 130
   integer, parameter :: uoutp = 131, uvtk  = 132, uwiel = 133
   integer :: bunit

   ! Card active/inactive indicator (active = 1, inactive = 0)
   integer :: bmode = 0, bxsec = 0, bgeom = 0, bcase = 0, besrc = 0
   integer :: biter = 0, bprnt = 0, badf = 0, bcrod = 0, bbcon = 0
   integer :: bftem = 0, bmtem = 0, bcden = 0, bcbcs = 0, bejct = 0
   integer :: bther = 0, bxtab = 0, bkern = 0, bextr = 0, bthet = 0
   integer :: boutp = 0, bvtk  = 0, bwiel = 0

   ! This declaration is to notify that the error in separated card file (in case so)
   integer, parameter :: ncard = 23                  ! Number of card
   integer, dimension(ncard) :: uarr = &             ! Array of buffer unit number
   (/umode, uxsec, ugeom, ucase, uesrc, uiter, uprnt, uadf, ucrod, ubcon, &
   uftem, umtem, ucden, ucbcs, uejct, uther, uxtab, ukern, uextr, uthet, &
   uoutp, uvtk , uwiel /)
   character(len = 4), dimension(ncard) :: carr = &    ! Array of card name
   (/'MODE', 'XSEC', 'GEOM', 'CASE', 'ESRC', 'ITER', 'PRNT', 'ADF ', 'CROD', 'BCON', &
   'FTEM', 'MTEM', 'CDEN', 'CBCS', 'EJCT', 'THER', 'XTAB', 'KERN', 'EXTR', 'THET', &
   'OUTP', 'VTK ', 'WIEL' /)
   character(len = 100), dimension(:), allocatable :: farr   ! Array of card file

   ! Geometry
   integer :: np                                           ! Number of planars
   integer, dimension(:), allocatable :: zpln              ! Planar assignment to z direction
   real(dp), dimension(:), allocatable :: xsize, ysize, zsize  !Assembly size for each direction
   type :: mat_asgn                                        ! Material assignment
      integer, dimension(:, :), allocatable :: asm         ! Material assignment into assembly
      integer, dimension(:, :), allocatable :: node        ! Material assignment into nodes
   end type
   type(mat_asgn), dimension(:), allocatable :: plnr       ! planar
   integer, dimension(:, :, :), allocatable :: mnum

   ! Output option
   type :: opt_data
      logical :: assembly = .true.
   end type
   type(opt_data) :: output_opt

   contains

   !******************************************************************************!

   subroutine inp_read()
      !
      ! Purpose:
      !    [Main subroutine in this module] To read input, echo the
      !    input and gives the description
      !    to the user about his/her input
      !


      use sdata, only: ng, nnod, mode, dc, exsrc

      implicit none

      integer :: g, i, n

      !Got this trick from: http://web.utah.edu/thorne/computing/Handy_Fortran_Tricks.pdf
      n = command_argument_count()
      if (n < 1) then
         write(*, *) '  NOTE : You can also write the input directly after the command'
         write(*, '(A,A100)', advance = 'NO') '  INPUT NAME : '
         read(*, *) iname
      else
         call get_command_argument(1, iname) !Grab the first command line argument
      endif

      iname = trim(iname)

      call openfile (iunit, iname, 'input', 'Input File Open Failed--status')

      oname = trim(iname) // '.out'
      oname = trim(oname)

      open (unit = ounit, file = oname, status = 'REPLACE', action = 'WRITE')

      open (unit = umode, status = 'SCRATCH', action = 'READWRITE')
      open (unit = uxsec, status = 'SCRATCH', action = 'READWRITE')
      open (unit = ugeom, status = 'SCRATCH', action = 'READWRITE')
      open (unit = ucase, status = 'SCRATCH', action = 'READWRITE')
      open (unit = uesrc, status = 'SCRATCH', action = 'READWRITE')
      open (unit = uiter, status = 'SCRATCH', action = 'READWRITE')
      open (unit = uwiel, status = 'SCRATCH', action = 'READWRITE')
      open (unit = uprnt, status = 'SCRATCH', action = 'READWRITE')
      open (unit = uadf, status = 'SCRATCH', action = 'READWRITE')
      open (unit = ucrod, status = 'SCRATCH', action = 'READWRITE')
      open (unit = ubcon, status = 'SCRATCH', action = 'READWRITE')
      open (unit = uftem, status = 'SCRATCH', action = 'READWRITE')
      open (unit = umtem, status = 'SCRATCH', action = 'READWRITE')
      open (unit = ucden, status = 'SCRATCH', action = 'READWRITE')
      open (unit = ucbcs, status = 'SCRATCH', action = 'READWRITE')
      open (unit = uejct, status = 'SCRATCH', action = 'READWRITE')
      open (unit = uther, status = 'SCRATCH', action = 'READWRITE')
      open (unit = uxtab, status = 'SCRATCH', action = 'READWRITE')
      open (unit = ukern, status = 'SCRATCH', action = 'READWRITE')
      open (unit = uextr, status = 'SCRATCH', action = 'READWRITE')
      open (unit = uthet, status = 'SCRATCH', action = 'READWRITE')
      open (unit = uoutp, status = 'SCRATCH', action = 'READWRITE')
      open (unit = uvtk, status = 'SCRATCH', action = 'READWRITE')

      ! By default, card file names are the input file name
      allocate (farr(ncard)) ;
      farr = adjustl(iname)

      ! Echo the input to the output file
      call inp_echo()

      ! Remove comments and write the entire input to the buffer file
      call inp_comments (iunit, buff, '!')

      ! Break the buffer file and re-write into different input card buffer
      call inp_rewrite(buff)

      ! Back to the first line for all input card buffer
      rewind(umode) ;
      rewind(uxsec) ;
      rewind(ugeom) ;
      rewind(ucase) ;
      rewind(uesrc)
      rewind(uiter) ;
      rewind(uwiel) ;
      rewind(uprnt) ;
      rewind(uadf) ;
      rewind(ucrod) ;
      rewind(ubcon)
      rewind(uftem) ;
      rewind(umtem) ;
      rewind(ucden) ;
      rewind(ucbcs) ;
      rewind(uejct)
      rewind(uther) ;
      rewind(uxtab) ;
      rewind(ukern) ;
      rewind(uextr) ;
      rewind(uthet)
      rewind(uoutp) ;
      rewind(uvtk)

      ! Start reading buffer files for each input card buffer

      write(ounit, *)
      write(ounit, *)
      write(ounit, 1008)
      write(ounit, *) &
      ' ******************************************************************************'

      ! Card MODE
      if (bmode == 1) then
         call inp_mode(umode)
      else
         write(ounit, 1021) '%MODE'
         stop
      end if

      ! Card CASE
      if (bcase == 1) call inp_case (ucase)

      ! Card %KERN
      if (bkern == 1) then
         call inp_kern(ukern)
      else
         if (scr) then
            write(*, *)
            write(*, *) ' NODAL KERNEL  : SEMI-ANALYTIC NODAL METHOD'
         end if
      end if

      ! Card XSEC
      if (bxsec == 1) then
         call inp_xsec(uxsec)
      else if (bxtab == 1) then
         call inp_xtab(uxtab)
      else
         write(ounit, 1021) '%XSEC OR %XTAB'
         write(*, 1021) '%XSEC OR %XTAB'
         stop
      end if

      ! Card GEOM
      if (bgeom == 1) then
         call inp_geom1(ugeom)
         call inp_geom2(ugeom)
      else
         write(ounit, 1021) '%GEOM'
         write(*, 1021) '%GEOM'
         stop
      end if

      ! Card PRNT
      if (bprnt == 1) call inp_prnt (uprnt)

      ! Card CBCS
      if (mode == 'BCSEARCH' .and. bcbcs == 1) then
         call inp_cbcs(ucbcs)
      else if (mode == 'BCSEARCH' .and. bcbcs == 0 .and. bxtab == 0) then
         write(ounit, *) '   ERROR: CALCULATION MODE IS CRITICAL BORON CONCENTRATION SEARCH'
         write(ounit, 1041) 'CBCS', 'CRITICAL BORON CONCENTRATION SEARCH'
         write(*, *) '   ERROR: CALCULATION MODE IS CRITICAL BORON CONCENTRATION SEARCH'
         write(*, 1041) 'CBCS', 'CRITICAL BORON CONCENTRATION SEARCH'
         stop
      else if (mode /= 'BCSEARCH' .and. bcbcs == 1) then
         write(ounit, *) '   ERROR: CBCS CARD IS PROHIBITED FOR THIS CALCULATION MODE'
         write(*, *) '   ERROR: CBCS CARD IS PROHIBITED FOR THIS CALCULATION MODE'
         stop
      else if (mode == 'BCSEARCH' .and. bbcon == 1) then
         write(ounit, *) '   ERROR: BCON CARD IS PROHIBITED FOR THIS CALCULATION MODE'
         write(*, *) '   ERROR: BCON CARD IS PROHIBITED FOR THIS CALCULATION MODE'
         stop
      else
         continue
      end if

      if (mode == 'BCSEARCH' .and. bxtab == 1 .and. bther == 0) then
         write(ounit, *) '   ERROR: THER CARD REQUIRED IN THIS PROBLEM'
         write(*, *) '   ERROR: THER CARD IS REQUIRED IN THIS PROBLEM'
         stop
      end if


      !CARD CROD
      if (bcrod == 1) call inp_crod (ucrod)

      ! Card EJCT (Rod Ejection)
      if (mode == 'RODEJECT' .and. bejct == 1 .and. bcrod == 1) then
         call inp_ejct(uejct)
      else if (mode == 'RODEJECT' .and. bejct /= 1) then
         write(ounit, *) '   CALCULATION MODE ROD EJECTION'
         write(ounit, 1041) 'EJCT', 'ROD EJECTION - TRANSIENT'
         write(*, *) '   CALCULATION MODE ROD EJECTION'
         write(*, 1041) 'EJCT', 'ROD EJECTION - TRANSIENT'
         stop
      else if (mode == 'RODEJECT' .and. bcrod /= 1) then
         write(ounit, *) '   CALCULATION MODE ROD EJECTION'
         write(ounit, 1041) 'CROD', 'CONTROL ROD'
         write(*, *) '   CALCULATION MODE ROD EJECTION'
         write(*, 1041) 'CROD', 'CONTROL ROD'
         stop
      else if (mode /= 'RODEJECT' .and. bejct == 1) then
         write(ounit, *) '   EJCT CARD IS NOT NECESSARY FOR THIS CALCULATION MODE'
         write(*, *) '   EJCT CARD IS NOT NECESSARY FOR THIS CALCULATION MODE'
         stop
      else
         continue
      end if

      if (boutp == 1) call inp_outp(uoutp)

      ! Miscellaneous things
      call misc()

      ! Card ITER
      if (biter == 1) call inp_iter (uiter)

      ! Card WIEL
      if (bwiel == 1) call inp_wiel (uwiel)

      ! Card THET
      if (bthet == 1) call inp_thet (uthet)

      !!CARD THER
      if (bther == 1 .and. mode == 'FIXEDSRC') then
         write(ounit, *) '   ERROR: %THER CARD NOT VALID FOR FIXED SOURCE CALCULATION MODE'
         write(*, *) '   ERROR: %THER CARD NOT VALID FOR FIXED SOURCE CALCULATION MODE'
         stop
      else if (bther == 1 .and. mode == 'ADJOINT') then
         write(ounit, *) '   ERROR: %THER CARD NOT VALID FOR ADJOINT CALCULATION MODE'
         write(*, *) '   ERROR: %THER CARD NOT VALID FOR ADJOINT CALCULATION MODE'
         stop
      else if (bther == 1 .and. bftem == 1 .and. (bmtem == 1 .or. bcden == 1)) then
         call inp_ther (uther)
      else if (bther == 1 .and. bxtab == 1) then
         call inp_ther (uther)
      else if (bther == 0) then
         continue
      else
         if (bxtab /= 1) then
            write(ounit, *) '   ERROR: WHEN %THER CARD PRESENT %FTEM AND,' // &
            'AT LEAST ONE OF THE FOLLOWING CARDS MUST PRESENT'
            write(ounit, *) '   1. %MTEM   2. %CDEN'
            write(*, *) '   ERROR: WHEN %THER CARD PRESENT %FTEM AND,' // &
            'AT LEAST ONE OF THE FOLLOWING CARDS MUST PRESENT'
            write(*, *) '   1. %MTEM   2. %CDEN'
            stop
         end if
      end if

      if (bextr == 1) call inp_extr()

      !!CARD BCON
      if (bbcon == 1) call inp_bcon (ubcon)

      !!CARD FTEM
      if (bftem == 1) call inp_ftem (uftem)

      !!CARD MTEM
      if (bmtem == 1) call inp_mtem (umtem)

      !!CARD CDEN
      if (bcden == 1) call inp_cden (ucden)


      !CARD ADF
      allocate(dc(nnod, ng, 6))
      do g = 1, ng
         do n = 1, nnod
            dc(n, g, 1:) = 1._dp     !by default, adf = 1
         end do
      end do

      if (badf == 1 .and. bxtab == 0) then
         call inp_adf (uadf)
      else if (badf == 1 .and. bxtab == 1) then
         write(ounit, *) '  BOTH %ADF AND %XTAB CARDS CANNOT PRESENT TOGETHER'
         write(*, *) '  BOTH %ADF AND %XTAB CARDS CANNOT PRESENT TOGETHER'
         stop
      else
         continue
      end if

      ! Card ESRC
      allocate(exsrc(nnod, ng))   ! For transient or rod ejection problem, this used to
      exsrc = 0._dp               ! store transient terms
      if (mode == 'FIXEDSRC' .and. besrc == 1) then
         call inp_esrc(uesrc)
      else if (mode == 'FIXEDSRC' .and. besrc /= 1) then
         write(ounit, *) '   CALCULATION MODE IS FIXED SOURCE'
         write(ounit, 1041) 'ESRC', 'FIXED SOURCE'
         write(*, *) '   CALCULATION MODE IS FIXED SOURCE'
         write(*, 1041) 'ESRC', 'FIXED SOURCE'
         stop
      else if (mode /= 'FIXEDSRC' .and. besrc == 1) then
         write(ounit, *) '   ESRC CARD IS NOT NECESSARY FOR THIS CALCULATION MODE'
         write(*, *) '   ESRC CARD IS NOT NECESSARY FOR THIS CALCULATION MODE'
         stop
      else
         continue
      end if

      ! card VTK
      if (bvtk == 1) call inp_vtk()

      deallocate(mnum)
      do i = 1, np
         deallocate(plnr(i)%asm)
         deallocate(plnr(i)%node)
      end do
      deallocate(plnr)
      deallocate(zpln)

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) &
      ' ************************', '  STOP READING INPUT  ', '********************************'

      1008 format (30x, 'START READING INPUT')
      1021 format(2x, 'CARD ', a, ' DOES NOT PRESENT. THIS CARD IS MANDATORY')
      1041 format(2x, 'CARD ', a, ' DOES NOT PRESENT. THIS CARD IS MANDATORY FOR ', a, ' CALCULATION MODE')


      close(unit = umode) ;
      close(unit = uxsec) ;
      close(unit = ugeom) ;
      close(unit = ucase)
      close(unit = uesrc) ;
      close(unit = uiter) ;
      close(unit = uwiel) ;
      close(unit = uprnt) ;
      close(unit = uadf)
      close(unit = ucrod) ;
      close(unit = ubcon) ;
      close(unit = uftem) ;
      close(unit = umtem)
      close(unit = ucden) ;
      close(unit = ucbcs) ;
      close(unit = uejct) ;
      close(unit = uther)
      close(unit = uxtab) ;
      close(unit = ukern) ;
      close(unit = uextr) ;
      close(unit = uthet)
      close(unit = uoutp) ;
      close(unit = uvtk)
      close(unit = buff)


   end subroutine inp_read

   !******************************************************************************!

   subroutine openfile(file_unit, iname, file, er_message)

      integer :: file_unit
      character(len = *) :: iname, file, er_message
      integer :: iost

      open (unit = file_unit, file = iname, status = 'OLD', action = 'READ', &
      iostat = iost)

      if (iost /= 0) then
         write(*, 1020) er_message, iost
         write(*, *) '  CANNOT FIND ' // file // ' FILE : ', iname
         1020 format (2x, a, i6)
         stop
      end if

   end subroutine openfile

   !******************************************************************************!

   subroutine inp_echo()
      !
      ! Purpose:
      !    To rewrite the input
      !

      implicit none

      integer :: eof
      integer :: nline

      write(ounit, 2409)
      write(ounit, 2411)
      write(ounit, 2412)
      write(ounit, 2409)
      write(ounit, *)
      write(ounit, *)

      if (scr) then
         write(*, *)
         write(*, *)
         write(*, 2409)
         write(*, 2411)
         write(*, 2412)
         write(*, 2409)
#ifdef __GIT
         write(*, *)
         write(*, *) "GIT COMMIT SHA    : ", __GIT_COMMIT_HASH
         write(*, *) "GIT COMMIT DATE   : ", __GIT_DATE
         write(*, *) "GIT COMMIT BRANCH : ", __GIT_BRANCH
         write(*, *)
#endif
         write(*, *)
         write(*, *)
      end if


      write(ounit, 1002) 'STARTS'

      nline = 0
      do
         read(iunit, '(A200)', iostat = eof) iline
         nline = nline + 1
         if (eof < 0) then
            write(ounit, 1002) 'ENDS'
            write(ounit, *)
            exit
         end if
         write(ounit, 1001) nline, iline
      end do

      2409 format(11x, '###########################################################')
      2411 format(11x, '#                 KOMODO Version: 0.2                     #')
      2412 format(11x, '#           OPEN NUCLEAR REACTOR SIMULATOR                #')

      1001 format (2x, i4, ': ', a200)
      1002 format (2x, '=============================INPUT DATA', a7, &
      ' HERE===========================')
      rewind (iunit)

   end subroutine inp_echo

   !******************************************************************************!

   subroutine inp_comments (inunit, buffer, mark)
      !
      ! Purpose:
      !    To remove the comments in input and rewrite the
      !    input into input buffer. Comments marked by !.
      !

      implicit none

      integer, intent(in) :: inunit, buffer
      character(len = *), intent(in) :: mark

      integer :: ln                  ! line number
      integer :: eof, comm

      open (unit = buffer, status = 'SCRATCH', action = 'READWRITE')

      ! Start removing comments and rewrite into one input buffer
      ln = 0
      do
         ln = ln + 1
         read (inunit, '(A200)', iostat = eof) iline
         if (eof < 0) exit              !Check end of file
         iline = trim(adjustl(iline))   ! Remove trailing blanks and adjust to left
         comm = index(iline, mark)       ! Find position '!' if any
         ! If there is no '!' and no first 20 blank spaces (in case line is blank)
         if (comm == 0 .and. iline(1:20) /= '                    ') then
            write(buffer, 1012) 'x ', ln, iline
         end if
         !If the first character is not '!'
         if (comm > 1) then
            iline = iline(1:comm - 1)       ! Take only part of input
            write(buffer, 1012) 'x ', ln, iline
         end if
      end do

      rewind(buffer)

      1012 format(a2, i5, ' ', a200)

   end subroutine inp_comments

   !******************************************************************************!

   subroutine inp_rewrite (buffer)

      ! Purpose:
      ! Read previous input buffer and rewrite and break into different input buffer
      ! for particular cards. Cards identfied by %

      implicit none

      integer, intent(in) :: buffer

      integer :: ln                  ! line number
      integer :: eof, per, comm
      character(len = 100) :: fname    ! Card File name
      character(len = 4) :: card    ! Card name
      integer, parameter :: cunit = 996  !XTAB file unit number
      integer, parameter :: xunit = 997  !XTAB Buffer unit number

      do
         per = 0
         read (buffer, '(A2,I5,A200)', iostat = eof) ind, ln, iline
         if (eof < 0) exit                                                   !Check end of file

         ! If the card is placed in a separated file as indicated by keyword FILE
         if (index(iline, 'FILE') > 0 .or. index(iline, 'file') > 0) then
            backspace(buffer)
            read (buffer, '(A2,I5,A200)') ind, ln, iline
            iline = adjustl(iline)                                             ! Adjust to left
            comm = index(iline, ' ')                                           ! Get space position
            fname = trim(adjustl(iline(comm + 1:200)))                           ! Get card file name
            farr(getloc(uarr, bunit)) = fname                                 ! Change default file name for error notification
            call openfile(xunit, fname, card, 'CARD File Open Failed--status') ! Open card file
            call inp_comments(xunit, cunit, '!') ! Remove comments in card file
            ! Begin read card in a separated file
            do
               read (cunit, '(A2,I5,A200)', iostat = eof) ind, ln, iline
               if (eof < 0) exit                                                !Check end of file
               write(bunit, 1019) 'x ', ln, iline                                ! 'x' used to prevent reading next line
            end do
            close(xunit) ;
            close(cunit)
         end if

         per = index(iline, '%')
         if (per > 0) then              ! IF %card detected
            iline = iline(per + 1:200)
            iline = trim(adjustl(iline))
            card = trim(iline)
            select case (card)
               case('MODE') ;
                  bunit = umode;
                  bmode = 1
               case('XSEC') ;
                  bunit = uxsec;
                  bxsec = 1
               case('GEOM') ;
                  bunit = ugeom;
                  bgeom = 1
               case('CASE') ;
                  bunit = ucase;
                  bcase = 1
               case('ESRC') ;
                  bunit = uesrc;
                  besrc = 1
               case('ITER') ;
                  bunit = uiter;
                  biter = 1
               case('WIEL') ;
                  bunit = uwiel;
                  bwiel = 1
               case('PRNT') ;
                  bunit = uprnt;
                  bprnt = 1
               case('ADF') ;
                  bunit = uadf ;
                  badf = 1
               case('CROD') ;
                  bunit = ucrod;
                  bcrod = 1
               case('EJCT') ;
                  bunit = uejct;
                  bejct = 1
               case('CBCS') ;
                  bunit = ucbcs;
                  bcbcs = 1
               case('FTEM') ;
                  bunit = uftem;
                  bftem = 1
               case('MTEM') ;
                  bunit = umtem;
                  bmtem = 1
               case('CDEN') ;
                  bunit = ucden;
                  bcden = 1
               case('BCON') ;
                  bunit = ubcon;
                  bbcon = 1
               case('THER') ;
                  bunit = uther;
                  bther = 1
               case('XTAB') ;
                  bunit = uxtab;
                  bxtab = 1
               case('KERN') ;
                  bunit = ukern;
                  bkern = 1
               case('EXTR') ;
                  bunit = uextr;
                  bextr = 1
               case('THET') ;
                  bunit = uthet;
                  bthet = 1
               case('OUTP') ;
                  bunit = uoutp;
                  boutp = 1
               case('VTK') ;
                  bunit = uvtk ;
                  bvtk = 1
               case default

                  write(ounit, 1014) ln, iline
                  write(*, 1014) ln, iline
                  stop
            end select
         end if
         ! Write input buffer for each card
         if (per == 0) write(bunit, 1019) 'x ', ln, iline  ! 'x' used to prevent reading next line
      end do

      1014 format(2x, 'AT LINE', i3, ' : THIS IS A WRONG INPUT CARD : ', a8)
      1019 format(a2, i5, ' ', a200)

   end subroutine inp_rewrite

   !******************************************************************************!

   subroutine inp_mode (xbunit)
      !
      ! Purpose:
      !    To read case mode in input
      !

      use sdata, only: mode

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status
      character(len = 60) :: mode_desc

      read(xbunit, '(A2, I5,A100)', iostat = ios) ind, ln, mode
      message = ' error in reading MODE'
      call er_message(ounit, ios, ln, message, buf = xbunit)


      mode = trim(adjustl(mode))  ! ADJUSTL = MOVE PRECEDING BLANK TO TRAILING

      select case(mode)
         case('FORWARD')
            mode_desc = trim(adjustl('FORWARD CALCULATION'))
         case('ADJOINT')
            mode_desc = trim(adjustl('ADJOINT CALCULATION'))
         case('FIXEDSRC')
            mode_desc = trim(adjustl('FIXED SOURCE CALCULATION'))
         case('RODEJECT')
            if (bther == 0) then
               mode_desc = trim(adjustl('ROD EJECTION CALCULATION WITHOUT T-H FEEDBACK'))
            else
               mode_desc = trim(adjustl('ROD EJECTION CALCULATION WITH T-H FEEDBACK'))
            end if
         case('BCSEARCH')
            if (bther == 0) then
               mode_desc = trim(adjustl('CRITICAL BORON CONCENTRATION SEARCH' &
               // ' WITHOUT T-H FEEDBACK'))
            else
               mode_desc = trim(adjustl('CRITICAL BORON CONCENTRATION SEARCH' &
               // ' WITH T-H FEEDBACK'))
            end if
         case default
            write(ounit, 1032) mode
            write(*, 1032) mode
            stop
      end select

      write(ounit, 1031) mode_desc
      write(ounit, *)
      if (scr) then
         write(*, 1031) mode_desc
         write(*, *)
      end if

      1031 format(2x, 'CALCULATION MODE : ', a60)
      1032 format(2x, 'MODE : ', a10, ' IS UNIDENTIFIED')

   end subroutine inp_mode

   !******************************************************************************!

   subroutine inp_case (xbunit)
      !
      ! Purpose:
      !    To read case card in input
      !

      implicit none

      integer, intent(in) :: xbunit

      character(len = 100) :: case_id
      character(len = 100) :: case_exp

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status

      read(xbunit, '(A2, I5,A100)', iostat = ios) ind, ln, case_id
      message = ' error in CASE ID'
      call er_message(ounit, ios, ln, message, buf = xbunit)
      read(xbunit, '(A2, I5,A100)', iostat = ios) ind, ln, case_exp
      message = ' error in case description'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      case_id = trim(adjustl(case_id))
      case_exp = trim(adjustl(case_exp))

      write(ounit, 1006) case_id
      write(ounit, 1007) case_exp
      if (scr) then
         write(*, 1006) case_id
         write(*, 1007) case_exp
      end if
      1006 format(2x, 'CASE ID : ', a100)
      1007 format(2x, a100)

   end subroutine inp_case

   !******************************************************************************!

   subroutine inp_xsec (xbunit)
      !
      ! Purpose:
      !    To read CROSS SECTIONS card in input
      !

      use sdata, only: nmat, ng, xsigtr, xsiga, xnuf, xsigf, &
      xsigs, xd, xsigr, chi, mode, ccnuf, ccsigf

      implicit none

      integer, intent(in) :: xbunit

      integer :: i, g, h
      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status
      real(dp) :: dum
      integer, dimension(:), allocatable :: group
      integer :: comm1, comm2

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>> READING MACROSCOPIC CROSS SECTIONS <<<<'
      write(ounit, *) '           --------------------------------------------'

      read(xbunit, *, iostat = ios) ind, ln, ng, nmat
      message = ' error in material number'
      call er_message(ounit, ios, ln, message, buf = xbunit)


      allocate(group(ng))
      do g = 1, ng
         group(g) = g
      end do

      allocate(xsigtr(nmat, ng))
      allocate(xsiga (nmat, ng))
      allocate(xnuf (nmat, ng))
      allocate(xsigf (nmat, ng))
      allocate(xsigs (nmat, ng, ng))
      allocate(xd (nmat, ng))
      allocate(xsigr (nmat, ng))
      allocate(chi (nmat, ng))

      ! To ancticipate users make mistake when they use %XTAB instead of %XSEC
      read(xbunit, '(A200)') iline
      comm1 = index(iline, '/')
      comm2 = index(iline, '\')
      if (comm1 > 0 .or. comm2 > 0) then
         write(ounit, *) '  ERROR: SLASH (/) OR BACKSLASH (\) NOT ACCEPTED IN %XSEC CARD '
         write(*, *) '  ERROR: SLASH (/) OR BACKSLASH (\) NOT ACCEPTED IN %XSEC CARD '
         stop
      end if
      backspace(xbunit)

      ! Reading MACROSCOPIC CXs
      do i = 1, nmat
         do g = 1, ng
            read(xbunit, *, iostat = ios) ind, ln, xsigtr(i, g), &
            xsiga(i, g), xnuf(i, g), xsigf(i, g), &
            chi(i, g), (xsigs(i, g, h), h = 1, ng)
            message = ' error in cross section data'
            call er_message(ounit, ios, ln, message, buf = xbunit)

            ! Check CXs values
            if (xsigtr(i, g) <= 0.0) then
               write(ounit, 1020) i, g
               stop
            end if
            if (xnuf(i, g) > 0.) ccnuf = .false.
            if (xsigf(i, g) > 0.) ccsigf = .false.

            xd(i, g) = 1. / (3. * xsigtr(i, g))
            dum = 0.0
            do h = 1, ng
               if (g /= h) dum = dum + xsigs(i, g, h)
            end do
            xsigr(i, g) = xsiga(i, g) + dum
         end do
      end do


      ! Writing output
      if (oxsec) then
         do i = 1, nmat
            write(ounit, *)
            write(ounit, 1009) i
            write(ounit, 1011) 'GROUP', 'TRANSPORT', 'DIFFUSION', 'ABSORPTION', &
            'REMOVAL', 'NU*FISS', 'KAP*FIS', 'FISS. SPCTR'
            do g = 1, ng
               write(ounit, 1010) g, xsigtr(i, g), xd(i, g), xsiga(i, g), &
               xsigr(i, g), xnuf(i, g), xsigf(i, g), chi(i, g)
            end do
            write(ounit, *) '  --SCATTERING MATRIX--'
            write(ounit, '(4X, A5, 20I11)') "G/G'", (group(g), g = 1, ng)
            do g = 1, ng
               write(ounit, 1015) g, (xsigs(i, g, h), h = 1, ng)
            end do
         end do
      end if

      if (ccnuf .and. mode /= 'FIXEDSRC') then
         write(ounit, *) "ERROR: The Problem has no fission material (nu*fission for all materials are zero)"
         stop
      end if
      if (ccsigf .and. mode /= 'FIXEDSRC') then
         write(ounit, *) "ERROR: The Problem has no fission material (fission xsec for all materials are zero)"
         stop
      end if

      write(ounit, *)
      write(ounit, *) ' ...Macroscopic CX Card is successfully read...'

      1009 format(5x, 'MATERIAL', i3)
      1011 format(2x, a7, a12, a11, a12, a9, a11, a11, a14)
      1010 format(2x, i6, f13.6, f10.6, f11.6, f12.6, f11.6, es12.5, f10.6)
      1015 format(4x, i3, f16.6, 20f12.6)
      1020 format(2x, 'ERROR: Transport cross section (sigtr)is zero or negative in material: ', i3, ' ;group: ', i3)

      deallocate(group)
      deallocate(xd, xsigr)

   end subroutine inp_xsec

   !******************************************************************************!

   subroutine inp_geom1 (xbunit)
      !
      ! Purpose:
      !    To read geometry card in input (1st part)
      !

      use sdata, only: nx, ny, nz, nxx, nyy, nzz, xdel, ydel, zdel, &
      xdiv, ydiv, zdiv

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln, ios

      integer :: i, j, k, lx, ly, lz, xtot, ytot, ztot
      real(dp) :: div

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>>>READING CORE GEOMETRY<<<<<'
      write(ounit, *) '           -------------------------------'

      ! Read number of assemblies in x, y and z directions
      read(xbunit, *, iostat = ios) ind, ln, nx, ny, nz
      message = ' error in reading number assemblies'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      ! Limit values of nx, ny and nz
      ! Limit values of nx, ny and nz
      if (nx < 2) then
         write(ounit, *) '  Error: nx shall be at least 2'
         write(*, *) '  Error: nx shall be at least 2'
         stop
      end if
      if (ny < 2) then
         write(ounit, *) '  Error: ny shall be at least 2'
         write(*, *) '  Error: ny shall be at least 2'
         stop
      end if
      if (nz < 2) then
         write(ounit, *) '  Error: nz shall be at least 2'
         write(*, *) '  Error: nz shall be at least 2'
         stop
      end if

      allocate(xsize(nx), ysize(ny), zsize(nz))
      allocate(xdiv(nx), ydiv(ny), zdiv(nz))

      ! Read assembly sizes and assembly division in x, y and z directions
      ! x-direction
      read(xbunit, *, iostat = ios) ind, ln, (xsize(i), i = 1, nx)
      message = ' error in reading x-direction assembly sizes'
      call er_message(ounit, ios, ln, message, buf = xbunit)
      read(xbunit, *, iostat = ios) ind, ln, (xdiv(i), i = 1, nx)
      message = ' error in reading x-direction assembly division'
      call er_message(ounit, ios, ln, message, buf = xbunit)
      ! y-direction
      read(xbunit, *, iostat = ios) ind, ln, (ysize(j), j = 1, ny)
      message = ' error in reading y-direction assembly sizes'
      call er_message(ounit, ios, ln, message, buf = xbunit)
      read(xbunit, *, iostat = ios) ind, ln, (ydiv(j), j = 1, ny)
      message = ' error in reading y-direction assembly division'
      call er_message(ounit, ios, ln, message, buf = xbunit)
      ! z-direction
      read(xbunit, *, iostat = ios) ind, ln, (zsize(k), k = 1, nz)
      message = ' error in reading z-direction assembly sizes'
      call er_message(ounit, ios, ln, message, buf = xbunit)
      read(xbunit, *, iostat = ios) ind, ln, (zdiv(k), k = 1, nz)
      message = ' error in reading z-direction assembly division'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      !Calculate number of nodes in x, y and z directions
      !x-direction
      nxx = 0
      do i = 1, nx
         nxx = nxx + xdiv(i)
      end do
      !y-direction
      nyy = 0
      do j = 1, ny
         nyy = nyy + ydiv(j)
      end do
      !z-direction
      nzz = 0
      do k = 1, nz
         nzz = nzz + zdiv(k)
      end do

      ! Limit values of nxx, nyy and nzz
      if (nxx > 800) then
         write(ounit, *) '  Error: nxx shall be maximum 800'
         write(*, *) '  Error: nxx shall be maximum 800'
         stop
      end if
      if (nyy > 800) then
         write(ounit, *) '  Error: nyy shall be maximum 800'
         write(*, *) '  Error: nyy shall be maximum 800'
         stop
      end if
      if (nzz > 800) then
         write(ounit, *) '  Error: nzz shall be maximum 800'
         write(*, *) '  Error: nzz shall be maximum 800'
         stop
      end if

      allocate(xdel(nxx), ydel(nyy), zdel(nzz))

      !Calculate delta x, y, and z (node sizes)
      !Delta x
      xtot = 0
      do i = 1, nx
         div = xsize(i) / real(xdiv(i))
         do lx = 1, xdiv(i)
            xtot = xtot + 1
            xdel(xtot) = div
         end do
      end do
      !Delta y
      ytot = 0
      do j = 1, ny
         div = ysize(j) / real(ydiv(j))
         do ly = 1, ydiv(j)
            ytot = ytot + 1
            ydel(ytot) = div
         end do
      end do
      !Delta z
      ztot = 0
      do k = 1, nz
         div = zsize(k) / real(zdiv(k))
         do lz = 1, zdiv(k)
            ztot = ztot + 1
            zdel(ztot) = div
         end do
      end do

   end subroutine inp_geom1

   !******************************************************************************!

   subroutine inp_geom2 (xbunit)
      !
      ! Purpose:
      !    To read geometry card in input (2nd part)
      !

      use sdata, only: nx, ny, nz, nxx, nyy, nzz, xdel, ydel, zdel, &
      xwest, xeast, ynorth, ysouth, zbott, ztop, nnod, &
      xstag, ystag, xdiv, ydiv, zdiv, nmat, coreh

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln, ios

      integer :: i, j, k, lx, ly, lz, xtot, ytot, ztot
      character(len = 2), dimension(nxx, nyy) :: mmap
      integer, parameter :: xm = 40
      integer :: ip, ipr, kp
      integer :: xs, xf

      ! Reading number of planar
      read(xbunit, *, iostat = ios) ind, ln, np
      message = ' error in reading number of planars'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      !Reading planar assignment into z-direction
      allocate(zpln(nz))
      read(xbunit, *, iostat = ios) ind, ln, (zpln(k), k = 1, nz)
      message = ' error in reading planar assignment'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      do k = 1, nz
         if (zpln(k) > np) then
            write(ounit, '(2X,A15,I3,A35)') 'ERROR: PLANAR ', &
            zpln(k), ' IS GREATER THAN NUMBER OF PLANAR'
            stop
         end if
         if (zpln(k) < 1) then
            write(ounit, '(2X,A)') 'ERROR: PLANAR SHOULD BE AT LEAST 1'
            stop
         end if
      end do

      ! Reading material assignment for each planar
      allocate(plnr(np))
      do k = 1, np
         allocate(plnr(k)%asm (nx, ny))
         allocate(plnr(k)%node(nxx, nyy))
         do j = ny, 1, -1
            read(xbunit, *, iostat = ios) ind, ln, (plnr(k)%asm(i, j), i = 1, nx)
            message = ' error in reading planar'
            call er_message(ounit, ios, ln, message, buf = xbunit)
         end do
      end do

      ! Material assignment into nodes (not part for geom output)
      allocate (mnum(nxx, nyy, nzz))

      ztot = 0
      do k = 1, nz
         do lz = 1, zdiv(k)
            ztot = ztot + 1
            ytot = 0
            do j = 1, ny
               do ly = 1, ydiv(j)
                  ytot = ytot + 1
                  xtot = 0
                  do i = 1, nx
                     do lx = 1, xdiv(i)
                        xtot = xtot + 1
                        mnum(xtot, ytot, ztot) = plnr(zpln(k))%asm(i, j)
                        if (mnum(xtot, ytot, ztot) > nmat) then
                           write(ounit, '(2X,A17,I3,A37)') 'ERROR: MATERIAL ', &
                           mnum(xtot, ytot, ztot), ' IS GREATER THAN NUMBER OF MATERIAL'
                           write(*, '(2X,A17,I3,A37)') 'ERROR: MATERIAL ', &
                           mnum(xtot, ytot, ztot), ' IS GREATER THAN NUMBER OF MATERIAL'
                           stop
                        end if
                        if (mnum(xtot, ytot, ztot) < 0) then
                           write(ounit, '(2X,A)') 'ERROR: NEGATIVE MATERIAL FOUND'
                           write(*, '(2X,A)') 'ERROR: NEGATIVE MATERIAL FOUND'
                           stop
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! Assign assembly wise planar to node wise planar for output
      ztot = 0
      do k = 1, np
         ytot = 0
         do j = 1, ny
            do ly = 1, ydiv(j)
               ytot = ytot + 1
               xtot = 0
               do i = 1, nx
                  do lx = 1, xdiv(i)
                     xtot = xtot + 1
                     plnr(k)%node(xtot, ytot) = plnr(k)%asm(i, j)
                  end do
               end do
            end do
         end do
      end do

      ! -Indexing non zero material for staggered mesh-
      allocate(ystag(nyy), xstag(nxx))
      !Indexing non zero material for staggered mesh along y direction
      do j = 1, nyy
         ystag(j)%smin = nxx
         do i = 1, nxx
            if (mnum(i, j, 1) /= 0) then
               ystag(j)%smin = i
               exit
            end if
         end do
      end do

      do j = 1, nyy
         ystag(j)%smax = 0
         do i = nxx, 1, -1
            if (mnum(i, j, 1) /= 0) then
               ystag(j)%smax = i
               exit
            end if
         end do
      end do

      !Indexing non zero material for staggered mesh along x direction
      do i = 1, nxx
         xstag(i)%smin = nyy
         do j = 1, nyy
            if (mnum(i, j, 1) /= 0) then
               xstag(i)%smin = j
               exit
            end if
         end do
      end do

      do i = 1, nxx
         xstag(i)%smax = 0
         do j = nyy, 1, -1
            if (mnum(i, j, 1) /= 0) then
               xstag(i)%smax = j
               exit
            end if
         end do
      end do

      ! Checking zero material between non-zero material
      do k = 1, nzz
         do j = 1, nyy
            do i = ystag(j)%smin, ystag(j)%smax
               if (mnum(i, j, k) == 0) then
                  write(ounit, *) 'Zero material found inside core. Check material assignment'
                  stop
               end if
            end do
         end do
      end do
      do k = 1, nzz
         do i = 1, nxx
            do j = xstag(i)%smin, xstag(i)%smax
               if (mnum(i, j, k) == 0) then
                  write(ounit, *) 'Zero material found inside core. Check material assignment'
                  stop
               end if
            end do
         end do
      end do

      !Reading Boundary Conditions
      read(xbunit, *, iostat = ios) ind, ln, xeast, xwest, ynorth, ysouth, zbott, ztop
      message = ' error in reading boundary conditions'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      if (xeast > 2 .or. xwest > 2 .or. ynorth > 2 .or. ysouth > 2 &
      .or. zbott > 2 .or. ztop > 2) then
         write(*, 1019) ln
         stop
      end if


      ! Wrting core geometry output
      if (ogeom) then
         write(ounit, *) ' Number of assembly in x, y and z directions respectively :'
         write(ounit, *) nx, ny, nz
         write(ounit, *) ' Number of nodes in x, y and z directions respectively :'
         write(ounit, *) nxx, nyy, nzz
         write(ounit, *)
         write(ounit, 1016) 'x', 'x'
         write(ounit, '(2X,10F7.2)') (xdel(i), i = 1, nxx)
         write(ounit, 1016) 'y', 'y'
         write(ounit, '(2X,10F7.2)') (ydel(j), j = 1, nyy)
         write(ounit, *)

         if (nxx < 100) then
            ip = nxx / xm
            ipr = mod(nxx, xm) - 1
            do k = 1, np

               do j = 1, nyy
                  do i = 1, nxx
                     if (plnr(k)%node(i, j) == 0) then
                        mmap(i, j) = '  '
                     else
                        write (mmap(i, j), '(I2)') plnr(k)%node(i, j)
                        mmap(i, j) = trim(adjustl(mmap(i, j)))
                     end if
                  end do
               end do

               write(ounit, 1017) k
               xs = 1;
               xf = xm
               do kp = 1, ip
                  write(ounit, '(6X,100I3)') (i, i = xs, xf)
                  do j = nyy, 1, -1
                     write(ounit, '(2X,I4,1X,100A3)') j, (mmap(i, j), i = xs, xf)
                  end do
                  xs = xs + xm
                  xf = xf + xm
               end do

               write(ounit, '(6X,100I3)') (i, i = xs, xs + ipr)
               if (xs + ipr > xs) then
                  do j = nyy, 1, -1
                     write(ounit, '(2X,I4,1X,100A3)') j, (mmap(i, j), i = xs, xs + ipr)
                  end do
               end if
            end do
         end if

         write(ounit, *)
         write(ounit, 1018)
         write(ounit, *) '--------------------------------------'
         write(ounit, *) '  Plane Number     Planar Region    delta-z'
         ztot = nzz
         do k = nz, 1, -1
            do lz = 1, zdiv(k)
               if (ztot == nzz) then
                  write(ounit, '(I9, A6, I13, F15.2)') ztot, ' (TOP)', zpln(k), zdel(ztot)
               else if (ztot == 1) then
                  write(ounit, '(I9, A9, I10, F15.2)') ztot, ' (BOTTOM)', zpln(k), zdel(ztot)
               else
                  write(ounit, '(I9, I19, F15.2)') ztot, zpln(k), zdel(ztot)
               end if
               ztot = ztot - 1
            end do
         end do
         write(ounit, *)
         write(ounit, *) '  Boundary conditions'

         if (xwest == 0) then
            write(ounit, *) ' X-directed West   : ZERO FLUX'
         else if (xwest == 1) then
            write(ounit, *) ' X-directed West   : ZERO INCOMING CURRENT'
         else
            write(ounit, *) ' X-directed West   : REFLECTIVE'
         end if

         if (xeast == 0) then
            write(ounit, *) ' X-directed East   : ZERO FLUX'
         else if (xeast == 1) then
            write(ounit, *) ' X-directed East   : ZERO INCOMING CURRENT'
         else
            write(ounit, *) ' X-directed East   : REFLECTIVE'
         end if

         if (ynorth == 0) then
            write(ounit, *) ' Y-directed North  : ZERO FLUX'
         else if (ynorth == 1) then
            write(ounit, *) ' Y-directed North  : ZERO INCOMING CURRENT'
         else
            write(ounit, *) ' Y-directed North  : REFLECTIVE'
         end if

         if (ysouth == 0) then
            write(ounit, *) ' Y-directed South  : ZERO FLUX'
         else if (ysouth == 1) then
            write(ounit, *) ' Y-directed South  : ZERO INCOMING CURRENT'
         else
            write(ounit, *) ' Y-directed South  : REFLECTIVE'
         end if

         if (zbott == 0) then
            write(ounit, *) ' Z-directed Bottom : ZERO FLUX'
         else if (zbott == 1) then
            write(ounit, *) ' Z-directed Bottom : ZERO INCOMING CURRENT'
         else
            write(ounit, *) ' Z-directed Bottom : REFLECTIVE'
         end if

         if (ztop == 0) then
            write(ounit, *) ' Z-directed Top    : ZERO FLUX'
         else if (ztop == 1) then
            write(ounit, *) ' Z-directed Top    : ZERO INCOMING CURRENT'
         else
            write(ounit, *) ' Z-directed Top    : REFLECTIVE'
         end if
      end if


      1016 format(2x, a, '-directed nodes division (delta-', a, ')')
      1017 format(3x, 'Material Map for Planar Region : ', i2)
      1018 format(2x, 'Planar Region Assignment to planes.')
      1019 format(2x, 'ERROR: Wrong boundary conditions in line : ', i4)

      write(ounit, *)
      write(ounit, *) ' ...Core geometry is successfully read...'

      ! Calculate core height
      coreh = 0._dp
      do k = 1, nzz
         coreh = coreh + zdel(k)
      end do

      ! set Number of nodes (nnod)
      nnod = 0
      do k = 1, nzz
         do j = 1, nyy
            do i = ystag(j)%smin, ystag(j)%smax
               nnod = nnod + 1
            end do
         end do
      end do

   end subroutine inp_geom2

   !******************************************************************************!

   subroutine misc ()
      !
      ! Purpose:
      !    To assign material xsec to nodes
      !    To arranges the nodes into 1-D array instead of 3-D array
      !    To allocate CMFD matrix index
      !    etc


      use sdata, only: nxx, nyy, nzz, ix, iy, iz, xyz, &
      nnod, sigtr, siga, nuf, sigf, &
      sigs, d, sigr, ng, ystag, xstag, &
      xdel, ydel, zdel, vdel, nupd, &
      mat, ind, a

      implicit none

      integer :: i, j, k, n, g, noz

      allocate(ix(nnod), iy(nnod), iz(nnod))
      allocate(xyz(nxx, nyy, nzz))
      allocate(mat(nnod))

      ! Set ix, iy, iz and xyz
      n = 0
      xyz = 0
      do k = 1, nzz      ! VERY IMPORTANT: NEVER CHANGE THE ORDERS OF THESE NESTED LOOPS
         do j = 1, nyy
            do i = ystag(j)%smin, ystag(j)%smax
               n = n + 1
               ix(n) = i
               iy(n) = j
               iz(n) = k
               xyz(i, j, k) = n
               mat(n) = mnum(i, j, k)
            end do
         end do
      end do

      allocate(sigtr(nnod, ng))
      allocate(siga (nnod, ng))
      allocate(nuf (nnod, ng))
      allocate(sigf (nnod, ng))
      allocate(sigs (nnod, ng, ng))
      allocate(d (nnod, ng))
      allocate(sigr (nnod, ng))

      ! Calculate nodes' volume
      allocate(vdel(nnod))
      do i = 1, nnod
         vdel(i) = xdel(ix(i)) * ydel(iy(i)) * zdel(iz(i))
      end do

      ! Calculate number of non-zero element in the CMFD Matrix
      noz = 0
      do n = 1, nnod
         i = ix(n) ;         ! Set i, j, k
         j = iy(n) ;
         k = iz(n)
         if (k /= 1) noz = noz + 1
         if (j /= xstag(i)%smin) noz = noz + 1
         if (i /= ystag(j)%smin) noz = noz + 1
         noz = noz + 1
         if (i /= ystag(j)%smax) noz = noz + 1
         if (j /= xstag(i)%smax) noz = noz + 1
         if (k /= nzz) noz = noz + 1
      end do

      allocate(a(ng))
      do g = 1, ng
         allocate(a(g)%elmn(noz))
      end do
      allocate(ind%col(noz))
      allocate(ind%row(nnod + 1))

      ! calaculate default nodal update interval
      nupd = ceiling((nxx + nyy + nzz) / 2.5)

   end subroutine misc

   !******************************************************************************!

   subroutine inp_esrc (xbunit)
      !
      ! Purpose:
      !    To read extra sources if any
      !

      use sdata, only: exsrc, ng, nx, ny, nz, &
      xdiv, ydiv, zdiv, xyz

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln, ios  ! line number and IOSTAT status
      integer :: g, i, j, k, n
      integer :: xt, yt, zt
      integer :: it, jt, kt
      integer :: nsrc
      real(dp) :: sden                                       ! Source density
      real(dp), dimension(:), allocatable :: spec               ! Source Spectrum
      character(len = 1), dimension(:, :), allocatable :: spos ! Source position
      integer :: xpos, ypos, zpos
      real(dp) :: summ

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>>>READING EXTRA SOURCE<<<<<'
      write(ounit, *) '           -------------------------------'

      ! Read number of source
      read(xbunit, *, iostat = ios) ind, ln, nsrc
      message = ' error in reading number of extra source'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      allocate(spec(ng), spos(nx, ny))

      do n = 1, nsrc
         ! Read source density
         read(xbunit, *, iostat = ios) ind, ln, sden
         message = ' error in reading source density'
         call er_message(ounit, ios, ln, message, buf = xbunit)

         if (sden <= 0.0) then
            write(ounit, *) '  ERROR: SOURCE DENSITY SHALL BE GREATER THAN ZERO'
            stop
         end if

         ! Read source spectrum
         read(xbunit, *, iostat = ios) ind, ln, (spec(g), g = 1, ng)
         message = ' error in reading source spectrum'
         call er_message(ounit, ios, ln, message, buf = xbunit)

         ! Is total spectrum = 1._DP?
         summ = 0._dp
         do g = 1, ng
            summ = summ + spec(g)
         end do
         ! Check total spectrum
         if (abs(summ - 1._dp) > 1.e-5_dp) then
            write(ounit, *) 'TOTAL SOURCE SPECTRUM AT LINE', ln, ' IS NOT EQUAL TO 1._DP'
            stop
         end if

         ! WRITE OUTPUT
         write(ounit, '(A12,I3)') '     SOURCE ', n
         write(ounit, *) '-----------------'
         write(ounit, '(A20,ES10.3, A11)') '  Source Density  : ', sden, '  n/(m^3*s)'
         write(ounit, '(A19,100F6.2)') '  Source Spectrum : ', (spec(g), g = 1, ng)
         write(ounit, *) ' Source Position '

         ! Read source position
         do
            read(xbunit, *, iostat = ios) ind, ln, zpos
            message = ' error in reading axial position (zpos) of extra source'
            call er_message(ounit, ios, ln, message, buf = xbunit)
            if (zpos < 1) exit
            if (zpos > nz) then
               write(ounit, *) '  ERROR: WRONG EXTRA SOURCES POSITION (ZPOS)'
               write(ounit, 2033) ln, zpos
               stop
            end if
            spos = '0'
            do
               read(xbunit, *, iostat = ios) ind, ln, xpos, ypos
               message = ' error in reading radial position (xpos and ypos) of extra source'
               call er_message(ounit, ios, ln, message, buf = xbunit)
               if (xpos < 1 .or. ypos < 1) exit

               if (xpos > nx) then
                  write(ounit, *) '  ERROR: WRONG EXTRA SOURCES POSITION (XPOS)'
                  write(ounit, 2033) ln, xpos, ypos
                  stop
               end if

               if (ypos > ny) then
                  write(ounit, *) '  ERROR: WRONG EXTRA SOURCES POSITION (YPOS)'
                  write(ounit, 2033) ln, xpos, ypos
                  stop
               end if

               spos(xpos, ypos) = 'X'

               zt = 0
               kt = 1
               do k = 1, zpos
                  if (k > 1) kt = zt + 1
                  zt = zt + zdiv(k)
               end do

               yt = 0
               jt = 1
               do j = 1, ypos
                  if (j > 1) jt = yt + 1
                  yt = yt + ydiv(j)
               end do

               xt = 0
               it = 1
               do i = 1, xpos
                  if (i > 1) it = xt + 1
                  xt = xt + xdiv(i)
               end do

               do k = kt, zt
                  do j = jt, yt
                     do i = it, xt
                        do g = 1, ng
                           exsrc(xyz(i, j, k), g) = exsrc(xyz(i, j, k), g) + &
                           sden * spec(g)
                        end do
                     end do
                  end do
               end do

            end do

            write(ounit, '(A18,I3)') '   Plane number : ', zpos
            write(ounit, '(7X,100I3)') (i, i = 1, nx)
            do j = ny, 1, -1
               write(ounit, '(4X,I3, 100A3 )') j, (spos(i, j), i = 1, nx)
            end do
            write(ounit, *)
         end do
      end do

      2033 format(2x, 'LINE', i4, ' : ', i3, i3)

      deallocate(spec, spos)

   end subroutine inp_esrc

   !******************************************************************************!

   subroutine inp_iter (xbunit)

      !
      ! Purpose:
      !    To read iteration control if any

      use sdata, only: nout, nin, serc, ferc, nac, nupd, th_niter, nth, kern, &
      matrix_solver, inner_atol, inner_rtol

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status

      read(xbunit, *, iostat = ios) ind, ln, nout, nin, serc, ferc, nac, nupd, &
      th_niter, nth, matrix_solver, inner_atol, inner_rtol
      message = ' error in reading iteration control'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>>>READING ITERATION CONTROL<<<<<'
      write(ounit, *) '           ------------------------------------'

      write(ounit, '(A,I5)') '  MAXIMUM NUMBER OF OUTER ITERATION                     : ', nout
      write(ounit, '(A,I5)') '  MAXIMUM NUMBER OF INNER ITERATION                     : ', nin
      write(ounit, '(A,ES12.3)') '  FISSION SOURCE ERROR CRITERIA                     : ', serc
      write(ounit, '(A,ES12.3)') '  FLUX ERROR CRITERIA                               : ', ferc
      write(ounit, '(A,I5)') '  OUTER ITERATION FISSION SOURCE EXTRAPOLATION INTERVAL : ', nac
      write(ounit, '(A,I5)') '  NODAL UPDATE INTERVAL                                 : ', nupd
      write(ounit, '(A,I5)') '  MAX. NUMBER OF T-H ITERATION                          : ', th_niter
      write(ounit, '(A,I5)') '  MAX. NUMBER OF OUTER ITERATION PER T-H ITERATION      : ', nth
      write(ounit, '(A,A)')  '  LINEAR SYSTEM SOLVER METHOD                           : ', trim(adjustl(matrix_solver))
      write(ounit, '(A,ES12.3)')  '  LINEAR SYSTEM INNER ABSOLUTE TOLERANCE           : ', inner_atol
      write(ounit, '(A,ES12.3)')  '  LINEAR SYSTEM INNER RELATIVE TOLERANCE           : ', inner_rtol

      if (nout < nupd .and. kern /= ' FDM') then
         write(*, *) "ERROR: MAX. NUMBER OF OUTER ITERATION SHOULD BE BIGGER THAN NODAL UPDATE INTERVAL"
         write(ounit, *) "ERROR: MAX. NUMBER OF OUTER ITERATION SHOULD BE BIGGER THAN NODAL UPDATE INTERVAL"
         stop
      end if

      if (nth < nupd .and. kern /= ' FDM' .and. bther == 1) then
         write(*, *) "ERROR: MAX. NUMBER OF OUTER ITERATION PER T-H ITERATION SHOULD BE BIGGER THAN NODAL UPDATE INTERVAL"
         write(ounit, *) "ERROR: MAX. NUMBER OF OUTER ITERATION PER T-H ITERATION SHOULD BE BIGGER THAN NODAL UPDATE INTERVAL"
         stop
      end if

      select case (matrix_solver)
         case ('prec_cg')
         case ('cg')
         case ('bicg')
         case default
            write(*, *) "ERROR: INVALID LINEAR SYSTEM SOLVER METHOD: ", trim(adjustl(matrix_solver))
            write(ounit, *) "ERROR: INVALID LINEAR SYSTEM SOLVER METHOD: ", trim(adjustl(matrix_solver))
            stop
      end select

   end subroutine inp_iter

   !******************************************************************************!

   subroutine inp_wiel (xbunit)

      !
      ! Purpose:
      !    To read iteration control if any

      use sdata, only: use_wielandt_shift, wielandt_init_iter, wielandt_shift

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status

      read(xbunit, *, iostat = ios) ind, ln, use_wielandt_shift, wielandt_init_iter, wielandt_shift
      message = ' error in reading wielandt shift control'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>>>READING WIELANDT   CONTROL<<<<<'
      write(ounit, *) '           ------------------------------------'

      write(ounit, '(A,L1)') '  USE WIELANDT SHIFT                                    : ', use_wielandt_shift
      write(ounit, '(A,I5)') '  WIELANDT SHIFT INITIAL ITERATIONS                     : ', wielandt_init_iter
      write(ounit, '(A,ES12.3)') '  WIELANDT SHIFT VALUE                              : ', wielandt_shift
      if (wielandt_shift < 0.0_dp) then
         write(ounit, '(A,A)')  '  WIELANDT SHIFT OPTION                                 : ', 'RELATIVE'
      else
         write(ounit, '(A,A)')  '  WIELANDT SHIFT OPTION                                 : ', 'ABSOLUTE'
      endif

   end subroutine inp_wiel

   !******************************************************************************!

   subroutine inp_extr ()

      !
      ! Purpose:
      !    To tell user the flux transformation is active

      implicit none

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>>>READING FLUX TRANSFORMATION<<<<<'
      write(ounit, *) '           -------------------------------------'
      write(ounit, 1571)

      1571 format (2x, 'EXPONENTIAL FLUX TRANSFORMATION IS ACTIVE')

   end subroutine inp_extr

   !******************************************************************************!

   subroutine inp_kern (xbunit)

      !
      ! Purpose:
      !    To determine nodal kernel

      use sdata, only: kern

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status
      character (len = 30) :: kern_desc


      read(xbunit, *, iostat = ios) ind, ln, kern
      message = ' error in reading iteration control'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      kern = trim(adjustr(kern))
      if (kern /= ' FDM' .and. kern /= ' PNM' .and. kern /= 'SANM') then
         write(*, 1646) kern
         write(*, *) ' NODAL KERNEL OPTIONS: FDM, PNM OR SANM'
         write(ounit, 1646) kern
         write(ounit, *) ' NODAL KERNEL OPTIONS: FDM, PNM OR SANM'
         stop
      end if

      if (kern == ' FDM') then
         kern_desc = "FINITE DIFFERENCE METHOD"
      else if (kern == ' PNM') then
         kern_desc = "POLYNOMIAL NODAL METHOD"
      else
         kern_desc = "SEMI-ANALYTIC NODAL METHOD"
      end if
      kern_desc = trim(adjustl(kern_desc))

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '             >>>>>READING NODAL KERNEL<<<<<'
      write(ounit, *) '           ------------------------------------'

      write(ounit, 1648) kern_desc
      if (scr) then
         write(*, *)
         write(*, 1648) kern_desc
      end if

      1646 format (2x, 'ERROR: COULD NOT RECOGNIZE NODAL KERNEL: ', a4)
      1648 format (2x, 'NODAL KERNEL  : ', a30)

   end subroutine inp_kern

   !******************************************************************************!

   subroutine inp_outp (xbunit)

      !
      ! Purpose:
      !    To read output option

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status
      character (len = 100) :: output_option


      read(xbunit, *, iostat = ios) ind, ln, output_option
      if (ios /= 0) output_opt % assembly = .true.

      if (trim(adjustl(output_option)) == "node-wise") then
         output_opt % assembly = .false.
      else
         output_opt % assembly = .true.
      end if

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '             >>>>>DETAILED OUTPUT OPTION<<<<<'
      write(ounit, *) '           ------------------------------------'

      write(ounit, 1650)
      if (scr) then
         write(*, *)
         write(*, 1650) 'THREE-DIMENSIONAL NODE POWER'
      end if

      1650 format (2x, 'OUTPUT OPTION   : ', a30)

   end subroutine inp_outp

   !******************************************************************************!

   subroutine inp_thet (xbunit)

      !
      ! Purpose:
      !    To read iteration theta value for transient problem

      use sdata, only: sth, bth

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status

      read(xbunit, *, iostat = ios) ind, ln, sth
      message = ' error in theta in %THETA card'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      if (sth < 0.001) then
         write(*, *)
         write(*, *) "ERROR: THETA VALUE IS TOO SMALL"
         stop
      end if
      if (sth > 1.0) then
         write(*, *)
         write(*, *) "ERROR: THETA VALUE SHALL NOT GREATER THAT 1.0"
         stop
      end if
      bth = (1._dp - sth) / sth

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '               >>>>>READING THETA CARD<<<<<'
      write(ounit, *) '           ------------------------------------'

      write(ounit, '(A,F6.2)') '  THETA IS  : ', sth


   end subroutine inp_thet

   !******************************************************************************!

   subroutine inp_prnt (xbunit)

      !
      ! Purpose:
      !    To read output print option if any

      use sdata, only: aprad, apaxi, afrad

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status

      character(len = 3) :: caprad = 'YES', capaxi = 'YES', cafrad = 'YES'

      read(xbunit, *, iostat = ios) ind, ln, aprad, apaxi, afrad
      message = ' error in reading output print option'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      if (aprad == 0) caprad = 'NO'
      if (apaxi == 0) capaxi = 'NO'
      if (afrad == 0) cafrad = 'NO'

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>>>READING OUTPUT PRINT OPTION<<<<<'
      write(ounit, *) '           -------------------------------------'

      write(ounit, '(A,A)') '  RADIAL ASSEMBLY POWER DISTRIBUTION : ', caprad
      write(ounit, '(A,A)') '  AXIAL ASSEMBLY POWER DISTRIBUTION  : ', capaxi
      write(ounit, '(A,A)') '  RADIAL FLUX POWER DISTRIBUTION     : ', cafrad


   end subroutine inp_prnt

   !******************************************************************************!

   subroutine inp_adf (xbunit)

      !
      ! Purpose:
      !    To read ADF values if any

      use sdata, only: ng, nmat, nx, ny, nz, nxx, nyy, nzz, &
      xdiv, ydiv, zdiv, xyz, ystag, dc

      implicit none

      integer, intent(in) :: xbunit

      type :: adf_type
         real(dp), dimension(6) :: dc
      end type
      type(adf_type), dimension(nmat, ng) :: mdc
      type(adf_type), dimension(nx, ny, nz, ng) :: xdc
      type(adf_type), dimension(nxx, nyy, nzz, ng) :: xxdc

      integer :: g, i, j, k, u
      integer :: rot, x1, x2, y1, y2, z1, z2
      integer :: xtot, ytot, ztot
      integer :: lz, ly, lx
      integer, dimension(nx + 1) :: tx
      integer, dimension(ny + 1) :: ty
      integer, dimension(nz + 1) :: tz
      integer :: zp
      character(len = 6), dimension(nx, ny) :: cadf
      integer, parameter :: xm = 20
      integer :: ip, ipr
      integer :: xs, xf

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '         >>>>>READING ASSEMBLY DISCONTINUITY FACTOR<<<<<'
      write(ounit, *) '         -----------------------------------------------'

      do i = 1, nmat
         do g = 1, ng
            read(xbunit, *, iostat = ios) ind, ln, (mdc(i, g)%dc(j), j = 1, 6)
            message = ' error in reading Assembly Discontinuity Factor (ADF)'
            call er_message(ounit, ios, ln, message, buf = xbunit)
         end do
      end do

      do g = 1, ng
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  if (plnr(zpln(k))%asm(i, j) /= 0) xdc(i, j, k, g) = mdc(plnr(zpln(k))%asm(i, j), g)
               end do
            end do
         end do
      end do

      !!! ADF ROTATION
      do
         read(xbunit, *, iostat = ios) ind, ln, rot
         message = ' error in reading ADF Rotation'
         call er_message(ounit, ios, ln, message, buf = xbunit)
         if (rot < 1) exit
         if (rot > 3) then
            write(ounit, *) '  ERROR: MAXIMUM ADF ROTATION IS 3 TIMES'
            write(ounit, 2030) ln, rot
         end if

         do
            read(xbunit, *, iostat = ios) ind, ln, x1, x2, y1, y2, z1, z2
            message = ' error in reading ADF Rotation'
            call er_message(ounit, ios, ln, message, buf = xbunit)
            if (x1 < 1 .or. x2 < 1 .or. y1 < 1 .or. y2 < 1 .or. z1 < 1 .or. z2 < 1) exit
            if (x1 > nx .or. x2 > nx .or. y1 > ny .or. y2 > ny .or. z1 > nz .or. z2 > nz) then
               write(ounit, *) '  ERROR: WRONG POSITION FOR ADF ROTATION. ' // &
               'OUT OF DIMENSION OF THE CORE'
               write(ounit, 2032) ln, x1, x2, y1, y2, z1, z2
               stop
            end if
            if (x2 < x1 .or. y2 < y1 .or. z2 < z1) then
               write(ounit, *) '  ERROR: WRONG POSITION FOR ADF ROTATION. ' // &
               'INITIAL POSITION IS SMALLER'
               write(ounit, 2032) ln, x1, x2, y1, y2, z1, z2
            end if

            do g = 1, ng
               do k = z1, z2
                  do j = y1, y2
                     do i = x1, x2
                        call rotate(rot, xdc(i, j, k, g)%dc(1), xdc(i, j, k, g)%dc(2), &
                        xdc(i, j, k, g)%dc(3), xdc(i, j, k, g)%dc(4))
                     end do
                  end do
               end do
            end do

         end do
      end do

      ! ADF PRINT OPTION
      read(xbunit, *, iostat = ios) ind, ln, zp
      if (ios == 0 .and. zp >= 1) then
         write(ounit, *)
         write(ounit, '(A,I3)') '  ADF VALUES ON PLANAR NUMBER : ', zp

         ip = nx / xm
         ipr = mod(nx, xm) - 1
         do g = 1, ng
            write(ounit, *)
            write(ounit, 1999) g
            write(ounit, *)

            write(ounit, *) '  EAST ADF'
            do j = 1, ny
               do i = 1, nx
                  !!! If ADF > 0, Convert ADF to character
                  if ((xdc(i, j, zp, g)%dc(1) - 0.) < 1.e-5_dp) then
                     cadf(i, j) = '      '
                  else
                     write (cadf(i, j), '(F6.4)') xdc(i, j, zp, g)%dc(1)
                     cadf(i, j) = trim(adjustl(cadf(i, j)))
                  end if
               end do
            end do

            xs = 1;
            xf = xm
            do k = 1, ip
               write(ounit, '(4X,100I8)') (i, i = xs, xf)
               do j = ny, 1, -1
                  write(ounit, '(2X,I4,100A8)') j, (cadf(i, j), i = xs, xf)
               end do
               write(ounit, *)
               xs = xs + xm
               xf = xf + xm
            end do

            write(ounit, '(4X,100I8)') (i, i = xs, xs + ipr)
            if (xs + ipr > xs) then
               do j = ny, 1, -1
                  write(ounit, '(2X,I4,100A8)') j, (cadf(i, j), i = xs, xs + ipr)
               end do
            end if


            write(ounit, *) '  WEST ADF'
            do j = 1, ny
               do i = 1, nx
                  !!! If ADF > 0, Convert ADF to character
                  if ((xdc(i, j, zp, g)%dc(2) - 0.) < 1.e-5_dp) then
                     cadf(i, j) = '      '
                  else
                     write (cadf(i, j), '(F6.4)') xdc(i, j, zp, g)%dc(2)
                     cadf(i, j) = trim(adjustl(cadf(i, j)))
                  end if
               end do
            end do

            xs = 1;
            xf = xm
            do k = 1, ip
               write(ounit, '(4X,100I8)') (i, i = xs, xf)
               do j = ny, 1, -1
                  write(ounit, '(2X,I4,100A8)') j, (cadf(i, j), i = xs, xf)
               end do
               write(ounit, *)
               xs = xs + xm
               xf = xf + xm
            end do

            write(ounit, '(4X,100I8)') (i, i = xs, xs + ipr)
            if (xs + ipr > xs) then
               do j = ny, 1, -1
                  write(ounit, '(2X,I4,100A8)') j, (cadf(i, j), i = xs, xs + ipr)
               end do
            end if


            write(ounit, *) '  NORTH ADF'
            do j = 1, ny
               do i = 1, nx
                  !!! If ADF > 0, Convert ADF to character
                  if ((xdc(i, j, zp, g)%dc(3) - 0.) < 1.e-5_dp) then
                     cadf(i, j) = '      '
                  else
                     write (cadf(i, j), '(F6.4)') xdc(i, j, zp, g)%dc(3)
                     cadf(i, j) = trim(adjustl(cadf(i, j)))
                  end if
               end do
            end do

            xs = 1;
            xf = xm
            do k = 1, ip
               write(ounit, '(4X,100I8)') (i, i = xs, xf)
               do j = ny, 1, -1
                  write(ounit, '(2X,I4,100A8)') j, (cadf(i, j), i = xs, xf)
               end do
               write(ounit, *)
               xs = xs + xm
               xf = xf + xm
            end do

            write(ounit, '(4X,100I8)') (i, i = xs, xs + ipr)
            if (xs + ipr > xs) then
               do j = ny, 1, -1
                  write(ounit, '(2X,I4,100A8)') j, (cadf(i, j), i = xs, xs + ipr)
               end do
            end if


            write(ounit, *) '  SOUTH ADF'
            do j = 1, ny
               do i = 1, nx
                  !!! If ADF > 0, Convert ADF to character
                  if ((xdc(i, j, zp, g)%dc(4) - 0.) < 1.e-5_dp) then
                     cadf(i, j) = '      '
                  else
                     write (cadf(i, j), '(F6.4)') xdc(i, j, zp, g)%dc(4)
                     cadf(i, j) = trim(adjustl(cadf(i, j)))
                  end if
               end do
            end do

            xs = 1;
            xf = xm
            do k = 1, ip
               write(ounit, '(4X,100I8)') (i, i = xs, xf)
               do j = ny, 1, -1
                  write(ounit, '(2X,I4,100A8)') j, (cadf(i, j), i = xs, xf)
               end do
               write(ounit, *)
               xs = xs + xm
               xf = xf + xm
            end do

            write(ounit, '(4X,100I8)') (i, i = xs, xs + ipr)
            if (xs + ipr > xs) then
               do j = ny, 1, -1
                  write(ounit, '(2X,I4,100A8)') j, (cadf(i, j), i = xs, xs + ipr)
               end do
            end if


            write(ounit, *) '  TOP ADF'
            do j = 1, ny
               do i = 1, nx
                  !!! If ADF > 0, Convert ADF to character
                  if ((xdc(i, j, zp, g)%dc(5) - 0.) < 1.e-5_dp) then
                     cadf(i, j) = '      '
                  else
                     write (cadf(i, j), '(F6.4)') xdc(i, j, zp, g)%dc(5)
                     cadf(i, j) = trim(adjustl(cadf(i, j)))
                  end if
               end do
            end do

            xs = 1;
            xf = xm
            do k = 1, ip
               write(ounit, '(4X,100I8)') (i, i = xs, xf)
               do j = ny, 1, -1
                  write(ounit, '(2X,I4,100A8)') j, (cadf(i, j), i = xs, xf)
               end do
               write(ounit, *)
               xs = xs + xm
               xf = xf + xm
            end do

            write(ounit, '(4X,100I8)') (i, i = xs, xs + ipr)
            if (xs + ipr > xs) then
               do j = ny, 1, -1
                  write(ounit, '(2X,I4,100A8)') j, (cadf(i, j), i = xs, xs + ipr)
               end do
            end if


            write(ounit, *) '  BOTTOM ADF'
            do j = 1, ny
               do i = 1, nx
                  !!! If ADF > 0, Convert ADF to character
                  if ((xdc(i, j, zp, g)%dc(6) - 0.) < 1.e-5_dp) then
                     cadf(i, j) = '      '
                  else
                     write (cadf(i, j), '(F6.4)') xdc(i, j, zp, g)%dc(6)
                     cadf(i, j) = trim(adjustl(cadf(i, j)))
                  end if
               end do
            end do

            xs = 1;
            xf = xm
            do k = 1, ip
               write(ounit, '(4X,100I8)') (i, i = xs, xf)
               do j = ny, 1, -1
                  write(ounit, '(2X,I4,100A8)') j, (cadf(i, j), i = xs, xf)
               end do
               write(ounit, *)
               xs = xs + xm
               xf = xf + xm
            end do

            write(ounit, '(4X,100I8)') (i, i = xs, xs + ipr)
            if (xs + ipr > xs) then
               do j = ny, 1, -1
                  write(ounit, '(2X,I4,100A8)') j, (cadf(i, j), i = xs, xs + ipr)
               end do
            end if
            write(ounit, *)
         end do
      end if

      1999 format (4x, 'GROUP : ', i3)

      write(ounit, *) '  ...Assembly Discontinuity Factors are successfully read...'


      tx(1) = 1
      do i = 2, nx + 1
         tx(i) = tx(i - 1) + xdiv(i - 1)
      end do
      ty(1) = 1
      do j = 2, ny + 1
         ty(j) = ty(j - 1) + ydiv(j - 1)
      end do
      tz(1) = 1
      do k = 2, nz + 1
         tz(k) = tz(k - 1) + zdiv(k - 1)
      end do

      do g = 1, ng
         ztot = 0
         do k = 1, nz
            do lz = 1, zdiv(k)
               ztot = ztot + 1
               ytot = 0
               do j = 1, ny
                  do ly = 1, ydiv(j)
                     ytot = ytot + 1
                     xtot = 0
                     do i = 1, nx
                        do lx = 1, xdiv(i)
                           xtot = xtot + 1
                           xxdc(xtot, ytot, ztot, g)%dc = 0._dp
                           if (mnum(xtot, ytot, ztot) /= 0) xxdc(xtot, ytot, ztot, g)%dc = 1._dp
                           if (xtot == tx(i)) xxdc(xtot, ytot, ztot, g)%dc(2) = xdc(i, j, k, g)%dc(2)
                           if (xtot == tx(i + 1) - 1) xxdc(xtot, ytot, ztot, g)%dc(1) = xdc(i, j, k, g)%dc(1)
                           if (ytot == ty(j)) xxdc(xtot, ytot, ztot, g)%dc(4) = xdc(i, j, k, g)%dc(4)
                           if (ytot == ty(j + 1) - 1) xxdc(xtot, ytot, ztot, g)%dc(3) = xdc(i, j, k, g)%dc(3)
                           if (ztot == tz(k)) xxdc(xtot, ytot, ztot, g)%dc(5) = xdc(i, j, k, g)%dc(5)
                           if (ztot == tz(k + 1) - 1) xxdc(xtot, ytot, ztot, g)%dc(6) = xdc(i, j, k, g)%dc(6)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      do g = 1, ng
         do k = 1, nzz
            do j = 1, nyy
               do i = ystag(j)%smin, ystag(j)%smax
                  do u = 1, 6
                     dc(xyz(i, j, k), g, u) = xxdc(i, j, k, g)%dc(u)
                  end do
               end do
            end do
         end do
      end do


      2030 format(3x, 'LINE', i4, ' : ', i3)
      2032 format(3x, 'LINE', i4, ' : ', i3, i3, i3, i3, i3, i3)


   end subroutine inp_adf

   !******************************************************************************!

   subroutine rotate(rot, a1, a2, a3, a4)

      ! Purpose:
      !           To rotate ADF values (necessary for BWR assemblies)

      integer, intent(in) :: rot
      real(dp), intent(inout) :: a1, a2, a3, a4
      real(dp) :: x1, x2, x3, x4


      x1 = a1
      x2 = a2
      x3 = a3
      x4 = a4

      if (rot == 1) then
         a1 = x4
         a2 = x3
         a3 = x1
         a4 = x2
      end if

      if (rot == 2) then
         a1 = x2
         a2 = x1
         a3 = x4
         a4 = x3
      end if

      if (rot == 3) then
         a1 = x3
         a2 = x4
         a3 = x2
         a4 = x1
      end if

   end subroutine rotate

   !******************************************************************************!

   subroutine inp_crod (xbunit)

      !
      ! Purpose:
      !    To read control rod position

      use sdata, only: nx, ny, nmat, ng, xdiv, ydiv, &
      nxx, nyy, bpos, nb, nstep, pos0, ssize, &
      dsigtr, dsiga, dnuf, dsigf, dsigs, ddc, nnod, coreh, fbmap

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios

      integer :: i, j, g, h
      integer, dimension(nx, ny) :: bmap       ! Radial control rod bank map (assembly wise)
      integer :: popt
      integer :: xtot, ytot, ly, lx
      integer, dimension(ng) :: group
      integer, dimension(:), allocatable :: bank

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>> READING CONTROL RODS INSERTION <<<<'
      write(ounit, *) '           --------------------------------------------'

      read(xbunit, *, iostat = ios) ind, ln, nb, nstep
      message = ' error in reading number of control rod bank and max. number of steps'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      read(xbunit, *, iostat = ios) ind, ln, pos0, ssize
      message = ' error in reading zeroth step rod position and step size'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      allocate(bpos(nb))

      !!! READ CONTROL ROD BANK POSITIONS
      read(xbunit, *, iostat = ios) ind, ln, (bpos(i), i = 1, nb)
      message = ' error in reading bank position'
      call er_message(ounit, ios, ln, message, buf = xbunit)


      !!! Check Control Rod Bank POSITION
      do i = 1, nb
         if (bpos(i) > real(nstep)) then
            write(ounit, 1999) 'ERROR: POSITION OF CONTROL ROD BANK ', i, ' IS ', bpos(i), ' WHICH IS HIGHER THAN NUMBER OF STEPS.'
            stop
         end if
         if (bpos(i) < 0.) then
            write(ounit, 1999) 'ERROR: POSITION OF CONTROL ROD BANK ', i, ' IS ', bpos(i), ' WHICH IS LOWER THAN 0.'
            stop
         end if
         if (coreh < bpos(i) * ssize) then
            write(ounit, 1998) 'ERROR: CORE HEIGHT ', coreh, ' IS LOWER THAN CONTROL ROD POSITION ', bpos(i) * ssize + pos0
            write(ounit, *) ' BANK NUMBER ', i
            stop
         end if
      end do
      1999 format (2x, a, i2, a, f5.1, a)
      1998 format (2x, a, f6.2, a, f6.2)

      !!! READ CONTROL ROD BANK MAP
      do j = ny, 1, -1
         read(xbunit, *, iostat = ios) ind, ln, (bmap(i, j), i = 1, nx)
         message = ' error in reading control rod bank map'
         call er_message(ounit, ios, ln, message, buf = xbunit)
         do i = 1, nx
            if (bmap(i, j) > nb) then
               write(ounit, *) '  ERROR: BANK NUMBER ON CR BANK MAP IS GREATER THAN NUMBER OF BANK'
               stop
            end if
         end do
      end do

      if (bxtab == 1) then  !IF XTAB FILE PRESENT
         allocate(dsigtr(nnod, ng))
         allocate(dsiga (nnod, ng))
         allocate(dnuf (nnod, ng))
         allocate(dsigf (nnod, ng))
         allocate(dsigs (nnod, ng, ng))
         allocate(ddc (nnod, 6, ng))
      else
         allocate(dsigtr(nmat, ng))
         allocate(dsiga (nmat, ng))
         allocate(dnuf (nmat, ng))
         allocate(dsigf (nmat, ng))
         allocate(dsigs (nmat, ng, ng))

         ! Reac CX changes due to control rod increment or dcrement
         do i = 1, nmat
            do g = 1, ng
               read(xbunit, *, iostat = ios) ind, ln, dsigtr(i, g), &
               dsiga(i, g), dnuf(i, g), dsigf(i, g), (dsigs(i, g, h), h = 1, ng)
               message = ' error in reading macro xs changes due to control rod insertion'
               call er_message(ounit, ios, ln, message, buf = xbunit)
            end do
         end do
      end if

      !! CROD PRINT OPTION
      read(xbunit, *, iostat = ios) ind, ln, popt
      if (ios == 0 .and. popt > 0) then

         write(ounit, 1201) nb
         write(ounit, 1216) nint(nstep)
         write(ounit, 1202) pos0
         write(ounit, 1203) ssize

         allocate(bank(nb))
         do i = 1, nb
            bank(i) = i
         end do
         write(ounit, *) ' INITIAL CONTROL ROD BANK POSITION (STEPS) : '
         write(ounit, *) ' (0 means fully inserted) '
         write(ounit, 1204) (bank(i), i = 1, nb)
         write(ounit, 1205) (bpos(i), i = 1, nb)

         write(ounit, *)
         write(ounit, *) ' CONTROL ROD BANK MAP : '
         do j = ny, 1, -1
            write(ounit, '(100I3)') (bmap(i, j), i = 1, nx)
         end do

         if (bxtab == 0) then  ! If xtab file present
            write(ounit, *)
            write(ounit, *) ' MATERIAL CX INCREMENT OR DECREMENT DUE TO CR INSERTION : '
            do i = 1, nmat
               write(ounit, 1209) i
               write(ounit, 1211) 'GROUP', 'TRANSPORT', 'ABSORPTION', &
               'NU*FISS', 'FISSION'
               do g = 1, ng
                  write(ounit, 1210) g, dsigtr(i, g), dsiga(i, g), &
                  dnuf(i, g), dsigf(i, g)
                  group(g) = g
               end do
               write(ounit, *) '  --SCATTERING MATRIX--'
               write(ounit, '(4X, A5, 20I9)') "G/G'", (group(g), g = 1, ng)
               do g = 1, ng
                  write(ounit, 1215) g, (dsigs(i, g, h), h = 1, ng)
               end do
            end do
         end if
         deallocate(bank)
      end if

      1201 format(2x, 'NUMBER OF CONTROL ROD BANK  :', i3)
      1216 format(2x, 'MAX. NUMBER OF STEPS        :', i4)
      1202 format(2x, 'FULLY INSERTED POSITION (cm): ', f4.1, ' (FROM BOTTOM OF THE CORE)')
      1203 format(2x, 'STEP SIZE (cm)              : ', f8.4)
      1204 format(2x, 10(:, 2x, 'Bank ', i2))
      1205 format(10(:, 2x, f7.1), /)
      1209 format(4x, 'MATERIAL', i3)
      1211 format(2x, a7, a12, a12, 2a13)
      1210 format(2x, i6, f13.6, f12.6, 2f13.6)
      1215 format(4x, i3, f14.6, 20f10.6)


      !!! Convert assembly wise CR bank map to node wise CR bank map
      allocate(fbmap(nxx, nyy))
      ytot = 0
      do j = 1, ny
         do ly = 1, ydiv(j)
            ytot = ytot + 1
            xtot = 0
            do i = 1, nx
               do lx = 1, xdiv(i)
                  xtot = xtot + 1
                  fbmap(xtot, ytot) = bmap(i, j)
               end do
            end do
         end do
      end do


      write(ounit, *)
      write(ounit, *) ' ...Control Rods Insertion card is successfully read...'

   end subroutine inp_crod

   !******************************************************************************!

   subroutine inp_ejct (xbunit)

      !
      ! Purpose:
      !    To read rod ejection input

      use sdata, only: nf, ng, lamb, ibeta, velo, nb, tbeta, nmat, &
      ttot, tstep1, tdiv, tstep2, bpos, fbpos, tmove, &
      bspeed, mdir, nstep

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status

      integer :: i, g
      integer :: popt
      integer, dimension(nb) :: bank
      character(len = 4) :: cnb         ! number of bank (character type)

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>>     READING ROD EJECTION DATA      <<<<'
      write(ounit, *) '           --------------------------------------------'

      allocate(tmove(nb), bspeed(nb), mdir(nb), fbpos(nb))
      allocate(tbeta(nmat))

      ! Read Final CR bank position, time to start, and speed
      do i = 1, nb
         read(xbunit, *, iostat = ios) ind, ln, fbpos(i), tmove(i), bspeed(i)
         write (cnb, '(I4)') nb
         cnb = trim(adjustl(cnb))
         message = ' error in reading Final CR Bank Position, time to move and speed for bank : ' // cnb
         call er_message(ounit, ios, ln, message, buf = xbunit)
         if (fbpos(i) > nstep) then
            write(ounit, 1889) ln
            write(*, 1889) ln
            stop
            1889 format(2x, "ERROR AT LINE ", i4, ": WRONG FINAL CONTROL ROD POSITION")
         end if
         if (abs(fbpos(i) - bpos(i)) < 1.e-5_dp) then
            mdir(i) = 0
         else if (fbpos(i) - bpos(i) > 1.e-5_dp) then
            mdir(i) = 2
         else
            mdir(i) = 1
         end if
      end do

      ! Read time for CR to be ejected
      read(xbunit, *, iostat = ios) ind, ln, ttot, tstep1, tdiv, tstep2
      message = ' error in time parameters'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      if (bxtab == 0) then  ! IF XTAB File does not present
         ! Read beta (delayed neutron fraction)
         read(xbunit, *, iostat = ios) ind, ln, (ibeta(i), i = 1, nf)
         message = ' error in reading delayed netron fraction (beta)'
         call er_message(ounit, ios, ln, message, buf = xbunit)

         ! Read precusor decay constant
         read(xbunit, *, iostat = ios) ind, ln, (lamb(i), i = 1, nf)
         message = ' error in reading precusor decay constant'
         call er_message(ounit, ios, ln, message, buf = xbunit)

         ! Read neutron velocity
         allocate(velo(ng))
         read(xbunit, *, iostat = ios) ind, ln, (velo(g), g = 1, ng)
         message = ' error in reading neutron velocity'
         call er_message(ounit, ios, ln, message, buf = xbunit)
      end if


      !! EJCT PRINT OPTION
      read(xbunit, *, iostat = ios) ind, ln, popt
      if (ios == 0 .and. popt > 0) then

         do i = 1, nb
            bank(i) = i
         end do
         write(ounit, 1294) (bank(i), i = 1, nb)
         write(ounit, 1295) (fbpos(i), i = 1, nb)
         write(ounit, 1281) (tmove(i), i = 1, nb)
         write(ounit, 1282) (bspeed(i), i = 1, nb)

         write(ounit, *)
         write(ounit, *) ' TIME PARAMETERS IN SECONDS : '
         write(ounit, 1297) ttot
         write(ounit, 1298) tstep1
         write(ounit, 1299) tstep2
         write(ounit, 1300) tdiv

         if (bxtab == 0) then  ! IF XTAB File does not present
            write(ounit, *)
            write(ounit, *) ' DELAYED NEUTRON FRACTION : '
            write(ounit, '(100F11.5)') (ibeta(i), i = 1, nf)

            write(ounit, *)
            write(ounit, *) ' PRECUSOR DECAY CONSTANT (1/s): '
            write(ounit, '(100F11.5)') (lamb(i), i = 1, nf)

            write(ounit, *)
            write(ounit, *) ' NEUTRON VELOCITY (cm/s) : '
            write(ounit, '(100ES15.5)') (velo(g), g = 1, ng)
         end if
      end if

      ! ttot must be bigger than tstep1 and tstep2
      if ((ttot < tstep1) .or. (ttot < tstep2)) then
         write(ounit, *) 'ERROR: TOTAL SIMULATION TIME SHALL BE GREATER THAN TIME STEPS'
         write(*, *) 'ERROR: TOTAL SIMULATION TIME SHALL BE GREATER THAN TIME STEPS'
         stop
      end if

      ! tdiv must be bigger than tstep1
      if (tdiv < tstep1) then
         write(ounit, *) 'ERROR: THE TIME WHEN SECOND TIME STEP STARTS SHALL BE GREATER THAN FIRST TIME STEP'
         write(*, *) 'ERROR: THE TIME WHEN SECOND TIME STEP STARTS SHALL BE GREATER THAN FIRST TIME STEP'
         stop
      end if

      ! tdiv must be less than ttot
      if (tdiv > ttot) then
         write(ounit, *) 'ERROR: THE TIME WHEN SECOND TIME STEP STARTS SHALL BE LESS THAN TOTAL TIME'
         write(*, *) 'ERROR: THE TIME WHEN SECOND TIME STEP STARTS SHALL BE LESS THAN TOTAL TIME'
         stop
      end if

      ! number of steps shall be less than 10,000
      if (nint(tdiv / tstep1) + nint((ttot - tdiv) / tstep2) > 10000) then
         write(ounit, *) 'ERROR: NUMBER OF TOTAL TIME STEPS ARE MORE THAN 10,000'
         write(*, *) 'ERROR: NUMBER OF TOTAL TIME STEPS ARE MORE THAN 10,000'
         stop
      end if

      write(ounit, *)
      write(ounit, *) ' ...Rod Ejection Card is successfully read...'

      1294 format(25x, 99(:, 2x, 'Bank ', i2))
      1295 format(2x, 'FINAL BANK POS. (STEP)', 99(:, 2x, f7.1), /)
      1281 format(2x, 'STARTS MOVE (SECOND)  ', 99(:, 2x, f7.1), /)
      1282 format(2x, 'SPEED (STEP/SECOND)   ', 99(:, 2x, f7.1), /)
      1297 format(4x, 'TOTAL SIMULATION TIME         : ', f6.2)
      1298 format(4x, 'FIRST TIME STEP               : ', f6.4)
      1299 format(4x, 'SECOND TIME STEP              : ', f6.4)
      1300 format(4x, 'WHEN SECOND TIME STEP APPLY?  : ', f6.2)

   end subroutine inp_ejct

   !******************************************************************************!

   subroutine inp_cbcs (xbunit)

      !
      ! Purpose:
      !    To read boron concentration for critical boron search

      use sdata, only: nmat, ng, rbcon, &
      csigtr, csiga, cnuf, csigf, csigs


      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status

      integer :: i, g, h
      integer :: popt
      integer, dimension(ng) :: group

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>> READING BORON CONCENTRATION FOR BC SEARCH <<<<'
      write(ounit, *) '           --------------------------------------------------'

      ! Read Boron Concentration
      read(xbunit, *, iostat = ios) ind, ln, rbcon
      message = ' error in reading bc guess and bc reference'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      if (bxtab == 0) then  ! IF XTAB File does not present
         allocate(csigtr(nmat, ng))
         allocate(csiga (nmat, ng))
         allocate(cnuf (nmat, ng))
         allocate(csigf (nmat, ng))
         allocate(csigs (nmat, ng, ng))

         ! Read CX changes per ppm born change
         do i = 1, nmat
            do g = 1, ng
               read(xbunit, *, iostat = ios) ind, ln, csigtr(i, g), &
               csiga(i, g), cnuf(i, g), csigf(i, g), (csigs(i, g, h), h = 1, ng)
               message = ' error in reading macro xs changes per ppm boron changes'
               call er_message(ounit, ios, ln, message, buf = xbunit)
            end do
         end do
      end if

      !! BCON PRINT OPTION
      read(xbunit, *, iostat = ios) ind, ln, popt
      if (ios == 0 .and. popt > 0) then

         write(ounit, 1422) rbcon
         if (bxtab == 0) then  ! IF XTAB File does not present
            write(ounit, *)
            write(ounit, *) ' MATERIAL CX CHANGES PER PPM BORON CHANGES : '
            do i = 1, nmat
               write(ounit, 1429) i
               write(ounit, 1431) 'GROUP', 'TRANSPORT', 'ABSORPTION', &
               'NU*FISS', 'FISSION'
               do g = 1, ng
                  write(ounit, 1430) g, csigtr(i, g), csiga(i, g), &
                  cnuf(i, g), csigf(i, g)
                  group(g) = g
               end do
               write(ounit, *) '  --SCATTERING MATRIX--'
               write(ounit, '(4X, A5, 20I11)') "G/G'", (group(g), g = 1, ng)
               do g = 1, ng
                  write(ounit, 1435) g, (csigs(i, g, h), h = 1, ng)
               end do
            end do
         end if
      end if

      1422 format(2x, 'BORON CONCENTRATION REFERENCE :', f8.2)
      1429 format(4x, 'MATERIAL', i3)
      1431 format(2x, a9, a12, a14, a13, a14)
      1430 format(2x, i6, e16.5, 3e14.5)
      1435 format(4x, i3, e17.5, 20e13.5)


      write(ounit, *)
      write(ounit, *) ' ...Critical Boron Search card is successfully read...'

   end subroutine inp_cbcs

   !******************************************************************************!

   subroutine inp_bcon (xbunit)

      !
      ! Purpose:
      !    To read boron concentration

      use sdata, only: nmat, ng, bcon, rbcon, &
      csigtr, csiga, cnuf, csigf, csigs

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status

      integer :: i, g, h
      integer :: popt
      integer, dimension(ng) :: group

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>>       READING BORON CONCENTRATION        <<<<'
      write(ounit, *) '           --------------------------------------------------'

      ! Read Boron Concentration
      read(xbunit, *, iostat = ios) ind, ln, bcon, rbcon
      message = ' error in reading boron concentration and boron concentration reference'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      if (bxtab == 0) then  ! If xtab file present
         ! Read CX changes per ppm born change
         allocate(csigtr(nmat, ng))
         allocate(csiga (nmat, ng))
         allocate(cnuf (nmat, ng))
         allocate(csigf (nmat, ng))
         allocate(csigs (nmat, ng, ng))
         do i = 1, nmat
            do g = 1, ng
               read(xbunit, *, iostat = ios) ind, ln, csigtr(i, g), &
               csiga(i, g), cnuf(i, g), csigf(i, g), (csigs(i, g, h), h = 1, ng)
               message = ' error in reading macro xs changes per ppm boron changes'
               call er_message(ounit, ios, ln, message, buf = xbunit)
            end do
         end do
      end if

      !! BCON PRINT OPTION
      read(xbunit, *, iostat = ios) ind, ln, popt
      if (ios == 0 .and. popt > 0) then

         write(ounit, 1221) bcon
         write(ounit, 1222) rbcon

         if (bxtab == 0) then  ! If xtab file present
            write(ounit, *)
            write(ounit, *) ' MATERIAL CX CHANGES PER PPM BORON CHANGES : '
            do i = 1, nmat
               write(ounit, 1229) i
               write(ounit, 1231) 'GROUP', 'TRANSPORT', 'ABSORPTION', &
               'NU*FISS', 'FISSION'
               do g = 1, ng
                  write(ounit, 1230) g, csigtr(i, g), csiga(i, g), &
                  cnuf(i, g), csigf(i, g)
                  group(g) = g
               end do
               write(ounit, *) '  --SCATTERING MATRIX--'
               write(ounit, '(4X, A5, 20I11)') "G/G'", (group(g), g = 1, ng)
               do g = 1, ng
                  write(ounit, 1235) g, (csigs(i, g, h), h = 1, ng)
               end do
            end do
         end if
      end if

      1221 format(2x, 'BORON CONCENTRATION SET       :', f8.2)
      1222 format(2x, 'BORON CONCENTRATION REFERENCE :', f8.2)
      1229 format(4x, 'MATERIAL', i3)
      1231 format(2x, a9, a12, a14, a13, a14)
      1230 format(2x, i6, e16.5, 3e14.5)
      1235 format(4x, i3, e17.5, 20e13.5)



      write(ounit, *)
      write(ounit, *) ' ...Boron Concentration card is successfully read...'

   end subroutine inp_bcon

   !******************************************************************************!

   subroutine inp_ftem (xbunit)

      !
      ! Purpose:
      !    To read fuel temperature

      use sdata, only: nmat, ng, nnod, ftem, rftem, &
      fsigtr, fsiga, fnuf, fsigf, fsigs

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status

      real(dp) :: cftem
      integer :: i, g, h
      integer :: popt
      integer, dimension(ng) :: group

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>>      READING FUEL TEMPERATURE      <<<<'
      write(ounit, *) '           --------------------------------------------'

      ! Read Fuel Temperature
      read(xbunit, *, iostat = ios) ind, ln, cftem, rftem
      message = ' error in reading average fuel temperature and fuel temperature reference'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      ! ASSIGN CFTEM to FTEM
      if (bther == 0) allocate (ftem(nnod))
      ftem = cftem   !Initial guess for fuel temperature

      if (bxtab == 0) then  ! IF XTAB File does not present
         ! Read CX changes fuel temperature change
         allocate(fsigtr(nmat, ng), fsiga(nmat, ng), fnuf(nmat, ng), fsigf(nmat, ng), fsigs(nmat, ng, ng))
         do i = 1, nmat
            do g = 1, ng
               read(xbunit, *, iostat = ios) ind, ln, fsigtr(i, g), &
               fsiga(i, g), fnuf(i, g), fsigf(i, g), (fsigs(i, g, h), h = 1, ng)
               message = ' error in reading macro xs changes per fuel temperature changes'
               call er_message(ounit, ios, ln, message, buf = xbunit)
            end do
         end do
      end if

      !! FTEM PRINT OPTION
      read(xbunit, *, iostat = ios) ind, ln, popt
      if (ios == 0 .and. popt > 0) then

         if (bther == 0) then
            write(ounit, 1241) cftem
         else
            write(ounit, 1256) cftem
         end if
         write(ounit, 1242) rftem

         if (bxtab == 0) then  ! IF XTAB File does not present
            write(ounit, *)
            write(ounit, *) ' MATERIAL CX CHANGES PER FUEL TEMPERATURE CHANGES : '
            do i = 1, nmat
               write(ounit, 1249) i
               write(ounit, 1251) 'GROUP', 'TRANSPORT', 'ABSORPTION', &
               'NU*FISS', 'FISSION'
               do g = 1, ng
                  write(ounit, 1250) g, fsigtr(i, g), fsiga(i, g), &
                  fnuf(i, g), fsigf(i, g)
                  group(g) = g
               end do
               write(ounit, *) '  --SCATTERING MATRIX--'
               write(ounit, '(4X, A5, 20I11)') "G/G'", (group(g), g = 1, ng)
               do g = 1, ng
                  write(ounit, 1255) g, (fsigs(i, g, h), h = 1, ng)
               end do
            end do
         end if
      end if

      1241 format(2x, 'AVERAGE FUEL TEMPERATURE   :', f6.2)
      1242 format(2x, 'FUEL TEMPERATURE REFERENCE :', f6.2)
      1249 format(4x, 'MATERIAL', i3)
      1251 format(2x, a9, a12, a14, a13, a14)
      1250 format(2x, i6, e16.5, 3e14.5)
      1255 format(4x, i3, e17.5, 20e13.5)
      1256 format(2x, 'AVERAGE FUEL TEMPERATURE   :', f6.2, '  (NOT USED)')



      write(ounit, *)
      write(ounit, *) ' ...Fuel Temperature is card successfully read...'

   end subroutine inp_ftem

   !******************************************************************************!

   subroutine inp_mtem (xbunit)

      !
      ! Purpose:
      !    To read moderator temperature

      use sdata, only: nmat, ng, nnod, mtem, rmtem, &
      msigtr, msiga, mnuf, msigf, msigs

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status

      real(dp) :: cmtem
      integer :: i, g, h
      integer :: popt
      integer, dimension(ng) :: group

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>>   READING MODERATOR TEMPERATURE    <<<<'
      write(ounit, *) '           --------------------------------------------'

      ! Read Moderator Temperature
      read(xbunit, *, iostat = ios) ind, ln, cmtem, rmtem
      message = ' error in reading Moderator temperature and Moderator temperature reference'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      ! ASSIGN CMTEM to MTEM
      if (bther == 0) allocate (mtem(nnod))
      mtem = cmtem

      if (bxtab == 0) then  ! IF XTAB File does not present
         ! Read CX changes per moderator temperature change
         allocate(msigtr(nmat, ng), msiga(nmat, ng), mnuf(nmat, ng), msigf(nmat, ng), msigs(nmat, ng, ng))
         do i = 1, nmat
            do g = 1, ng
               read(xbunit, *, iostat = ios) ind, ln, msigtr(i, g), &
               msiga(i, g), mnuf(i, g), msigf(i, g), (msigs(i, g, h), h = 1, ng)
               message = ' error in reading macro xs changes per Moderator temperature changes'
               call er_message(ounit, ios, ln, message, buf = xbunit)
            end do
         end do
      end if

      !! MTEM PRINT OPTION
      read(xbunit, *, iostat = ios) ind, ln, popt
      if (ios == 0 .and. popt > 0) then

         if (bther == 0) then
            write(ounit, 1261) cmtem
         else
            write(ounit, 1276) cmtem
         end if
         write(ounit, 1262) rmtem
         if (bxtab == 0) then  ! IF XTAB File does not present
            write(ounit, *)
            write(ounit, *) ' MATERIAL CX CHANGES PER MODERATOR TEMPERATURE CHANGES : '
            do i = 1, nmat
               write(ounit, 1269) i
               write(ounit, 1271) 'GROUP', 'TRANSPORT', 'ABSORPTION', &
               'NU*FISS', 'FISSION'
               do g = 1, ng
                  write(ounit, 1270) g, msigtr(i, g), msiga(i, g), &
                  mnuf(i, g), msigf(i, g)
                  group(g) = g
               end do
               write(ounit, *) '  --SCATTERING MATRIX--'
               write(ounit, '(4X, A5, 20I11)') "G/G'", (group(g), g = 1, ng)
               do g = 1, ng
                  write(ounit, 1275) g, (msigs(i, g, h), h = 1, ng)
               end do
            end do
         end if
      end if

      1261 format(2x, 'AVERAGE MODERATOR TEMPERATURE   :', f6.2)
      1262 format(2x, 'MODERATOR TEMPERATURE REFERENCE :', f6.2)
      1269 format(4x, 'MATERIAL', i3)
      1271 format(2x, a9, a12, a14, a13, a14)
      1270 format(2x, i6, e16.5, 3e14.5)
      1275 format(4x, i3, e17.5, 20e13.5)
      1276 format(2x, 'AVERAGE MODERATOR TEMPERATURE   :', f6.2, '  (NOT USED)')



      write(ounit, *)
      write(ounit, *) ' ...Moderator Temperature Card is successfully read...'


   end subroutine inp_mtem

   !******************************************************************************!

   subroutine inp_cden (xbunit)

      !
      ! Purpose:
      !    To read Coolant Density

      use sdata, only: nmat, ng, nnod, cden, rcden, &
      lsigtr, lsiga, lnuf, lsigf, lsigs

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status

      real(dp) :: ccden
      integer :: i, g, h
      integer :: popt
      integer, dimension(ng) :: group

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>>       READING COOLANT DENSITY      <<<<'
      write(ounit, *) '           --------------------------------------------'

      ! Read Coolant Density
      read(xbunit, *, iostat = ios) ind, ln, ccden, rcden
      message = ' error in reading Coolant Density and Coolant Density reference'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      !ASSIGN CCDEN TO CDEN
      if (bther == 0) allocate (cden(nnod))
      cden = ccden

      if (bxtab == 0) then  ! IF XTAB File does not present
         ! Read CX changes per Coolant Density change
         allocate(lsigtr(nmat, ng), lsiga(nmat, ng), lnuf(nmat, ng), lsigf(nmat, ng), lsigs(nmat, ng, ng))
         do i = 1, nmat
            do g = 1, ng
               read(xbunit, *, iostat = ios) ind, ln, lsigtr(i, g), &
               lsiga(i, g), lnuf(i, g), lsigf(i, g), (lsigs(i, g, h), h = 1, ng)
               message = ' error in reading macro xs changes per Coolant Density changes'
               call er_message(ounit, ios, ln, message, buf = xbunit)
            end do
         end do
      end if

      !! CDEN PRINT OPTION
      read(xbunit, *, iostat = ios) ind, ln, popt
      if (ios == 0 .and. popt > 0) then

         if (bther == 0) then
            write(ounit, 1361) ccden
         else
            write(ounit, 1376) ccden
         end if
         write(ounit, 1362) rcden
         if (bxtab == 0) then  ! IF XTAB File does not present
            write(ounit, *)
            write(ounit, *) ' MATERIAL CX CHANGES PER COOLANT DENSITY CHANGES : '
            do i = 1, nmat
               write(ounit, 1369) i
               write(ounit, 1371) 'GROUP', 'TRANSPORT', 'ABSORPTION', &
               'NU*FISS', 'FISSION'
               do g = 1, ng
                  write(ounit, 1370) g, lsigtr(i, g), lsiga(i, g), &
                  lnuf(i, g), lsigf(i, g)
                  group(g) = g
               end do
               write(ounit, *) '  --SCATTERING MATRIX--'
               write(ounit, '(4X, A5, 20I11)') "G/G'", (group(g), g = 1, ng)
               do g = 1, ng
                  write(ounit, 1375) g, (lsigs(i, g, h), h = 1, ng)
               end do
            end do
         end if
      end if


      1361 format(2x, 'AVERAGE COOLANT DENSITY   :', f8.4)
      1362 format(2x, 'COOLANT DENSITY REFERENCE :', f8.4)
      1369 format(4x, 'MATERIAL', i3)
      1371 format(2x, a9, a12, a14, a13, a14)
      1370 format(2x, i6, e16.5, 3e14.5)
      1375 format(4x, i3, e17.5, 20e13.5)
      1376 format(2x, 'AVERAGE COOLANT DENSITY   :', f8.4, '  (USED AS GUESS)')



      write(ounit, *)
      write(ounit, *) ' ...Coolant Density Card is successfully read...'


   end subroutine inp_cden

   !******************************************************************************!

   subroutine inp_ther (xbunit)

      !
      ! Purpose:
      !    To read thermalhydraulics parameters input

      use sdata, only: pow, tin, nx, ny, nxx, nyy, ystag, &
      rf, tg, tc, rg, rc, ppitch, cf, dia, cflow, dh, pi, &
      farea, xdiv, ydiv, ystag, node_nf, nm, nt, rdel, rpos, &
      nnod, tfm, ppow, ent, heatf, ntem, stab, ftem, mtem, cden, &
      thunit

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln   !Line number
      integer :: ios  ! IOSTAT status

      integer :: i, j
      integer :: nfpin, ngt                              ! Number of fuel pin and guide tubes

      integer :: ly, lx, ytot, xtot
      real(dp) :: cmflow
      real(dp), dimension(nx, ny) :: area
      real(dp) :: barea, div

      real(dp) :: dum

      integer :: popt

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>>   READING THERMAL-HYDRAULIC DATA   <<<<'
      write(ounit, *) '           --------------------------------------------'

      allocate (ftem(nnod), mtem(nnod), cden(nnod))

      ! Read Percent Power
      read(xbunit, *, iostat = ios) ind, ln, ppow
      message = ' error in reading percent power'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      ! Read reactor Power
      read(xbunit, *, iostat = ios) ind, ln, pow
      message = ' error in reading reactor full thermal power'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      ! Read inlet coolant temp. (Kelvin) and  FA flow rate (kg/s)
      read(xbunit, *, iostat = ios) ind, ln, tin, cmflow
      message = ' error in reading coolant inlet temp. and Fuel Assembly mass flow rate'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      ! Read fuel pin geometry in meter
      read(xbunit, *, iostat = ios) ind, ln, rf, tg, tc, ppitch
      message = ' error in reading fuel meat rad., gap thickness, clad thickness and pin pitch'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      ! Check gap and clad thickness
      if (tg > 0.25 * rf) then
         write(ounit, *) '  ERROR: GAP THICKNESS IS TO LARGE (> 0.25*rf)'
         stop
      end if
      if (tc > 0.25 * rf) then
         write(ounit, *) '  ERROR: CLADDING THICKNESS IS TO LARGE (> 0.25*rf)'
         stop
      end if

      ! Read Number of fuel pins and guide tubes
      read(xbunit, *, iostat = ios) ind, ln, nfpin, ngt
      message = ' error in reading number of fuel pins and guide tubes'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      ! Read Fraction of heat deposited in the coolant
      read(xbunit, *, iostat = ios) ind, ln, cf
      message = ' error in reading fraction of heat deposited in the coolant'
      call er_message(ounit, ios, ln, message, buf = xbunit)

      if (cf < 0. .or. cf > 1.0) stop "The value of the fraction of heat " &
      // "deposited in the coolant is incorrect"

      ! Calculate outer radius of gap and cladding
      rg = rf + tg
      rc = rg + tc

      ! Calculate fuel pin diameter
      dia = 2. * rc

      ! Calculate hydraulic diameter
      dh = dia * ((4. / pi) * (ppitch / dia)**2 - 1.)

      ! Calculate sub-channel area
      farea = ppitch**2 - 0.25 * pi * dia**2

      ! Calculate sub-channel mass flow rate
      cflow = cmflow / real(nfpin)

      ! Calculate total coolant mass flow rate and number of fuel pins per node
      barea = 0.
      do j = 1, ny
         do i = 1, nx
            area(i, j) = xsize(i) * ysize(j)             ! assembly area
            if (area(i, j) > barea) barea = area(i, j)  ! barea => largest assembly area for ref.
         end do
      end do

      allocate(node_nf(nxx, nyy))
      node_nf = 0.

      ytot = 0
      do j = 1, ny
         do ly = 1, ydiv(j)
            ytot = ytot + 1
            xtot = 0
            do i = 1, nx
               do lx = 1, xdiv(i)
                  xtot = xtot + 1
                  if ((xtot >= ystag(ytot)%smin) .and. (xtot <= ystag(ytot)%smax)) then
                     div = real (ydiv(j) * xdiv(i))            ! Number of nodes in current assembly
                     node_nf(xtot, ytot) = area(i, j) * real(nfpin) / (barea * div)   ! Number of fuel pin for this node
                  end if
               end do
            end do
         end do
      end do

      ! Calculate fuel pin mesh delta and position
      dum = rf / real(nm)

      allocate(rdel(nt), rpos(nt))
      do i = 1, nm
         rdel(i) = dum
      end do

      ! Fuel pin mesh size
      rdel(nm + 1) = tg
      rdel(nm + 2) = tc

      ! Fuel pin mesh position
      rpos(1) = 0.5 * rdel(1)
      do i = 2, nt
         rpos(i) = rpos(i - 1) + 0.5 * (rdel(i - 1) + rdel(i))
      end do

      ! Guess fuel and moderator temperature
      allocate(tfm(nnod, nt + 1)) ! Allocate fuel pin mesh temperature
      tfm = 900.                !Initial guess for radial fuel pin temperature distribution
      if (bxtab == 1) then
         ftem = 900. ;
         cden = 0.711;
         mtem = 500.
      end if

      allocate(ent(nnod))
      allocate(heatf(nnod))

      ! Initial heat-flux rate
      heatf = 0.

      ! Open steam table file
      open (unit = thunit, file = 'st155bar', status = 'OLD', action = 'READ', iostat = ios)
      if (ios == 0) then  !If steam table present in the current directory
         ! Save steam table data to stab
         do i = 1, ntem
            read(thunit, *) stab(i, 1), stab(i, 2), stab(i, 3), &
            stab(i, 4), stab(i, 5), stab(i, 6)
         end do
         close(unit = thunit)
      else  ! if steam table not present, use steam table data at 15.5 MPa
         stab(1, 1) = 543.15;
         stab(1, 2) = 0.78106745;
         stab(1, 3) = 1182595.0;
         stab(1, 4) = 0.820773;
         stab(1, 5) = 0.128988;
         stab(1, 6) = 0.60720
         stab(2, 1) = 553.15;
         stab(2, 2) = 0.76428125;
         stab(2, 3) = 1232620.0;
         stab(2, 4) = 0.829727;
         stab(2, 5) = 0.126380;
         stab(2, 6) = 0.59285
         stab(3, 1) = 563.15;
         stab(3, 2) = 0.74619716;
         stab(3, 3) = 1284170.0;
         stab(3, 4) = 0.846662;
         stab(3, 5) = 0.124118;
         stab(3, 6) = 0.57710
         stab(4, 1) = 573.15;
         stab(4, 2) = 0.72650785;
         stab(4, 3) = 1337630.0;
         stab(4, 4) = 0.863597;
         stab(4, 5) = 0.121856;
         stab(4, 6) = 0.55970
         stab(5, 1) = 583.15;
         stab(5, 2) = 0.70475081;
         stab(5, 3) = 1393570.0;
         stab(5, 4) = 0.915035;
         stab(5, 5) = 0.120105;
         stab(5, 6) = 0.54045
         stab(6, 1) = 593.15;
         stab(6, 2) = 0.68018488;
         stab(6, 3) = 1452895.0;
         stab(6, 4) = 0.966472;
         stab(6, 5) = 0.118354;
         stab(6, 6) = 0.51880
         stab(7, 1) = 603.15;
         stab(7, 2) = 0.65150307;
         stab(7, 3) = 1517175.0;
         stab(7, 4) = 1.166745;
         stab(7, 5) = 0.143630;
         stab(7, 6) = 0.49420
         stab(8, 1) = 613.15;
         stab(8, 2) = 0.61590149;
         stab(8, 3) = 1589770.0;
         stab(8, 4) = 1.515852;
         stab(8, 5) = 0.195931;
         stab(8, 6) = 0.46550
         stab(9, 1) = 617.91;
         stab(9, 2) = 0.59896404;
         stab(9, 3) = 1624307.1;
         stab(9, 4) = 1.681940;
         stab(9, 5) = 0.220813;
         stab(9, 6) = 0.45185
      end if

      !! THER PRINT OPTION
      read(xbunit, *, iostat = ios) ind, ln, popt
      if (ios == 0 .and. popt > 0) then
         write(ounit, 1309) ppow
         write(ounit, 1301) pow
         write(ounit, 1302) tin
         write(ounit, 1303) cmflow
         write(ounit, 1304) rf
         write(ounit, 1305) tg
         write(ounit, 1306) tc
         write(ounit, 1310) ppitch
         write(ounit, 1307) cf
      end if

      write(ounit, *)
      write(ounit, *) ' ...Thermal-hydraulic Card is successfully read...'

      1309 format(2x, 'REACTOR PERCENT POWER (%)                : ', f12.5)
      1301 format(2x, 'REACTOR POWER (Watt)                     : ', es12.4)
      1302 format(2x, 'COOLANT INLET TEMPERATURE (Kelvin)       : ', es12.4)
      1303 format(2x, 'FUEL ASSEMBLY MASS FLOW RATE (Kg/s)      : ', es12.4)
      1304 format(2x, 'FUEL MEAT RADIUS (m)                     : ', es12.4)
      1305 format(2x, 'GAP THICKNESS (m)                        : ', es12.4)
      1306 format(2x, 'CLAD THICKNESS (m)                       : ', es12.4)
      1310 format(2x, 'PIN PITCH(m)                             : ', es12.4)
      1307 format(2x, 'FRACTION OF HEAT DEPOSITED IN COOL.      : ', es12.4)

      deallocate(xsize, ysize, zsize)

   end subroutine inp_ther

   !******************************************************************************!

   subroutine er_message (funit, iost, ln, mess, xtab, buf)
      !
      ! Purpose:
      !    To provide error message
      !

      implicit none

      integer, intent(in) :: funit, iost, ln
      character(len = *), intent(in) :: mess
      integer, optional, intent(in) :: xtab, buf

      if (iost < 0) then
         write(funit, *)
         write(*, *)
         write(*, *) '' // achar(27) // '[31m OOPS! WE FOUND AN ERROR.' // achar(27) // '[0m'

         if (present(xtab)) then
            write(funit, 1014) ln, xtab
         else
            write(funit, 1006) carr(getloc(uarr, buf))
            write(funit, 1013) ln, farr(getloc(uarr, buf))
         end if
         write(funit, *) mess
         if (present(xtab)) then
            write(*, 1014) ln, xtab
         else
            write(*, 1006) carr(getloc(uarr, buf))
            write(*, 1013) ln, farr(getloc(uarr, buf))
         end if
         write(*, *) mess
         1013 format(2x, 'THIS LINE NEEDS MORE INPUT DATA. LINE', i4, &
         ' IN FILE : ', a100)
         1014 format(2x, 'ERROR: LINE', i4, &
         'IN XTAB FILE FOR MATERIAL NUMBER', i4, '. IT NEEDS MORE DATA')
         stop
      end if

      if (iost > 0) then
         write(funit, *)
         write(*, *)
         write(*, *) '' // achar(27) // '[31m OOPS! WE FOUND AN ERROR.' // achar(27) // '[0m'

         if (present(xtab)) then
            write(funit, 1005) ln, xtab
         else
            write(funit, 1006) carr(getloc(uarr, buf))
            write(funit, 1004) ln, farr(getloc(uarr, buf))
         end if
         write(funit, *) mess
         if (present(xtab)) then
            write(*, 1005) ln, xtab
         else
            write(*, 1006) carr(getloc(uarr, buf))
            write(*, 1004) ln, farr(getloc(uarr, buf))
         end if
         write(*, *) mess
         1006 format(2x, 'ERROR: THERE IS AN ERROR IN CARD %', a4)
         1004 format(2x, 'PLEASE CHECK LINE NUMBER', i4, ' IN FILE : ', a100)
         1005 format(2x, 'ERROR: PLEASE CHECK LINE NUMBER', i4, &
         ' IN XTAB FILE FOR MATERIAL NUMBER ', i4)
         stop
      end if

   end subroutine er_message

   !****************************************************************************!

   function getloc(arr, elm) result (loc)

      !Purpose: To  get location of an element in an array

      integer, dimension(:), intent(in) :: arr    ! input array
      integer, intent(in) :: elm                  ! element whose location is wanted to find
      integer :: loc

      integer :: i

      do i = 1, size(arr)
         if (arr(i) == elm) then
            loc = i
            exit
         end if
      end do

   end function getloc

   !******************************************************************************!

   subroutine print_outp(fn)

      !
      ! Purpose:
      !    To print detailed output
      !

      use sdata, only: nx, ny, nxx, nyy, nzz, zdel, &
      nnod, ix, iy, iz, &
      ftem, mtem, cden, tfm, pow, ppow, vdel

      implicit none

      real(dp), dimension(:), intent(in) :: fn

      real(dp), dimension(nxx, nyy, nzz) :: fx
      integer :: i, j, k, n
      real(dp), dimension(:, :, :), allocatable :: fdist
      real(dp) :: totp
      real(dp) :: avg
      character(len = 12) :: avg_char
      integer :: nfuel
      integer :: nxc, nyc
      character(len = 6), dimension(:, :, :), allocatable :: cpow
      logical, dimension(nzz) :: is_print_planar

      fx = 0._dp
      do n = 1, nnod
         fx(ix(n), iy(n), iz(n)) = fn(n)
      end do

      call get_distribution(fx, fdist, is_power = .true., is_max = .false.)

      !NORMALIZE
      if (output_opt % assembly) then
         nxc = nx
         nyc = ny
      else
         nxc = nxx
         nyc = nyy
      end if

      nfuel = 0
      totp = 0._dp
      do k = 1, nzz
         do j = 1, nyc
            do i = 1, nxc
               if (fdist(i, j, k) > 0) then
                  nfuel = nfuel + 1
                  totp = totp + fdist(i, j, k)
               end if
            end do
         end do
      end do
      do k = 1, nzz
         do j = 1, nyc
            do i = 1, nxc
               fdist(i, j, k) = fdist(i, j, k) * nfuel / totp
            end do
         end do
      end do

      ! Get average power density
      if (bther == 0) then
         avg = 1.0
      else
         avg = pow * ppow * 1.e-5_dp / get_active_volume(fn)  ! Power density (kw/cm3)
         avg = avg * 1.0e6  ! to kw/m3
      end if
      write (avg_char, '(ES12.5)') avg

      open (unit = 102, file = trim(iname) // '_3d_power.out', status = 'REPLACE', action = 'WRITE')
      call dist_to_char(fdist, 'F6.4', is_print_planar, cpow)
      if (bther == 0) then
         call print_3d(102, '   3-D Power Distribution (Relative)', is_print_planar, cpow)
      else
         call print_3d(102, '   3-D Power Distribution (Relative. Average =' &
         // avg_char // ' kw/m3)', is_print_planar, cpow)
      end if
      close(102)

      if (bther == 1) then
         fx = 0._dp
         do n = 1, nnod
            fx(ix(n), iy(n), iz(n)) = ftem(n)
         end do
         call get_distribution(fx, fdist, is_power = .false., is_max = .false.)
         call dist_to_char(fdist, 'F6.1', is_print_planar, cpow)
         open (unit = 102, file = trim(iname) // '_fuel_temp.out', status = 'REPLACE', action = 'WRITE')
         call print_3d(102, '   3D Fuel Temperature Distribution (K)', is_print_planar, cpow)
         close(102)

         fx = 0._dp
         do n = 1, nnod
            fx(ix(n), iy(n), iz(n)) = mtem(n)
         end do
         call get_distribution(fx, fdist, is_power = .false., is_max = .false.)
         call dist_to_char(fdist, 'F6.1', is_print_planar, cpow)
         open (unit = 102, file = trim(iname) // '_cool_temp.out', status = 'REPLACE', action = 'WRITE')
         call print_3d(102, '   3D Coolant Temperature Distribution (K)', is_print_planar, cpow)
         close(102)

         fx = 0._dp
         do n = 1, nnod
            fx(ix(n), iy(n), iz(n)) = cden(n)
         end do
         call get_distribution(fx, fdist, is_power = .false., is_max = .false.)
         call dist_to_char(fdist, 'F6.3', is_print_planar, cpow)
         open (unit = 102, file = trim(iname) // '_cool_dens.out', status = 'REPLACE', action = 'WRITE')
         call print_3d(102, '   3D Coolant Density Distribution (g/cm3)', is_print_planar, cpow)
         close(102)

         fx = 0._dp
         do n = 1, nnod
            fx(ix(n), iy(n), iz(n)) = tfm(n, 1)
         end do
         call get_distribution(fx, fdist, is_power = .false., is_max = .true.)
         call dist_to_char(fdist, 'F6.1', is_print_planar, cpow)
         open (unit = 102, file = trim(iname) // '_max_centerline.out', status = 'REPLACE', action = 'WRITE')
         call print_3d(102, '   3D Max Centerline Fuel Temperature (K)', is_print_planar, cpow)
         close(102)
      end if

   end subroutine print_outp

   !******************************************************************************!

   function get_active_volume(fn) result(vol)

      use sdata, only: nnod, vdel

      implicit none

      real(dp), dimension(:), intent(in) :: fn
      real(dp) :: vol

      integer :: n

      vol = 0.
      do n = 1, nnod
         if (fn(n) > 0.0) vol = vol + vdel(n)
      end do

   end function get_active_volume

   !******************************************************************************!

   subroutine get_distribution(node_dist, fdist, is_power, is_max)

      !
      ! Purpose:
      !    convert node-wise distribution to assembly-wise or node-wise distribution
      !

      use sdata, only: nx, ny, nxx, nyy, nzz

      implicit none


      real(dp), dimension(:, :, :), intent(in) :: node_dist
      logical, intent(in) :: is_power
      logical, intent(in) :: is_max
      real(dp), dimension(:, :, :), allocatable, intent(out) :: fdist

      if (output_opt % assembly) then
         allocate(fdist(nx, ny, nzz))
         call node_to_asm_dist(node_dist, fdist, is_power, is_max)
      else
         allocate(fdist(nxx, nyy, nzz))
         fdist = node_dist
      end if


   end subroutine get_distribution

   !******************************************************************************!

   subroutine node_to_asm_dist(node_dist, asm_dist, is_power, is_max)

      !
      ! Purpose:
      !    convert node-wise distribution to assembly-wise distribution
      !

      use sdata, only: nx, ny, nzz, ydiv, xdiv, xdel, ydel, zdel

      implicit none


      real(dp), dimension(:, :, :), intent(in) :: node_dist
      logical, intent(in) :: is_power
      logical, intent(in) :: is_max
      real(dp), dimension(:, :, :), intent(out) :: asm_dist

      integer :: i, j, k
      integer :: ly, lx, ys, xs, yf, xf
      real(dp) :: tmp, vsumm

      do k = 1, nzz
         ys = 1
         yf = 0
         do j = 1, ny
            yf = yf + ydiv(j)
            xf = 0
            xs = 1
            do i = 1, nx
               xf = xf + xdiv(i)
               tmp = 0._dp
               vsumm = 0._dp

               do ly = ys, yf
                  do lx = xs, xf
                     ! Get assembly average
                     if (.not. is_max) then
                        if (is_power) then
                           tmp = tmp + node_dist(lx, ly, k)
                        else
                           tmp = tmp + node_dist(lx, ly, k) * xdel(lx) * ydel(ly) * zdel(k)
                        end if
                        vsumm = vsumm + xdel(lx) * ydel(ly) * zdel(k)
                     else
                        ! Get max value
                        if (tmp < node_dist(lx, ly, k)) tmp = node_dist(lx, ly, k)
                        vsumm = 1.0
                     end if
                  end do
               end do
               asm_dist(i, j, k) = tmp / vsumm
               xs = xs + xdiv(i)
            end do
            ys = ys + ydiv(j)
         end do
      end do


   end subroutine node_to_asm_dist

   !******************************************************************************!

   subroutine dist_to_char(val, format, is_print_planar, cpow)

      !
      ! Purpose:
      !    convert distribution to character
      !

      use sdata, only: nx, ny, nxx, nyy, nzz

      implicit none


      real(dp), dimension(:, :, :), intent(in) :: val
      character(len = *), intent(in) :: format
      logical, dimension(:), intent(out) :: is_print_planar
      character(len = 6), dimension(:, :, :), allocatable, intent(out) :: cpow
      integer :: i, j, k

      if (output_opt % assembly) then
         allocate(cpow(nx, ny, nzz))
      else
         allocate(cpow(nxx, nyy, nzz))
      end if

      ! CONVERT TO CHARACTER
      is_print_planar = .false.
      do k = 1, size(val, dim = 3)
         if (sum(val(:, :, k)) > 0.0) is_print_planar(k) = .true.
         do j = 1, size(val, dim = 2)
            do i = 1, size(val, dim = 1)
               ! Convert power to character (If power == 0 convert to blank spaces)
               if ((val(i, j, k) - 0.) < 1.e-5_dp) then
                  cpow(i, j, k) = '     '
               else
                  write (cpow(i, j, k), '(' // trim(format) // ')') val(i, j, k)
                  cpow(i, j, k) = trim(cpow(i, j, k))
               end if
            end do
         end do
      end do


   end subroutine dist_to_char

   !******************************************************************************!

   subroutine print_3d(funit, msg, is_print_planar, cpow)

      !
      ! Purpose:
      !    To print 3d distribution
      !

      use sdata, only: nx, ny, nxx, nyy, nzz

      implicit none

      integer, intent(in) :: funit
      character(len = *), intent(in) :: msg
      logical, dimension(:), intent(in) :: is_print_planar
      character(len = 6), dimension(:, :, :), intent(in) :: cpow
      integer :: n, i, j, k
      integer :: xs, xf

      integer, parameter :: xm = 20
      integer :: ip, ipr
      integer :: nxc, nyc

      if (output_opt % assembly) then
         nxc = nx
         nyc = ny
      else
         nxc = nxx
         nyc = nyy
      end if


      ! Print assembly power distribution
      write(funit, *) '   ' // trim(adjustl(msg))
      write(funit, *) '  =============================='

      do k = nzz, 1, -1
         if (is_print_planar(k)) then
            ip = nxc / xm
            ipr = mod(nxc, xm) - 1
            xs = 1;
            xf = xm
            write(funit, 100) k
            do n = 1, ip
               write(funit, '(4X,100I8)') (i, i = xs, xf)
               do j = nyc, 1, -1
                  write(funit, '(2X,I4,100A8)') j, (cpow(i, j, k), i = xs, xf)
               end do
               write(funit, *)
               xs = xs + xm
               xf = xf + xm
            end do


            write(funit, '(4X,100I8)') (i, i = xs, xs + ipr)
            if (xs + ipr > xs) then
               do j = nyc, 1, -1
                  write(funit, '(2X,I4,100A8)') j, (cpow(i, j, k), i = xs, xs + ipr)
               end do
            end if
         end if
      end do

      write(funit, *)

      100 format(2x, 'z = ', i3)


   end subroutine print_3d

   !******************************************************************************!

   subroutine asmpow(fn)

      !
      ! Purpose:
      !    To print axially averaged assembly-wise power distribution
      !

      use sdata, only: nx, ny, nxx, nyy, nzz, zdel, &
      xdel, ydel, ystag, nnod, ix, iy, iz, &
      xdiv, ydiv

      implicit none

      real(dp), dimension(:), intent(in) :: fn

      real(dp), dimension(nxx, nyy, nzz) :: fx
      integer :: i, j, k, n
      integer :: ly, lx, ys, xs, yf, xf
      real(dp) :: summ, vsumm
      real(dp), dimension(nxx, nyy) :: fnode
      real(dp), dimension(nx, ny) :: fasm
      real(dp) :: totp
      integer :: nfuel
      real(dp) :: fmax
      integer :: xmax, ymax
      character(len = 6), dimension(nx, ny) :: cpow

      integer, parameter :: xm = 20
      integer :: ip, ipr

      fx = 0._dp
      do n = 1, nnod
         fx(ix(n), iy(n), iz(n)) = fn(n)
      end do

      !Calculate axially averaged node-wise distribution
      fnode = 0._dp
      do j = 1, nyy
         do i = ystag(j)%smin, ystag(j)%smax
            summ = 0._dp
            vsumm = 0._dp
            do k = 1, nzz
               summ = summ + fx(i, j, k) * zdel(k)
               vsumm = vsumm + zdel(k)
            end do
            fnode(i, j) = summ / vsumm
         end do
      end do

      !Calculate assembly power
      nfuel = 0
      totp = 0._dp
      ys = 1
      yf = 0
      do j = 1, ny
         yf = yf + ydiv(j)
         xf = 0
         xs = 1
         do i = 1, nx
            xf = xf + xdiv(i)
            summ = 0._dp
            vsumm = 0._dp
            do ly = ys, yf
               do lx = xs, xf
                  summ = summ + fnode(lx, ly) * xdel(lx) * ydel(ly)
                  vsumm = vsumm + xdel(lx) * ydel(ly)
               end do
            end do
            fasm(i, j) = summ / vsumm
            xs = xs + xdiv(i)
            if (fasm(i, j) > 0._dp) nfuel = nfuel + 1
            if (fasm(i, j) > 0._dp) totp = totp + fasm(i, j)
         end do
         ys = ys + ydiv(j)
      end do


      ! Normalize assembly power to 1._DP
      xmax = 1;
      ymax = 1
      fmax = 0._dp
      do j = 1, ny
         do i = 1, nx
            if (totp > 0.) fasm(i, j) = real(nfuel) / totp * fasm(i, j)
            if (fasm(i, j) > fmax) then     ! Get max position
               xmax = i
               ymax = j
               fmax = fasm(i, j)
            end if
            ! Convert power to character (If power == 0 convert to blank spaces)
            if ((fasm(i, j) - 0.) < 1.e-5_dp) then
               cpow(i, j) = '     '
            else
               write (cpow(i, j), '(F6.3)') fasm(i, j)
               cpow(i, j) = trim(cpow(i, j))
            end if
         end do
      end do


      ! Print assembly power distribution
      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '    Radial Power Distribution'
      write(ounit, *) '  =============================='

      ip = nx / xm
      ipr = mod(nx, xm) - 1
      xs = 1;
      xf = xm
      do k = 1, ip
         write(ounit, '(4X,100I8)') (i, i = xs, xf)
         do j = ny, 1, -1
            write(ounit, '(2X,I4,100A8)') j, (cpow(i, j), i = xs, xf)
         end do
         write(ounit, *)
         xs = xs + xm
         xf = xf + xm
      end do


      write(ounit, '(4X,100I8)') (i, i = xs, xs + ipr)
      if (xs + ipr > xs) then
         do j = ny, 1, -1
            write(ounit, '(2X,I4,100A8)') j, (cpow(i, j), i = xs, xs + ipr)
         end do
      end if



      write(ounit, *)

      write(ounit, *) '  MAX POS.       Maximum Value'
      write(ounit, 1101) ymax, xmax, fasm(xmax, ymax)

      1101 format(2x, '(', i3, ',', i3, ')', f15.3)


   end subroutine asmpow

   !******************************************************************************!

   subroutine axipow(fn)

      !
      ! Purpose:
      !    To print radially averaged  power distribution
      !

      use sdata, only: nxx, nyy, nzz, nz, zdiv, &
      vdel, ystag, nnod, ix, iy, iz, xyz, zdel, ystag, coreh

      implicit none

      real(dp), dimension(:), intent(in) :: fn

      real(dp), dimension(nxx, nyy, nzz) :: fx
      integer :: i, j, k, n, ztot
      integer :: lz
      real(dp) :: summ, vsumm
      real(dp), dimension(nz) :: faxi
      real(dp) :: totp
      integer :: nfuel
      real(dp) :: fmax
      integer :: amax

      fx = 0._dp
      do n = 1, nnod
         fx(ix(n), iy(n), iz(n)) = fn(n)
      end do

      ! Calculate Axial Power
      nfuel = 0
      totp = 0._dp
      ztot = 0
      do k = 1, nz
         summ = 0._dp
         vsumm = 0._dp
         do lz = 1, zdiv(k)
            ztot = ztot + 1
            do j = 1, nyy
               do i = ystag(j)%smin, ystag(j)%smax
                  summ = summ + fx(i, j, ztot)
                  vsumm = vsumm + vdel(xyz(i, j, ztot))
               end do
            end do
         end do
         faxi(k) = summ / vsumm
         if (faxi(k) > 0._dp) nfuel = nfuel + 1
         if (faxi(k) > 0._dp) totp = totp + faxi(k)
      end do

      ! Normalize Axial power to 1._DP
      fmax = 0._dp
      amax = 1
      do k = 1, nz
         faxi(k) = real(nfuel) / totp * faxi(k)
         if (faxi(k) > fmax) then
            amax = k   ! Get max position
            fmax = faxi(k)
         end if
      end do

      ! Print Axial power distribution
      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '    Axial Power Density Distribution'
      write(ounit, *) '  ===================================='
      write(ounit, *)
      write(ounit, *) '    Plane Number        Power      Height'
      write(ounit, *) '   -----------------------------------------'
      summ = 0.
      ztot = nzz
      do k = nz, 1, -1
         if (k == nz) then
            write(ounit, '(2X,I8,A7,F13.3, F12.2)') k, ' (TOP)', faxi(k), coreh - summ
         else if (k == 1) then
            write(ounit, '(2X,I8,A10,F10.3, F12.2)') k, ' (BOTTOM)', faxi(k), coreh - summ
         else
            write(ounit, '(2X,I8,F20.3, F12.2)') k, faxi(k), coreh - summ
         end if
         do lz = 1, zdiv(k)
            summ = summ + zdel(ztot)
            ztot = ztot - 1
         end do
      end do
      write(ounit, *)
      write(ounit, *) '  MAX POS.       Maximum Value'
      write(ounit, 1102) amax, faxi(amax)

      1102 format(4x, '(', i3, ')', f18.3)


   end subroutine axipow

   !******************************************************************************!

   subroutine asmflux(norm)

      !
      ! Purpose:
      !    To print axially averaged assembly-wise flux distribution
      !

      use sdata, only: f0, ng, nx, ny, nxx, nyy, nzz, zdel, &
      xdel, ydel, ystag, nnod, ix, iy, iz, &
      xdiv, ydiv

      implicit none

      real(dp), optional, intent(in) :: norm

      real(dp), dimension(nxx, nyy, nzz, ng) :: fx
      integer :: g, i, j, k, n
      integer :: ly, lx, ys, xs, yf, xf
      real(dp) :: summ, vsumm
      real(dp), dimension(nxx, nyy, ng) :: fnode
      real(dp), dimension(nnod, ng) :: flux
      real(dp), dimension(nx, ny, ng) :: fasm
      real(dp), dimension(ng) :: totp
      character(len = 10), dimension(nx, ny, ng) :: cflx

      integer, parameter :: xm = 20
      integer :: ip, ipr
      integer :: negf

      if (bther == 1) then
         call flux_atpower(flux)  !calculate actual flux at given power level
      else
         flux = f0
      end if

      fx = 0._dp
      do g = 1, ng
         do n = 1, nnod
            fx(ix(n), iy(n), iz(n), g) = flux(n, g)
         end do
      end do

      !Calculate axially averaged node-wise distribution
      fnode = 0._dp
      do g = 1, ng
         do j = 1, nyy
            do i = ystag(j)%smin, ystag(j)%smax
               summ = 0._dp
               vsumm = 0._dp
               do k = 1, nzz
                  summ = summ + fx(i, j, k, g) * zdel(k)
                  vsumm = vsumm + zdel(k)
               end do
               fnode(i, j, g) = summ / vsumm
            end do
         end do
      end do

      !Calculate Radial Flux (assembly wise)
      negf = 0
      do g = 1, ng
         totp(g) = 0._dp
         ys = 1
         yf = 0
         do j = 1, ny
            yf = yf + ydiv(j)
            xf = 0
            xs = 1
            do i = 1, nx
               xf = xf + xdiv(i)
               summ = 0._dp
               vsumm = 0._dp
               do ly = ys, yf
                  do lx = xs, xf
                     summ = summ + fnode(lx, ly, g) * xdel(lx) * ydel(ly)
                     vsumm = vsumm + xdel(lx) * ydel(ly)
                  end do
               end do
               fasm(i, j, g) = summ / vsumm
               xs = xs + xdiv(i)
               if (fasm(i, j, g) > 0._dp) then
                  totp(g) = totp(g) + fasm(i, j, g)
               end if
               ! Check if there is negative flux
               if (fasm(i, j, g) < 0._dp) negf = 1
            end do
            ys = ys + ydiv(j)
         end do
      end do

      ! Normalize Flux to norm
      if (present(norm) .and. bther /= 1) then
         do g = 1, ng
            do j = 1, ny
               do i = 1, nx
                  fasm(i, j, g) = norm / totp(g) * fasm(i, j, g) * norm
               end do
            end do
         end do
      end if


      ! Print assembly Flux distribution
      write(ounit, *)
      if (negf > 0) write(ounit, *) '    ....WARNING: NEGATIVE FLUX ENCOUNTERED....'
      write(ounit, *)
      if (bther == 1) then
         write(ounit, *) '    Radial Flux Distribution [#neutrons/(cm2.second)]'
         write(ounit, *) '  ===================================================='
      else
         write(ounit, *) '    Radial Flux Distribution'
         write(ounit, *) '  ============================'
      end if

      ip = nx / xm
      ipr = mod(nx, xm) - 1

      !!! Convert to character (zero flux convert to blank spaces)
      do g = 1, ng
         do j = 1, ny
            do i = 1, nx
               if (fasm(i, j, g) > 0.) then
                  write (cflx(i, j, g), '(ES10.3)') fasm(i, j, g)
                  cflx(i, j, g) = trim(adjustl(cflx(i, j, g)))
               else
                  cflx(i, j, g) = '         '
               end if
            end do
         end do
      end do

      ! Print flux
      xs = 1;
      xf = xm
      do k = 1, ip
         write(ounit, '(3X,100I11)') (i, i = xs, xf)
         do j = ny, 1, -1
            write(ounit, '(2X,I4,2X,100A11)') j, (cflx(i, j, 1), i = xs, xf)
            do g = 2, ng
               write(ounit, '(8X,100A11)') (cflx(i, j, g), i = xs, xf)
            end do
            write(ounit, *)
         end do
         write(ounit, *)
         xs = xs + xm
         xf = xf + xm
      end do

      write(ounit, '(3X,100I11)') (i, i = xs, xs + ipr)
      if (xs + ipr > xs) then
         do j = ny, 1, -1
            write(ounit, '(2X,I4,2X,100A11)') j, (cflx(i, j, 1), i = xs, xs + ipr)
            do g = 2, ng
               write(ounit, '(8X,100A11)') (cflx(i, j, g), i = xs, xs + ipr)
            end do
            write(ounit, *)
         end do
      end if
      write(ounit, *)


   end subroutine asmflux

   !****************************************************************************!

   subroutine inp_vtk()

      !
      ! Purpose:
      !    To read VTK card in input.
      !

      implicit none

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '                  >>>>>READING VTK<<<<<'
      write(ounit, *) '           -------------------------------------'
      write(ounit, 1581)

      oname = trim(iname) // '.vtk'
      oname = trim(oname)

      open (unit = vunit, file = oname, status = 'REPLACE', action = 'WRITE')
      call get_vertices()
      call init_vtk()

      1581 format (2x, 'VTK FILE WILL BE GENERATED')

   end subroutine inp_vtk

   !****************************************************************************!

   subroutine flux_atpower(flux)

      !
      ! Purpose:
      !    To calculate flux at given power level
      !


      use sdata, only: ng, nnod, sigf, f0, vdel, pow, ppow

      implicit none

      real(dp), intent(out) :: flux(:, :)
      integer :: g, n
      real(dp) :: tpow

      ! Get total power
      tpow = 0._dp
      do g = 1, ng
         do n = 1, nnod
            tpow = tpow + f0(n, g) * sigf(n, g) * vdel(n)
         end do
      end do

      !Get flux at given power level
      do g = 1, ng
         do n = 1, nnod
            flux(n, g) = f0(n, g) * pow * ppow * 0.01_dp / tpow
         end do
      end do

   end subroutine flux_atpower

   !******************************************************************************!

   subroutine inp_xtab(xbunit)

      !
      ! Purpose:
      !    To read tabular xsec file
      !

      use sdata, only: ng, nmat, nf, m, chi

      implicit none

      integer, intent(in) :: xbunit

      integer :: ln !Line number
      integer :: iost  ! IOSTAT status

      integer, dimension(nf) :: precf
      integer :: i, j, g, h, s, t, u, v

      ! XTAB TYPE
      type :: xfile
         character(len = 100) :: fname             ! XTAB File name
         integer :: cnum                        ! Composition number in the XTAB Files
      end type
      type(xfile), dimension(:), allocatable :: xtab

      logical, dimension(:), allocatable :: noty   !to check if this buffer was read or not?
      integer, parameter :: xunit = 998  !XTAB file unit number
      integer, parameter :: tunit = 999  !XTAB Buffer unit number
      integer :: comm  ! Position of comment mark
      integer :: popt
      integer, dimension(:), allocatable :: group
      integer :: nskip  !Number of lines to skip

      write(ounit, *)
      write(ounit, *)
      write(ounit, *) '           >>>> READING TABULAR XSEC FILE <<<<'
      write(ounit, *) '           --------------------------------------------'

      read(xbunit, *, iostat = iost) ind, ln, ng, nmat  !Read numbef of group and material
      message = ' error in material number'
      call er_message(ounit, iost, ln, message)

      allocate(chi(nmat, ng))
      allocate(xtab(nmat), noty(nmat), m(nmat))
      noty = .true.

      ! Reading input file for XTAB file names and composition number
      do i = 1, nmat
         read(xbunit, '(A2,I5,A200)', iostat = iost) ind, ln, iline
         iline = adjustl(iline)                !Adjust to left
         comm = index(iline, ' ')              ! Get space position
         xtab(i)%fname = iline(1:comm - 1)       !Get xtab file name
         read(iline(comm:100), *) xtab(i)%cnum  ! Get composition number (convert to integer)
         message = ' error in reading XTAB files'
         call er_message(ounit, iost, ln, message)
      end do

      ! Starting to read XTAB File, remove comments and write to buffer file
      do i = 1, nmat
         if (noty(i)) then  !If this composition was not written in buffer
            ! Open XTAB File
            call openfile(xunit, xtab(i)%fname, 'XTAB', 'XTAB File Open Failed--status')

            ! Start removing comments and rewrite into one input XTAB buffer
            call inp_comments(xunit, tunit, '*')

            !This loop to read another composition in the same XTAB file
            do j = i, nmat
               ! If next material has the same file name
               if (trim(adjustl(xtab(i)%fname)) == trim(adjustl(xtab(j)%fname))) then

                  !Read buffer file and saved the xsec and transient data
                  rewind(tunit)
                  read(tunit, *, iostat = iost) ind, ln, m(j)%tadf, m(j)%trod     ! Read input control
                  message = ' ERROR IN XTAB FILE ' // trim(adjustl(xtab(j)%fname)) &
                  // ': CANNOT READ CONTROL PARAMETERS'
                  call er_message(ounit, iost, ln, message, xtab = j)

                  read(tunit, *, iostat = iost) ind, ln, m(j)%nd, m(j)%nb, m(j)%nf, m(j)%nm    ! Read number of branch
                  message = ' ERROR IN XTAB FILE ' // trim(adjustl(xtab(j)%fname))&
                  // ': CANNOT READ BRANCH DIMENSION'
                  call er_message(ounit, iost, ln, message)

                  ! Check branch dimension
                  if (m(j)%nd < 1 .or. m(j)%nb < 1 .or. m(j)%nf < 1 .or. m(j)%nm < 1) then
                     write(ounit, *) ' ERROR: MINIMUM NUMBER OF BRANCH IS 1'
                     write(*, *) ' ERROR: MINIMUM NUMBER OF BRANCH IS 1'
                     stop
                  end if

                  ! Allocate and read branch paramaters
                  call branchpar(tunit, m(j)%nd, j, m(j)%pd, xtab(j)%fname, 'COOLANT DENSITY')  ! Allocate and read coolant dens. branc paramaters
                  call branchpar(tunit, m(j)%nb, j, m(j)%pb, xtab(j)%fname, 'BORON CONCENTRATION')  ! Allocate and read Boron conc. branc paramaters
                  call branchpar(tunit, m(j)%nf, j, m(j)%pf, xtab(j)%fname, 'FUEL TEMPERATURE')  ! Allocate and read fule temp. branc paramaters
                  call branchpar(tunit, m(j)%nm, j, m(j)%pm, xtab(j)%fname, 'MODERATOR TEMPERATURE')  ! Allocate and read moderator temp. branc paramaters

                  ! ALLOCATE XSEC DATA
                  allocate(m(j)%velo(ng))
                  allocate(m(j)%xsec(m(j)%nd, m(j)%nb, m(j)%nf, m(j)%nm))
                  if (m(j)%trod == 1) allocate(m(j)%rxsec(m(j)%nd, m(j)%nb, m(j)%nf, m(j)%nm))
                  do s = 1, m(j)%nd
                     do t = 1, m(j)%nb
                        do u = 1, m(j)%nf
                           do v = 1, m(j)%nm
                              allocate(m(j)%xsec(s, t, u, v)%sigtr(ng))
                              allocate(m(j)%xsec(s, t, u, v)%siga(ng))
                              allocate(m(j)%xsec(s, t, u, v)%sigf(ng))
                              allocate(m(j)%xsec(s, t, u, v)%nuf(ng))
                              allocate(m(j)%xsec(s, t, u, v)%dc(ng, 6))
                              allocate(m(j)%xsec(s, t, u, v)%sigs(ng, ng))
                              if (m(j)%trod == 1) then
                                 allocate(m(j)%rxsec(s, t, u, v)%sigtr(ng))
                                 allocate(m(j)%rxsec(s, t, u, v)%siga(ng))
                                 allocate(m(j)%rxsec(s, t, u, v)%sigf(ng))
                                 allocate(m(j)%rxsec(s, t, u, v)%nuf(ng))
                                 allocate(m(j)%rxsec(s, t, u, v)%dc(ng, 6))
                                 allocate(m(j)%rxsec(s, t, u, v)%sigs(ng, ng))
                              end if
                           end do
                        end do
                     end do
                  end do

                  ! Skip lines to read desired composition in the xtab file
                  nskip = ng * m(j)%nb * m(j)%nf * m(j)%nm
                  if (m(j)%tadf == 1) then  ! IF dc present
                     if (m(j)%trod == 1) then
                        call skipread(tunit, j, (xtab(j)%cnum - 1) * (10 * nskip + 2 * ng * nskip + 4))
                     else
                        call skipread(tunit, j, (xtab(j)%cnum - 1) * (5 * nskip + ng * nskip + 4))
                     end if
                  else if (m(j)%tadf == 2) then
                     if (m(j)%trod == 1) then
                        call skipread(tunit, j, (xtab(j)%cnum - 1) * (20 * nskip + 2 * ng * nskip + 4))
                     else
                        call skipread(tunit, j, (xtab(j)%cnum - 1) * (10 * nskip + ng * nskip + 4))
                     end if
                  else
                     if (m(j)%trod == 1) then
                        call skipread(tunit, j, (xtab(j)%cnum - 1) * (8 * nskip + 2 * ng * nskip + 4))
                     else
                        call skipread(tunit, j, (xtab(j)%cnum - 1) * (4 * nskip + ng * nskip + 4))
                     end if
                  end if

                  ! Read unrodded XSEC
                  call readxs (tunit, j, 0, m(j)%xsec)
                  ! Read rodded XSEC
                  if (m(j)%trod == 1) call readxs (tunit, j, 1, m(j)%rxsec)
                  !Read fission spectrum
                  read(tunit, *, iostat = iost) ind, ln, (chi(j, g), g = 1, ng)
                  message = ' ERROR IN XTAB FILE: CANNOT READ FISSION SPECTRUM'
                  call er_message(ounit, iost, ln, message, xtab = j)
                  !Read neutron Inverse velocity
                  read(tunit, *, iostat = iost) ind, ln, (m(j)%velo(g), g = 1, ng)
                  message = ' ERROR IN XTAB FILE: CANNOT READ Inverse Velocity'
                  call er_message(ounit, iost, ln, message, xtab = j)
                  do g = 1, ng
                     m(j)%velo(g) = 1._dp / m(j)%velo(g)  !COnvert to velocity
                  end do
                  ! Read decay constant
                  read(tunit, *, iostat = iost) ind, ln, (m(j)%lamb(t), t = 1, nf)
                  message = ' ERROR IN XTAB FILE: CANNOT READ DECAY CONSTANT'
                  call er_message(ounit, iost, ln, message, xtab = j)
                  ! Read beta
                  read(tunit, *, iostat = iost) ind, ln, (m(j)%ibeta(t), t = 1, nf)
                  message = ' ERROR IN XTAB FILE: CANNOT READ DELAYED NEUTRON FRACTION'
                  call er_message(ounit, iost, ln, message, xtab = j)

                  ! If read, indicate that it has been read
                  noty(j) = .false.
               end if
            end do

            ! Close XTAB File and buffer file
            close(unit = xunit) ;
            close(unit = tunit)
         end if
      end do

      ! XTAB PRINT OPTION
      read(xbunit, *, iostat = iost) ind, ln, popt
      if (iost == 0 .and. popt > 0) then
         s = 1;
         t = 1;
         u = 1;
         v = 1
         allocate(group(ng))
         do g = 1, ng
            group(g) = g
         end do
         do j = 1, nf
            precf(j) = j
         end do
         do i = 1, nmat
            write(ounit, *)
            write(ounit, 1709) i
            write(ounit, '(A,I3)') '     XTAB FILE ' // trim(adjustl(xtab(i)%fname)) &
            // '. COMPOSITION NUMBER', xtab(i)%cnum
            write(ounit, 1707) 'GROUP', 'TRANSPORT', 'DIFFUSION', 'ABSORPTION', &
            'NU*FISS', 'KAP*FIS', 'FISS. SPCTR', 'NEUTRON VELOCITY'
            do g = 1, ng
               write(ounit, 1706) g, m(i)%xsec(s, t, u, v)%sigtr(g), &
               1. / (3. * m(i)%xsec(s, t, u, v)%sigtr(g)), m(i)%xsec(s, t, u, v)%siga(g), &
               m(i)%xsec(s, t, u, v)%nuf(g), m(i)%xsec(s, t, u, v)%sigf(g), chi(i, g), &
               m(i)%velo(g), m(i)%xsec(s, t, u, v)%dc(g, 1)
            end do
            write(ounit, *) '  --SCATTERING MATRIX--'
            write(ounit, '(4X, A5, 20I11)') "G/G'", (group(g), g = 1, ng)
            do g = 1, ng
               write(ounit, 1015) g, (m(i)%xsec(s, t, u, v)%sigs(g, h), h = 1, ng)
            end do
            write(ounit, *) '  --BETA AND LAMBDA--'
            write(ounit, '(6X, A5, I9, 20I13)') "GROUP", (precf(j), j = 1, nf)
            write(ounit, 1016) 'BETA ', (m(i)%ibeta(j), j = 1, nf)
            write(ounit, 1016) 'LAMBDA ', (m(i)%lamb(j), j = 1, nf)
         end do
      end if

      write(ounit, *)
      write(ounit, *) ' ...XTAB FILE Card is successfully read...'

      1709 format(5x, 'MATERIAL', i3)
      1707 format(2x, a7, a12, a13, a12, a11, a13, a15, a18)
      1706 format(2x, i6, f13.6, 2f12.6, f13.6, es14.5, f12.6, 2es16.5)
      1015 format(4x, i3, f16.6, 20f12.6)
      1016 format(4x, a9, 20es13.5)

      deallocate(xtab, noty)

   end subroutine inp_xtab

   !******************************************************************************!

   subroutine branchpar (tunit, dim, matnum, par, fname, messpar)

      !Purpose: To allocate and read branch paramaters

      implicit none

      integer, intent(in) :: tunit, dim, matnum
      character(len = *), intent(in) :: messpar, fname
      real(dp), dimension(:), allocatable, intent(inout) :: par

      integer :: k
      integer :: ln, iost

      if (dim > 1) then                         ! If branch DIMENSION > 1
         allocate(par(dim))
         read(tunit, *, iostat = iost) ind, ln, par(1:dim)
         message = ' ERROR IN XTAB FILE ' // trim(adjustl(fname))&
         // ': CANNOT READ BRANCH PARAMETERS ' // messpar
         call er_message(ounit, iost, ln, message, xtab = matnum)
         do k = 2, dim
            if (par(k - 1) > par(k)) then
               write(ounit, *) "  ERROR IN XTAB FILE  ", fname
               write(ounit, *) "  ", messpar, " PARAMETER SHALL BE IN ORDER, SMALL to BIG"
               write(*, *) "  ERROR IN XTAB FILE  ", fname
               write(*, *) "  ", messpar, " PARAMETER SHALL BE IN ORDER, SMALL to BIG"
               stop
            end if
         end do
      else
         allocate(par(1))
         par(1) = 0.0    !Arbitrary
      end if

   end subroutine branchpar

   !******************************************************************************!

   subroutine readxs (tunit, mnum, rod, xsec)
      !Purpose: To read xsec in XTAB file

      use sdata, only: m, ng, xbranch
      implicit none

      integer, intent(in) :: tunit, mnum, rod  ! file unit number, material number, and rod indicator
      type(xbranch), dimension(:, :, :, :), intent(inout) :: xsec  !Set INOUT, see: http://www.cs.rpi.edu/~szymansk/OOF90/bugs.html#2
      integer :: iost, ln
      integer :: g, h, s, t, u, v, k

      !Read sigtr
      do g = 1, ng
         do v = 1, m(mnum)%nm
            do u = 1, m(mnum)%nf
               do t = 1, m(mnum)%nb
                  read(tunit, *, iostat = iost) ind, ln, &
                  (xsec(s, t, u, v)%sigtr(g), s = 1, m(mnum)%nd)
                  if (rod == 0) then
                     message = ' ERROR IN XTAB FILE: CANNOT READ TRANSPORT XSEC'
                  else
                     message = ' ERROR IN XTAB FILE: CANNOT READ RODDED TRANSPORT XSEC'
                  end if
                  call er_message(ounit, iost, ln, message, xtab = mnum)
               end do
            end do
         end do
      end do
      !Read siga
      do g = 1, ng
         do v = 1, m(mnum)%nm
            do u = 1, m(mnum)%nf
               do t = 1, m(mnum)%nb
                  read(tunit, *, iostat = iost) ind, ln, &
                  (xsec(s, t, u, v)%siga(g), s = 1, m(mnum)%nd)
                  if (rod == 0) then
                     message = ' ERROR IN XTAB FILE: CANNOT READ ABSORPTION XSEC'
                  else
                     message = ' ERROR IN XTAB FILE: CANNOT READ RODDED ABSORPTION XSEC'
                  end if
                  call er_message(ounit, iost, ln, message, xtab = mnum)
               end do
            end do
         end do
      end do
      !Read nu*sigf
      do g = 1, ng
         do v = 1, m(mnum)%nm
            do u = 1, m(mnum)%nf
               do t = 1, m(mnum)%nb
                  read(tunit, *, iostat = iost) ind, ln, &
                  (xsec(s, t, u, v)%nuf(g), s = 1, m(mnum)%nd)
                  if (rod == 0) then
                     message = ' ERROR IN XTAB FILE: CANNOT READ NU*SIGF XSEC'
                  else
                     message = ' ERROR IN XTAB FILE: CANNOT READ RODDED NU*SIGF XSEC'
                  end if
                  call er_message(ounit, iost, ln, message, xtab = mnum)
               end do
            end do
         end do
      end do
      !Read kappa*sigf
      do g = 1, ng
         do v = 1, m(mnum)%nm
            do u = 1, m(mnum)%nf
               do t = 1, m(mnum)%nb
                  read(tunit, *, iostat = iost) ind, ln, &
                  (xsec(s, t, u, v)%sigf(g), s = 1, m(mnum)%nd)
                  if (rod == 0) then
                     message = ' ERROR IN XTAB FILE: CANNOT READ KAPPA*SIGF XSEC'
                  else
                     message = ' ERROR IN XTAB FILE: CANNOT READ RODDED KAPPA*SIGF XSEC'
                  end if
                  call er_message(ounit, iost, ln, message, xtab = mnum)
               end do
            end do
         end do
      end do
      !Read sigs
      do g = 1, ng
         do h = 1, ng
            do v = 1, m(mnum)%nm
               do u = 1, m(mnum)%nf
                  do t = 1, m(mnum)%nb
                     read(tunit, *, iostat = iost) ind, ln, &
                     (xsec(s, t, u, v)%sigs(g, h), s = 1, m(mnum)%nd)
                     if (rod == 0) then
                        message = ' ERROR IN XTAB FILE: CANNOT READ SCATTERING XSEC'
                     else
                        message = ' ERROR IN XTAB FILE: CANNOT READ RODDED SCATTERING XSEC'
                     end if
                     call er_message(ounit, iost, ln, message, xtab = mnum)
                  end do
               end do
            end do
         end do
      end do
      !Read dc
      if (m(mnum)%tadf == 1) then  ! IF dc present
         do g = 1, ng
            do v = 1, m(mnum)%nm
               do u = 1, m(mnum)%nf
                  do t = 1, m(mnum)%nb
                     read(tunit, *, iostat = iost) ind, ln, &
                     (xsec(s, t, u, v)%dc(g, 1), s = 1, m(mnum)%nd)
                     do k = 1, 6
                        do s = 1, m(mnum)%nd
                           xsec(s, t, u, v)%dc(g, k) = xsec(s, t, u, v)%dc(g, 1)
                        end do
                     end do
                     if (rod == 0) then
                        message = ' ERROR IN XTAB FILE: CANNOT READ ADFs'
                     else
                        message = ' ERROR IN XTAB FILE: CANNOT READ RODDED ADFs'
                     end if
                     call er_message(ounit, iost, ln, message, xtab = mnum)
                  end do
               end do
            end do
         end do
      else if (m(mnum)%tadf == 2) then
         do g = 1, ng
            do k = 1, 6
               do v = 1, m(mnum)%nm
                  do u = 1, m(mnum)%nf
                     do t = 1, m(mnum)%nb
                        read(tunit, *, iostat = iost) ind, ln, &
                        (xsec(s, t, u, v)%dc(g, k), s = 1, m(mnum)%nd)
                        if (rod == 0) then
                           message = ' ERROR IN XTAB FILE: CANNOT READ ADFs'
                        else
                           message = ' ERROR IN XTAB FILE: CANNOT READ RODDED TRANSPORT ADFs'
                        end if
                        call er_message(ounit, iost, ln, message, xtab = mnum)
                     end do
                  end do
               end do
            end do
         end do
      else
         continue
      end if


   end subroutine readxs

   !******************************************************************************!

   subroutine skipread (iunit, matnum, nskip)
      !Purpose: To allocate and read branch paramaters

      implicit none

      integer, intent(in) :: iunit, matnum, nskip
      integer :: i, eof

      do i = 1, nskip
         read (iunit, *, iostat = eof)
         if (eof < 0) then              !Check end of file
            write(ounit, 1131) matnum
            write(ounit, 1132)
            write(*, 1131) matnum
            write(*, 1132)
            stop
         end if
      end do

      1131 format(2x, 'ERROR: END OF FILE REACHED FOR XTAB FILE IN MATERIAL NUMBER ', i3)
      1132 format(2x, 'KOMODO IS STOPPING')

   end subroutine skipread

   !****************************************************************************!

   subroutine get_vertices()

      !
      ! Purpose:
      !    To calculate vertices for vtk output
      !

      use sdata, only: nxx, nyy, nzz, cell, xyz, nver, &
      zdel, xdel, ydel, nnod, ystag, xstag

      implicit none

      integer :: i, j, k, n
      integer, allocatable :: uvw(:, :, :)
      real(dp) :: x, y, z

      nver = nnod * 8

      allocate(cell(nnod))

      n = -1
      z = 0.0
      do k = 1, nzz
         y = 0.0
         do j = 1, nyy
            x = 0.0
            do i = 1, nxx
               if (i >= ystag(j)%smin .and. i <= ystag(j)%smax) then
                  n = n + 1
                  cell(xyz(i, j, k))%vertice(1) % x = x
                  cell(xyz(i, j, k))%vertice(1) % y = y
                  cell(xyz(i, j, k))%vertice(1) % z = z
                  cell(xyz(i, j, k))%vertice(1) % index = n

                  n = n + 1
                  cell(xyz(i, j, k))%vertice(2) % x = x + xdel(i)
                  cell(xyz(i, j, k))%vertice(2) % y = y
                  cell(xyz(i, j, k))%vertice(2) % z = z
                  cell(xyz(i, j, k))%vertice(2) % index = n

                  n = n + 1
                  cell(xyz(i, j, k))%vertice(3) % x = x
                  cell(xyz(i, j, k))%vertice(3) % y = y + ydel(j)
                  cell(xyz(i, j, k))%vertice(3) % z = z
                  cell(xyz(i, j, k))%vertice(3) % index = n

                  n = n + 1
                  cell(xyz(i, j, k))%vertice(4) % x = x + xdel(i)
                  cell(xyz(i, j, k))%vertice(4) % y = y + ydel(j)
                  cell(xyz(i, j, k))%vertice(4) % z = z
                  cell(xyz(i, j, k))%vertice(4) % index = n

                  n = n + 1
                  cell(xyz(i, j, k))%vertice(5) % x = x
                  cell(xyz(i, j, k))%vertice(5) % y = y
                  cell(xyz(i, j, k))%vertice(5) % z = z + zdel(k)
                  cell(xyz(i, j, k))%vertice(5) % index = n

                  n = n + 1
                  cell(xyz(i, j, k))%vertice(6) % x = x + xdel(i)
                  cell(xyz(i, j, k))%vertice(6) % y = y
                  cell(xyz(i, j, k))%vertice(6) % z = z + zdel(k)
                  cell(xyz(i, j, k))%vertice(6) % index = n

                  n = n + 1
                  cell(xyz(i, j, k))%vertice(7) % x = x
                  cell(xyz(i, j, k))%vertice(7) % y = y + ydel(j)
                  cell(xyz(i, j, k))%vertice(7) % z = z + zdel(k)
                  cell(xyz(i, j, k))%vertice(7) % index = n

                  n = n + 1
                  cell(xyz(i, j, k))%vertice(8) % x = x + xdel(i)
                  cell(xyz(i, j, k))%vertice(8) % y = y + ydel(j)
                  cell(xyz(i, j, k))%vertice(8) % z = z + zdel(k)
                  cell(xyz(i, j, k))%vertice(8) % index = n

               end if

               x = x + xdel(i)
            end do
            y = y + ydel(j)
         end do
         z = z + zdel(k)
      end do

   end subroutine

   !******************************************************************************!

   subroutine init_vtk()

      !Purpose: To init vtk output

      use sdata, only: nnod, nver, cell
      implicit none

      integer :: n, i

      write(vunit, '(A26)') '# vtk DataFile Version 4.0'
      write(vunit, '(A33)') 'VTK data file generated by KOMODO'
      write(vunit, '(A5)') 'ASCII'
      write(vunit, '(A25)') 'DATASET UNSTRUCTURED_GRID'

      write(vunit, 3000) 'POINTS', nver, 'double'

      do n = 1, nnod
         do i = 1, 8
            write(vunit, 3001) cell(n) % vertice(i) % x, &
            cell(n) % vertice(i) % y, &
            cell(n) % vertice(i) % z
         end do
      end do

      write(vunit, 3002) 'CELLS', nnod, nnod * 9

      do n = 1, nnod
         write(vunit, 3003) '8', &
         cell(n)%vertice(1)%index, cell(n)%vertice(2)%index, &
         cell(n)%vertice(3)%index, cell(n)%vertice(4)%index, &
         cell(n)%vertice(5)%index, cell(n)%vertice(6)%index, &
         cell(n)%vertice(7)%index, cell(n)%vertice(8)%index
      end do

      write(vunit, 3004) 'CELL_TYPES', nnod

      do i = 1, nnod
         write(vunit, *) '11'
      end do

      write(vunit, 3005) 'CELL_DATA', nnod

      deallocate(cell)

      3000 format(a6, 1x, i8, 1x, a6)
      3001 format(3es13.5)
      3002 format(a5, 1x, i8, 1x, i8)
      3003 format(a1, 8(1x, i8))
      3004 format(a10, 1x, i8)
      3005 format(a9, 1x, i8)

   end subroutine

   !******************************************************************************!

   subroutine print_vtk(time_step)

      !Purpose: To init vtk output

      use sdata, only: nnod, ng, f0, npow, pow, ppow, ftem, mtem, cden, vdel
      implicit none

      integer, intent(in) :: time_step

      integer :: g
      character(5) :: step
      real(dp) :: flux(nnod, ng)
      real(dp) :: total

      step = int_to_char(time_step)

      ! Power
      if (bther == 1) then
         write(vunit, *) 'SCALARS Power_Density_' // trim(adjustl(step)) &
         // '_kW/cm3_ double'
         total = sum(npow)
         npow = npow / total
         npow = npow * pow * ppow * 1.e-5 / vdel
      else
         write(vunit, *) 'SCALARS Relative_Power_' // trim(adjustl(step)) &
         // ' double'
      end if
      write(vunit, *) 'LOOKUP_TABLE default'

      write(vunit, '(ES13.5)') npow

      ! FLux
      if (bther == 1) then
         call flux_atpower(flux)
      else
         flux = f0
      end if
      do g = 1, ng
         if (bther == 1) then
            write(vunit, *) 'SCALARS flux_g_' // int_to_char(g) &
            // '_' // trim(adjustl(step)) // '_n/cm^2*s_ double'
         else
            write(vunit, *) 'SCALARS flux_g_' // int_to_char(g) &
            // '_' // trim(adjustl(step)) // '(Relative) double'
         end if
         write(vunit, *) 'LOOKUP_TABLE default'
         write(vunit, '(ES13.5)') flux(:, g)
      end do

      ! TH PARAMETERS
      if (bther == 1) then
         write(vunit, *) 'SCALARS Fuel_temp_' // trim(adjustl(step)) &
         // '_C_ double'
         write(vunit, *) 'LOOKUP_TABLE default'
         write(vunit, '(ES13.5)') ftem - 273.15

         write(vunit, *) 'SCALARS Moderator_temp_' // trim(adjustl(step)) &
         // '_C_ double'
         write(vunit, *) 'LOOKUP_TABLE default'
         write(vunit, '(ES13.5)') mtem - 273.15

         write(vunit, *) 'SCALARS Coolant_density_' // trim(adjustl(step)) &
         // '_g/cm3_ double'
         write(vunit, *) 'LOOKUP_TABLE default'
         write(vunit, '(ES13.5)') cden
      end if

   end subroutine

   !******************************************************************************!

   function int_to_char(num, digit) result(char_out)
      implicit none
      character(len = 32) :: char
      character(len = 32) :: form
      character(len = :), allocatable :: char_out
      integer :: n_out
      integer :: num
      integer :: i
      integer, optional :: digit
      integer :: digit0
      integer :: digit1

      form = "(I30)"
      write(char, form) num
      char = trim(adjustl(char))

      if (present(digit)) then
         digit0 = abs(digit)
         digit1 = get_digit(num)
         do i = 1, digit0 - digit1
            if (digit > 0) then
               char = "0" // trim(char)
            else
               char = " " // trim(char)
            end if
         end do
      end if

      n_out = len(trim(char))
      allocate(character(n_out) :: char_out)
      char_out = char(1:n_out)

   end function

   !******************************************************************************!
   function get_digit(num) result(digit0)
      implicit none
      integer :: digit0
      integer :: num
      integer :: t

      t = num
      digit0 = 1
      do
         if (t < 10) then
            exit
         end if
         t = t / 10
         digit0 = digit0 + 1
      end do

   end function

end module
