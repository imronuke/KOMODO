module util

    implicit none
  
    save

    contains

    !===============================================================================================!
    ! To open file                                                                                  !
    !===============================================================================================!

    function open_file(rw, file_name) result(file_unit)

        character(len=*), intent(in)            :: rw
        character(len=*), optional, intent(in)  :: file_name
        integer                                 :: file_unit

        integer  :: iost

        if (.not. present(file_name) .and. rw .ne. "scratch") then
            print *, '  File name must be present in write or read mode'
            stop '  Stop in open_file function'
        end if
    
        if (rw == 'read') then
            open (newunit=file_unit, file=file_name, status='old', action='read', &
            iostat = iost)
        elseif (rw == 'scratch') then
            open (newunit=file_unit, status='scratch', action='readwrite', &
            iostat = iost)
        elseif (rw == 'write') then
            open (newunit=file_unit, file=file_name, status='replace', action='write', &
            iostat = iost)
        else
            print *, 'mode ' // rw // ' is not known'
            stop '  Stop in open_file function'
        end if

        if (iost .ne. 0) then
            print *, '  Cannot find file : ' // file_name
            stop '  Stop in open_file function'
        end if
        
    end function

    !===============================================================================================!
    ! to convert character to lower case                                                            !
    !===============================================================================================!

    function lower_case(iword) result(word)
        ! convert a word to lower CASE
        character (len=*), intent(in)   :: iword
        character (len=200)             ::  word
        integer                         :: i,ic,nlen
     
        word = iword
        nlen = len(word)
        do i=1,nlen
           ic = ichar(word(i:i))
           if (ic >= 65 .and. ic < 91) word(i:i) = char(ic+32)
        end do

    end function

    !===============================================================================================!
    ! to convert character to upper case                                                            !
    !===============================================================================================!
     
    function upper_case(iword) result(word)
        
        character (len=*), intent(in)   :: iword
        character (len=200)             ::  word
        integer                         :: i,ic,nlen
     
        word = iword
        nlen = len(word)
        do i=1,nlen
           ic = ichar(word(i:i))
           if (ic >= 97 .and. ic < 123) word(i:i) = char(ic-32)
        end do
        
    end function

end module