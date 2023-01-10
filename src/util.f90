module util

    use iso_fortran_env, only: real64, int64

    implicit none
  
    save

    public

    integer, private, parameter :: dp = real64
    integer, private, parameter :: i8 = int64

    interface n2c
    module procedure real_to_char
    module procedure int_to_char
    end interface n2c

    private :: LU_decompose, LU_subst

    contains

    !================================================================!
    ! To solve Tridiagonal matrix
    !================================================================!

    subroutine TridiaSolve(a,b,c,d,x)

        ! a vector is upper diagnal, b vector is diagonal, and c vector is lower diagonal
        ! d vector is the RHS

        real(dp), dimension(:), intent(inout) :: a, b, c, d
        real(dp),               intent(out)   :: x(:)
        
        integer :: i, n
        real(dp) :: pivot
        
        n = size(d)
        
        ! Gauss Elimination
        c(1) = c(1)/b(1)
        d(1) = d(1)/b(1)
        do i = 2, n
            pivot = (b(i) - a(i) * c(i-1))
            c(i) = c(i) / pivot
            d(i) = (d(i) - a(i) * d(i-1)) / pivot
        end do
        
        ! Back Substitution
        x(n) = d(n)
        do i = n-1, 1, -1
            x(i) = d(i) - c(i) * x(i+1)
        end do
        
    end subroutine

    !===============================================================================================!
    ! To factorize A to L, U. L and U matrices combined into LU
    !===============================================================================================!

    function LU_solve(node, nsize, A, b) result(x)
    
        integer, intent(in)   :: node, nsize
        real(dp), intent(in)  :: A(:,:)   ! LHS matrix
        real(dp), intent(in)  :: b(:)     ! RHS vector
        real(dp)              :: x(nsize)

        integer  :: err
        real(dp) :: LU(nsize, nsize)

        call LU_decompose(A, LU, err)
        if (err < 0) then
            print *, "LU solve failed for node " // n2c(node)
            stop
        end if

        x = LU_subst(LU, b)

    end function

    !===============================================================================================!
    ! To factorize A to L, U. L and U matrices combined into LU
    !===============================================================================================!

    pure subroutine LU_decompose(A, LU, err)
    
        real(dp), intent(in)  :: A(:,:)        ! the matrix A
        real(dp), intent(out) :: LU(:,:)
        integer, intent(out)  :: err

        real(dp) :: pivot
        integer  :: msize
        integer  :: i, j, k
    
        err = 0
        msize = size(A, dim=1)
        LU = A
    
        ! Start matrix decomposition
        do i= 1, msize-1
            if (ABS(LU(i,i)) < 1.0e-5) then
                err = -1
                return
            end if
            do j= i+1, msize
                pivot = LU(j,i)/LU(i,i)
                LU(j,i) = pivot
                do k= i+1, msize
                    LU(j,k) = LU(j,k) - pivot*LU(i,k)
                end do
            end do
        end do

    end subroutine

    !===============================================================================================!
    ! solve decomposed LU matrix by forward-bacward substitution
    !===============================================================================================!

    pure function LU_subst(LU, b) result(x)
  
        real(dp), intent(in)  :: LU(:,:)
        real(dp), intent(in)  :: b(:)            ! the vector b
        real(dp), allocatable :: x(:)            ! the vector x
  

        real(dp), allocatable :: y(:)
        integer               :: msize
        integer               :: i, j
        real(dp)              :: isum
  
        msize = size(b)        
        allocate(x(msize))
        allocate(y(msize))
  
        !Solve y in Ly = b (Forward substitution)
        y(1) = b(1)
        do i=2,msize
            isum = 0._dp
            do j =1, i-1
                isum = isum + LU(i,j)*y(j)
            end do
            y(i) = b(i)-isum
        end do
  
        ! Solve x in Ux=y(Backward substitution)
        x(msize) = y(msize)/LU(msize,msize)
        do i = msize-1,1,-1
            isum = 0._dp
            do j =i+1,msize
                isum = isum + LU(i,j)*x(j)
            end do
            x(i) = (y(i)-isum) / LU(i,i)
        end do
  
    end function

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

    !===============================================================================================!
    ! terminate program and give error message                                                            !
    !===============================================================================================!
     
    subroutine fatal_error(fid, message)
        
        integer, intent(in)             :: fid
        character (len=*), intent(in)   :: message

        write(*,*)
        write(*,*)''//achar(27)//'[31m OOPS! WE FOUND AN ERROR.'//achar(27)//'[0m'
        write(*,*)'  ' // message

        write(fid,*)
        write(fid,*)'OOPS! WE FOUND AN ERROR.'
        write(fid,*)'  ' // message

        stop
        
    end subroutine

    !==============================================================================!
    ! convert integer to character
    !==============================================================================!

    function int_to_char(num, digit) result(char_out)

        integer                       :: num
        integer, optional             :: digit
        character(len=20)             :: char

        character(len=:), allocatable :: char_out
        integer                       :: n_out
        integer                       :: i
        integer                       :: digit0
        integer                       :: digit1
     
        write(char,*) num
        char=trim(adjustl(char))
        
        if(present(digit)) then
           digit0 = abs(digit)
           digit1 = get_digit(num)
           do i = 1, digit0 - digit1
              if(digit > 0) then
                 char = "0"//trim(char)
              else
                 char = " "//trim(char)
              end if
           end do
        end if
        
        n_out = len(trim(char))
        allocate(character(n_out) :: char_out)
        char_out = char(1:n_out)
     
    end function

    function get_digit(num) result(digit0)
        
        integer, intent(in)  :: num
        integer  :: digit0

        integer  :: t
        
        t = num
        digit0 = 1
        do
           if(t < 10) then
                exit
           end if
           t = t/10
           digit0 = digit0 + 1
        end do
        
    end function get_digit

    !==============================================================================!
    ! convert real to character
    !==============================================================================!

    function real_to_char(num, digit_int) result(char_out)
        implicit none
        character(len=30)    :: char
        character(len=30)    :: form
        character(len=30)    :: digit
        character(len=:), allocatable :: char_out
        integer                       :: n_out
        real(dp)                      :: num
        integer                       :: digit_int
     
        write(digit, *) abs(digit_int)
     
        form = ""
     
        if(digit_int >= 0) then
           form="(f30."//trim(adjustl(digit))//")"
        else
           form="(Es30."//trim(adjustl(digit))//")"
        end if
     
        write(char,form) num
     
        char = trim(adjustl(char))
        
        n_out = len(trim(char))
        allocate(character(n_out) :: char_out)
        char_out = char(1:n_out)
     
    end function real_to_char

end module