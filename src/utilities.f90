module utilities

    use iso_fortran_env, only: real64, int64, error_unit

    implicit none
  
    save

    public

    integer, private, parameter :: dp = real64
    integer, private, parameter :: i8 = int64

    interface n2c
    module procedure real_to_char
    module procedure int_to_char
    end interface n2c

    interface c2n
    module procedure char_to_real
    end interface c2n

    contains

    !================================================================!
    ! To solve Tridiagonal matrix
    !================================================================!

    subroutine TridiaSolve(a,b,c,d,x,error_flag)

        ! a vector is upper diagnal, b vector is diagonal, and c vector is lower diagonal
        ! d vector is the RHS

        real(dp), dimension(:), intent(inout) :: a, b, c, d
        real(dp),               intent(out)   :: x(:)
        integer,                intent(out)   :: error_flag
        
        integer :: i, n
        real(dp) :: pivot
        
        n = size(d)
        error_flag = 0
        
        ! Gauss Elimination
        c(1) = c(1)/b(1)
        d(1) = d(1)/b(1)
        do i = 2, n
            pivot = (b(i) - a(i) * c(i-1))
            c(i) = c(i) / pivot
            if (abs(pivot) < 1.e-20) then
                error_flag = 1
                return
            end if
            d(i) = (d(i) - a(i) * d(i-1)) / pivot
        end do
        
        ! Back Substitution
        x(n) = d(n)
        do i = n-1, 1, -1
            x(i) = d(i) - c(i) * x(i+1)
        end do
        
    end subroutine

    !===============================================================================================!
    ! To solve LU
    !===============================================================================================!

    function mat_solve(msize,mat,b) result(x)
  
        implicit none
  
        integer, intent(in)                  :: msize  ! node and and matrix size
        real(dp), dimension(:,:), intent(in) :: mat           ! the matrix A
        real(dp), dimension(:), intent(in)   :: b            ! the vector b
        real(dp), dimension(msize)           :: x            ! the vector b
  
        real(dp), dimension(msize,msize) :: L, U
        real(dp), dimension(msize) :: y
        real(dp) :: piv, isum
        integer :: i, j, k
  
        U = mat
        L = 0._dp
  
        ! Start matrix decomposition
        do i= 1, msize
            if (ABS(mat(i,i)) < 10e-5) stop "fail"
            L(i,i) = 1._DP
            do j= i+1, msize
                piv = U(j,i)/U(i,i)
                L(j,i) = piv
                do k= i, msize
                    U(j,k) = U(j,k) - piv*U(i,k)
                end do
                U(j,i) = 0._dp
            end do
        end do  
  
        !Solve y in Ly = b (Forward substitution)
        y(1) = b(1)
        do i=2,msize
            isum = 0._dp
            do k =1, i-1
                isum = isum + L(i,k)*y(k)
            end do
            y(i) = b(i)-isum
        end do
  
        ! Solve x in Ux=y(Backward substitution)
        x(msize) = y(msize)/U(msize,msize)
        do i = msize-1,1,-1
            isum = 0._dp
            do k =i+1,msize
                isum = isum + U(i,k)*x(k)
            end do
            x(i) = (y(i)-isum) / U(i,i)
        end do
  
    end function

    !===============================================================================================!
    ! To factorize A to L, U. L and U matrices combined into LU
    !===============================================================================================!

    pure subroutine LU_decompose(A, LU, error_flag)
    
        real(dp), intent(in)  :: A(:,:)        ! the matrix A
        real(dp), intent(out) :: LU(:,:)
        integer, intent(out)  :: error_flag

        real(dp) :: pivot
        integer  :: msize
        integer  :: i, j, k
    
        error_flag = 0
        msize = size(A, dim=1)
        LU = A
    
        ! Start matrix decomposition
        do i= 1, msize-1
            if (ABS(LU(i,i)) < 1.0e-5) then
                error_flag = -1
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
            write(error_unit,*) '  File name must be present in write or read mode'
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
            write(error_unit,*) 'mode ' // rw // ' is not known'
            stop '  Stop in open_file function'
        end if

        if (iost .ne. 0) then
            write(error_unit,*) '  Cannot find file : ' // file_name
            stop '  Stop in open_file function'
        end if
        
    end function

    !===============================================================================================!
    ! To open file                                                                                  !
    !===============================================================================================!

    subroutine write_vtk(file_name, x_min, x_max, y_min, y_max, &
        z_min, z_max, dat)

        character(len=*), optional, intent(in)  :: file_name
        real(dp), intent(in) :: x_min, x_max, y_min, y_max, z_min, z_max
        real(dp)             :: dat(:,:,:)
        
        integer               :: fid
        integer               :: nx, ny, nz
        real(dp), allocatable :: x_coord(:)
        real(dp), allocatable :: y_coord(:)
        real(dp), allocatable :: z_coord(:)
        real(dp)              :: del
        integer               :: i, j, k

        nx = size(dat, dim=1)
        ny = size(dat, dim=2)
        nz = size(dat, dim=3)

        allocate(x_coord(nx+1))
        x_coord(1) = x_min
        del = (x_max-x_min) / real(nx)
        do i = 2, nx+1
            x_coord(i) = x_coord(i-1) + del
        end do

        allocate(y_coord(ny+1))
        y_coord(1) = y_min
        del = (y_max-y_min) / real(ny)
        do i = 2, ny+1
            y_coord(i) = y_coord(i-1) + del
        end do

        allocate(z_coord(nz+1))
        z_coord(1) = z_min
        del = (z_max-z_min) / real(nz)
        do i = 2, nz+1
            z_coord(i) = z_coord(i-1) + del
        end do

        fid = open_file('write', file_name)

        write(fid, '(A)') "# vtk DataFile Version 2.0"
        write(fid, '(A)') "VTK data file generated Imron"
        write(fid, '(A)') "ASCII"
        write(fid, '(A)') "DATASET RECTILINEAR_GRID"
        write(fid, '(A)') "DIMENSIONS " // n2c(nx+1) // " " // n2c(ny+1) // " " // n2c(nz+1)
        write(fid, '(A)') "X_COORDINATES "// n2c(nx+1) // " float"
        do i = 1, nx+1
            write(fid, '(ES14.5)') x_coord(i)
        end do
        write(fid, '(A)') "Y_COORDINATES "// n2c(ny+1) // " float"
        do i = 1, ny+1
            write(fid, '(ES14.5)') y_coord(i)
        end do
        write(fid, '(A)') "Z_COORDINATES "// n2c(nz+1) // " float"
        do i = 1, nz+1
            write(fid, '(ES14.5)') x_coord(i)
        end do
        write(fid, '(A)') "CELL_DATA "// n2c(nx*ny*nz)
        write(fid, '(A)') "SCALARS average float"
        write(fid, '(A)') "LOOKUP_TABLE default"
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    write(fid,'(ES13.5)') dat(i,j,k)
                end do
            end do
        end do
        
    end subroutine

    !===============================================================================================!
    ! to convert character to lower case                                                            !
    !===============================================================================================!

    pure function lower_case(iword) result(word)

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
     
    pure function upper_case(iword) result(word)
        
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

        write(error_unit,*)
        write(error_unit,*)''//achar(27)//'[31m OOPS! WE FOUND AN ERROR.'//achar(27)//'[0m'
        write(error_unit,*)'  ' // message

        write(fid,*)
        write(fid,*)'OOPS! WE FOUND AN ERROR.'
        write(fid,*)'  ' // message

        stop
        
    end subroutine

    !==============================================================================!
    ! convert integer to character
    !==============================================================================!

    pure function int_to_char(num, digit) result(char_out)

        integer, intent(in)             :: num
        integer, intent(in),  optional  :: digit
        character(len=:), allocatable   :: char_out
        
        character(len=20)             :: char
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

    !==============================================================================!
    ! get number of digit
    !==============================================================================!

    pure function get_digit(num) result(digit0)
        
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

    pure function real_to_char(num, digit_int) result(char_out)
        
        real(dp), intent(in)   :: num
        integer , intent(in)   :: digit_int

        character(len=30)    :: char
        character(len=30)    :: form
        character(len=30)    :: digit
        character(len=:), allocatable :: char_out
        integer                       :: n_out
     
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

    !==============================================================================!
    ! convert character to real
    !==============================================================================!

    function char_to_real(char) result(number)
        
        character(len=*), intent(in) :: char
        real(dp)                     :: number

        character(len(char)) :: char2
        integer    :: ios
     
        char2 = trim(adjustl(char))
        
        read(char2,*,iostat=ios) number
        if (ios .ne. 0) then
            write(error_unit,*) "failed to convert character " // char2 // " to real"
            stop
        end if
     
    end function

    !==============================================================================!
    ! to skip reading some lines
    !==============================================================================!

    subroutine skip_lines(fid, n_line)
        
        integer, intent(in) :: fid, n_line

        integer :: i, eof

        do i = 1, n_line
            read(fid,*, iostat=eof)
            if (eof < 0) exit
        end do
     
    end subroutine

    !==============================================================================!
    ! divide a 2D matrix into several smaller blocks
    !==============================================================================!

    function matrix_block(matrix, px, py) result(res)
        
        real(dp), intent(in)  :: matrix(:,:)
        integer, intent(in)   :: px, py
        real(dp), allocatable :: res(:,:,:,:)

        integer:: nxx, nyy, ny, nx
        integer :: i, j, ii, jj
        integer :: is, js

        nxx = size(matrix, dim=1)
        nyy = size(matrix, dim=2)
     
        if (mod(nxx, px) .ne. 0) then
            write(error_unit,*) "matrix_block: matrix dim 1 cannot be divided by nx"
            stop
        end if
        if (mod(nyy, py) .ne. 0) then
            write(error_unit,*) "matrix_block: matrix dim 2 cannot be divided by ny"
            stop
        end if

        nx = nxx / px
        ny = nyy / py
        allocate(res(nx, ny, px, py))

        js = 0
        do j = 1, ny
            do jj = 1, py
                js = js + 1
                is = 0
                do i = 1, nx
                    do ii = 1, px
                        is = is + 1
                        res(i,j,ii,jj) = matrix(is,js)
                    end do
                end do
            end do
        end do
     
    end function

    !==============================================================================!
    ! flip a 2D matrix in x axis
    !==============================================================================!

    function flipy(matrix) result(res)
        
        real(dp), intent(in)      :: matrix(:,:)
        real(dp)                  :: res(size(matrix, dim=1), size(matrix, dim=2))

        integer :: n, i

        n = size(matrix, dim=1)
        do i = 1, size(matrix, dim=1)
            res(i,:) = matrix(n, :)
            n = n - 1
        end do

    end function

    !==============================================================================!
    ! flip a 2D matrix in y axis
    !==============================================================================!

    function flipx(matrix) result(res)
        
        real(dp), intent(in)      :: matrix(:,:)
        real(dp)                  :: res(size(matrix, dim=1), size(matrix, dim=2))

        integer :: n, j

        n = size(matrix, dim=2)
        do j = 1, size(matrix, dim=2)
            res(:,j) = matrix(:, n)
            n = n - 1
        end do

    end function

    !==============================================================================!
    ! rotate clock-wise for 90 degree
    !==============================================================================!

    function rotate90(matrix) result(res)
        
        real(dp), intent(in)   :: matrix(:,:)
        real(dp)               :: res(size(matrix, dim=2),size(matrix, dim=1))
        
        integer :: nx
        integer :: i

        nx = size(matrix, dim=1)

        do i = 1, nx
            res(:,i) = matrix(i,:)
        end do
        res = flipy(res)

    end function

    !==============================================================================!
    ! rotate clock-wise for 180 degree
    !==============================================================================!

    function rotate180(matrix) result(res)
        
        real(dp), intent(in)   :: matrix(:,:)
        real(dp)               :: res(size(matrix, dim=1),size(matrix, dim=2))
        
        real(dp)    :: tmp(size(matrix, dim=2),size(matrix, dim=1))

        tmp = rotate90(matrix)
        res = rotate90(tmp)

    end function

    !==============================================================================!
    ! rotate clock-wise for 270 degree
    !==============================================================================!

    function rotate270(matrix) result(res)
        
        real(dp), intent(in)   :: matrix(:,:)
        real(dp)               :: res(size(matrix, dim=2),size(matrix, dim=1))
        
        real(dp)  :: tmp(size(matrix, dim=1),size(matrix, dim=2))

        tmp = rotate180(matrix)
        res = rotate90(tmp)

    end function
    
    !==============================================================================!
    ! find average of an 2d array (zero excluded)
    !==============================================================================!

    function mean(matrix) result(res)
        
        real(dp), intent(in)   :: matrix(:,:)
        real(dp)               :: res
        
        integer :: i, j, n
        
        n = 0
        do j = 1, size(matrix, dim=2)
            do i = 1, size(matrix, dim=1)
                if (matrix(i,j) > 0.0) n = n + 1
            end do
        end do
        
        res = sum(matrix) / real(n)

    end function

end module