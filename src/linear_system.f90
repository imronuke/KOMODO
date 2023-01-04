module linear_system

    use iso_fortran_env, only: real64

    implicit none

    private

    save

    integer, parameter :: dp = real64
    
    ! Linear system Matrix (Stored in Compressed Sparse Row aka CSR)
    type :: matrix_elements
        real(dp), allocatable     :: elmn(:)               ! Non-zero elements of FDM matrix for a row
    end type

    ! Matrix indexing
    type :: matrix_index
        integer, allocatable :: row(:)      ! Row pointer
        integer, allocatable :: col(:)      ! Column index for the non-zero element of the FDM matrix
    end type
    ! type(matrix_index) :: ind                          ! Index of the FDM matrix

    ! The matrix and index
    type, public :: matrix
        type(matrix_elements), public, allocatable :: mtx(:)
        type(matrix_index), public                 :: ind
    end type

    real(dp), allocatable  :: r(:), rs(:), v(:), p(:), s(:), t(:), tmp(:)

    public :: allocate_matrix, linear_solver

    contains

    !===============================================================================================!
    ! Allocate matrix data                                                                          !
    !===============================================================================================!

    subroutine allocate_matrix(m, ng, nnod, n_non_zero)

        type(matrix), intent(inout) :: m
        integer, intent(in)         :: ng, nnod, n_non_zero
        logical                     :: first = .true.

        integer  :: g

        allocate(m % mtx(ng))
        do g = 1, ng
            allocate(m % mtx(g) % elmn(n_non_zero))
        end do

        allocate(m % ind % row(nnod+1))
        allocate(m % ind % col(n_non_zero))

        if (first) then
            allocate(r  (nnod))
            allocate(rs (nnod))
            allocate(v  (nnod))
            allocate(p  (nnod))
            allocate(t  (nnod))
            allocate(s  (nnod))
            allocate(tmp(nnod))
            first = .false.
        end if

    end subroutine

    !===============================================================================================!
    !solve linear of system equation with BiCGSTAB method                                           !
    ! adapted from : https://www.cfd-online.com/Wiki/Sample_code_for_BiCGSTAB_-_Fortran_90          !
    !===============================================================================================!

    subroutine linear_solver(mtx, ind, max_iter, b, x)

        type(matrix_elements), intent(in)    :: mtx
        type(matrix_index), intent(in)       :: ind
        integer, intent(in)      :: max_iter     ! Max. number of iteration and group number
        real(dp), intent(in)     :: b(:)         ! RHS
        real(dp), intent(inout)  :: x(:)
    
        real(dp) :: rho, rho_prev
        real(dp) :: alpha, omega, beta, theta
        integer  :: i, n

        n = size(ind % row) - 1
    
        call sp_matvec(n, mtx, ind, x, r)
    
        call xpby(n, b, -1._dp, r, rs)
        call xew(n, rs, r)
        
        rho   = 1.0_dp
        alpha = 1.0_dp
        omega = 1.0_dp
        v     = 0.0_dp
        p     = 0.0_dp
    
        do i = 1, max_iter
            rho_prev = rho
            rho      = dproduct(n, rs,r)
            beta     = (rho / rho_prev) * (alpha / omega)
            call xpby(n, p, -omega, v, tmp)
            call xpby(n, r, beta, tmp, p)
            call sp_matvec(n, mtx, ind, p, v)
            alpha    = rho/dproduct(n, rs, v)
            call xpby(n, r, -alpha, v, s)
            call sp_matvec(n, mtx, ind, s, t)
            theta    = dproduct(n, t, t)
            omega    = dproduct(n, t, s) / theta
            call axpby(n, alpha, p, omega, s, tmp)
            call xpby(n, x, 1.0_dp, tmp, r)
            call xew(n, r, x)
            call xpby(n, s, -omega, t, r)
        end do

    end subroutine

    !===============================================================================================!
    !to perform sparse matrix vector multiplication Axb. A is a sparse square matrixin CSR format   !
    !===============================================================================================!

    pure subroutine sp_matvec(n, mtx, ind, x, rx)

        integer, intent(in)                  :: n        ! length of diagonal
        type(matrix_elements), intent(in)    :: mtx
        type(matrix_index), intent(in)       :: ind
        real(dp), dimension(:), intent(in)   :: x       ! vector
        real(dp), intent(out)                :: rx(:)   ! resulting vector
     
        integer   :: i, j
        integer   :: row_start, row_end
        real(dp)  :: tmpsum
     

        do concurrent (i = 1:n)
            tmpsum = 0.
            row_start = ind % row(i)
            row_end   = ind % row(i+1) - 1
            do j = row_start, row_end
                tmpsum = tmpsum + mtx % elmn(j) * x(ind % col(j)) 
            end do
            rx(i) = tmpsum
        end do
     
    end subroutine

    !===============================================================================================!
    !to calculate dot product of vectors a and b          !
    !===============================================================================================!

    pure function dproduct(n, a, b) result(rx)

        integer, intent(in)                  :: n
        real(dp), dimension(:), intent(in)   :: a, b   ! vector
        real(dp)                             :: rx     ! resulting vector

        real(dp)  :: x
        integer   :: i
        
        x = 0._dp
        do concurrent (i = 1:n)
            x = x + a(i) * b(i)
        end do
     
        rx = x
     
    end function

    !===============================================================================================!
    !to perform a*\vec{x}+b*\vec{y} calculation                                                     !
    !===============================================================================================!

    pure subroutine axpby(n, alpha, x, beta, y, w)

        integer, intent(in)   :: n
        real(dp), intent(in)  :: alpha, beta, x(:), y(:)
        real(dp), intent(out) :: w(:)

        integer               :: i
      
        do concurrent (i = 1:n)
            w(i) = alpha*x(i) + beta*y(i)
        enddo
     
    end subroutine

    !===============================================================================================!
    !to perform \vec{x}+b*\vec{y} calculation                                                     !
    !===============================================================================================!

    pure subroutine xpby(n, x, beta, y, w)

        integer, intent(in)   :: n
        real(dp), intent(in)  :: beta, x(:), y(:)
        real(dp), intent(out) :: w(:)

        integer               :: i
    
        do concurrent (i = 1:n)
            w(i) = x(i) + beta*y(i)
        enddo
     
    end subroutine

    !===============================================================================================!
    !to perform copy value x to y (x=y)                                                   !
    !===============================================================================================!

    pure subroutine xew(n, x, w)

        integer, intent(in)    :: n
        real(dp), intent(in)   :: x(:)
        real(dp), intent(out)  :: w(n)
        integer                :: i
    
        do concurrent (i = 1:n)
            w(i) = x(i)
        enddo
     
    end subroutine
    
end module
    