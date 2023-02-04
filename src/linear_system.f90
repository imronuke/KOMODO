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
        integer, allocatable :: diag(:)     ! row pointer of the diagonal element 
    end type
    ! type(matrix_index) :: ind                          ! Index of the FDM matrix

    ! The matrix and index
    type, public :: matrix
        type(matrix_elements), public, allocatable :: mtx(:)
        type(matrix_index), public                 :: ind
    end type

    real(dp), allocatable  :: r(:), r0(:), v(:), p(:), s(:), t(:), z(:)

    public :: allocate_matrix, linear_solver!, p_linear_solver

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
        allocate(m % ind % diag(nnod))
        allocate(m % ind % col(n_non_zero))

        if (first) then
            allocate(r  (nnod))
            allocate(r0 (nnod))
            allocate(v  (nnod))
            allocate(p  (nnod))
            allocate(t  (nnod))
            allocate(s  (nnod))
            allocate(z  (nnod))
            first = .false.
        end if

    end subroutine

    !===============================================================================================!
    ! solve linear of system equation with BiCGSTAB method                                           !
    ! adapted from : https://www.cfd-online.com/Wiki/Sample_code_for_BiCGSTAB_-_Fortran_90          !
    !===============================================================================================!

    subroutine linear_solver(mtx, ind, max_iter, b, x)

        type(matrix_elements), intent(inout) :: mtx
        type(matrix_index), intent(in)       :: ind

        integer, intent(in)      :: max_iter     ! Max. number of iteration and group number
        real(dp), intent(in)     :: b(:)         ! RHS
        real(dp), intent(inout)  :: x(:)
    
        real(dp) :: rho, rho_prev
        real(dp) :: alpha, omega, beta, theta
        integer  :: i, n

        n = size(ind % row) - 1

        z   = sp_matvec(n, mtx, ind, x)
        r   = xpby(n, b, -1._dp, z)
        r0  = copy(n, r)
        
        rho   = 1.0_dp
        alpha = 1.0_dp
        omega = 1.0_dp
        v     = 0.0_dp
        p     = 0.0_dp

        do i = 1, max_iter
            rho_prev = rho
            rho      = dproduct(n, r0,r)
            beta     = (rho / rho_prev) * (alpha / omega)
            v        = xpby(n, p, -omega, v)
            p        = xpby(n, r, beta, v)
            v        = sp_matvec(n, mtx, ind, p)
            alpha    = rho/dproduct(n, r0, v)
            s        = xpby(n, r, -alpha, v)
            t        = sp_matvec(n, mtx, ind, s)
            theta    = dproduct(n, t, t)
            omega    = dproduct(n, t, s) / theta
            r        = axpby(n, alpha, p, omega, s)
            x        = xpby(n, x, 1.0_dp, r)
            r        = xpby(n, s, -omega, t)
        end do

    end subroutine

    ! subroutine p_linear_solver(mtx, ind, max_iter, b, x)

    !     type(matrix_elements), intent(in)    :: mtx
    !     type(matrix_index), intent(in)       :: ind
    !     integer, intent(in)      :: max_iter     ! Max. number of iteration and group number
    !     real(dp), intent(in)     :: b(:)         ! RHS
    !     real(dp), intent(inout)  :: x(:)
    
    !     real(dp) :: rho, rho_prev
    !     real(dp) :: alpha, omega, beta, theta
    !     integer  :: i, n
    !     real(dp) :: norm_r

    !     n = size(ind % row) - 1

    !     ! v   = sgs_precond_multiply(n, mtx, ind, x)
    !     z   = sgs_precond_solver(n, mtx, ind, x)
    !     z   = sp_matvec(n, mtx, ind, z)
    !     r   = xpby(n, b, -1._dp, z)
    !     r0  = equal(n, r)
        
    !     rho   = 1.0_dp
    !     alpha = 1.0_dp
    !     omega = 1.0_dp
    !     v     = 0.0_dp
    !     p     = 0.0_dp

    !     norm_r = norm2(r)

    !     i = 0
    !     do while(norm_r > 1.e-8)
    !         i = i + 1
    !         rho_prev = rho
    !         rho      = dproduct(n, r0,r)
    !         beta     = (rho / rho_prev) * (alpha / omega)
    !         v        = xpby(n, p, -omega, v)
    !         p        = xpby(n, r, beta, v)
    !         z        = sgs_precond_solver(n, mtx, ind, p)
    !         v        = sp_matvec(n, mtx, ind, z)
    !         alpha    = rho/dproduct(n, r0, v)
    !         s        = xpby(n, r, -alpha, v)
    !         z        = sgs_precond_solver(n, mtx, ind, s)
    !         t        = sp_matvec(n, mtx, ind, z)
    !         theta    = dproduct(n, t, t)
    !         omega    = dproduct(n, t, s) / theta
    !         r        = axpby(n, alpha, p, omega, s)
    !         x        = xpby(n, x, 1.0_dp, r)
    !         r        = xpby(n, s, -omega, t)
    !         norm_r   = norm2(r)
    !     end do

    !     x = sgs_precond_solver(n, mtx, ind, x)
    !     write(*,*) i, sum(x)
    !     stop

    !     ! do i = 1, max_iter
    !     !     rho_prev = rho
    !     !     rho      = dproduct(n, r0,r)
    !     !     beta     = (rho / rho_prev) * (alpha / omega)
    !     !     v        = xpby(n, p, -omega, v)
    !     !     p        = xpby(n, r, beta, v)
    !     !     z        = sgs_precond_solver(n, mtx, ind, p)
    !     !     v        = sp_matvec(n, mtx, ind, z)
    !     !     alpha    = rho/dproduct(n, r0, v)
    !     !     s        = xpby(n, r, -alpha, v)
    !     !     z        = sgs_precond_solver(n, mtx, ind, s)
    !     !     t        = sp_matvec(n, mtx, ind, z)
    !     !     theta    = dproduct(n, t, t)
    !     !     omega    = dproduct(n, t, s) / theta
    !     !     r        = axpby(n, alpha, p, omega, s)
    !     !     x        = xpby(n, x, 1.0_dp, r)
    !     !     r        = xpby(n, s, -omega, t)
    !     ! end do

    !     ! x = sgs_precond_solver(n, mtx, ind, x)



    ! end subroutine

    !===============================================================================================!
    !to perform sparse matrix vector multiplication Axb. A is a sparse square matrixin CSR format   !
    !===============================================================================================!

    pure function sp_matvec(nnod, mtx, ind, x) result(rx)

        integer, intent(in)                  :: nnod        ! length of diagonal
        type(matrix_elements), intent(in)    :: mtx
        type(matrix_index), intent(in)       :: ind
        real(dp), dimension(:), intent(in)   :: x       ! vector
        real(dp)                             :: rx(nnod)   ! resulting vector
     
        integer   :: n, j
        integer   :: row_start, row_end
        real(dp)  :: tmpsum
     

        do concurrent (n = 1:nnod)
            tmpsum = 0.
            row_start = ind % row(n)
            row_end   = ind % row(n+1) - 1
            do j = row_start, row_end
                tmpsum = tmpsum + mtx % elmn(j) * x(ind % col(j)) 
            end do
            rx(n) = tmpsum
        end do
     
    end function

    !===============================================================================================!
    !to calculate dot product of vectors a and b          !
    !===============================================================================================!

    pure function dproduct(nnod, a, b) result(rx)

        integer, intent(in)                  :: nnod
        real(dp), dimension(:), intent(in)   :: a, b   ! vector
        real(dp)                             :: rx     ! resulting vector

        real(dp)  :: x
        integer   :: n
        
        x = 0._dp
        do concurrent (n = 1:nnod)
            x = x + a(n) * b(n)
        end do
     
        rx = x
     
    end function

    !===============================================================================================!
    !to perform a*\vec{x}+b*\vec{y} calculation                                                     !
    !===============================================================================================!

    pure function axpby(nnod, alpha, x, beta, y) result(w)

        integer, intent(in)   :: nnod
        real(dp), intent(in)  :: alpha, beta, x(:), y(:)
        real(dp)              :: w(nnod)

        integer               :: n
      
        do concurrent (n = 1:nnod)
            w(n) = alpha*x(n) + beta*y(n)
        end do
     
    end function

    !===============================================================================================!
    !to perform \vec{x}+b*\vec{y} calculation                                                     !
    !===============================================================================================!

    pure function xpby(nnod, x, beta, y) result(w)

        integer, intent(in)   :: nnod
        real(dp), intent(in)  :: beta, x(:), y(:)
        real(dp)              :: w(nnod)

        integer               :: n
    
        do concurrent (n = 1:nnod)
            w(n) = x(n) + beta*y(n)
        end do
     
    end function

    !===============================================================================================!
    !to perform copy value x to y (x=y)                                                   !
    !===============================================================================================!

    pure function copy(nnod, x) result(y)

        integer, intent(in)    :: nnod
        real(dp), intent(in)   :: x(:)
        real(dp)               :: y(nnod)

        integer :: n
    
        do concurrent (n = 1:nnod)
            y(n) = x(n)
        end do
     
    end function

    !===============================================================================================!
    ! SGS preconditioner solver
    !===============================================================================================!

    pure function sgs_precond_solver(nnod, mtx, ind, y) result(x)

        integer, intent(in)                :: nnod
        type(matrix_elements), intent(in)  :: mtx
        type(matrix_index), intent(in)     :: ind
        real(dp), intent(in)               :: y(:)
        real(dp)                           :: x(nnod)

        real(dp) :: tmp(nnod)
        real(dp) :: sump
        integer  :: n, j
        integer   :: row_start, row_end
     

        ! tmp(1) = y(1) / mtx % elmn(ind % diag(1))
        ! do concurrent (n = 2:nnod)
        !     sump = 0.
        !     row_start   = ind % row(n)
        !     row_end     = ind % diag(n) - 1
        !     do j = row_start, row_end
        !         sump = sump + mtx % elmn(j) * y(ind % col(j))
        !     end do
        !     tmp(n) = (y(n) - sump) / mtx % elmn(ind % diag(n))
        ! enddo

        ! x(nnod) = tmp(nnod)
        ! do n = nnod-1,1,-1
        !     sump = 0.
        !     row_start   = ind % diag(n) + 1
        !     row_end     = ind % row(n+1) - 1
        !     do j = row_start, row_end
        !         sump = sump + mtx % elmn(j) * tmp(ind % col(j)) / mtx % elmn(ind % diag(n))
        !     end do
        !     x(n) = tmp(n) - sump
        ! enddo

        do n = 1, nnod
            x(n) = y(n) / mtx % elmn(ind % diag(n))
        end do
     
    end function


    ! pure function sgs_precond_multiply(nnod, mtx, ind, y) result(x)

    !     integer, intent(in)                :: nnod
    !     type(matrix_elements), intent(in)  :: mtx
    !     type(matrix_index), intent(in)     :: ind
    !     real(dp), intent(in)               :: y(:)
    !     real(dp)                           :: x(nnod)

    !     real(dp) :: tmp(nnod)
    !     real(dp) :: sump
    !     integer  :: n, j
    !     integer   :: row_start, row_end
     

        ! tmp(1) = y(1) / mtx % elmn(ind % diag(1))
        ! do concurrent (n = 2:nnod)
        !     sump = 0.
        !     row_start   = ind % row(n)
        !     row_end     = ind % diag(n) - 1
        !     do j = row_start, row_end
        !         sump = sump + mtx % elmn(j) * y(ind % col(j))
        !     end do
        !     tmp(n) = (y(n) - sump) / mtx % elmn(ind % diag(n))
        ! enddo

        ! x(nnod) = tmp(nnod)
        ! do n = nnod-1,1,-1
        !     sump = 0.
        !     row_start   = ind % diag(n) + 1
        !     row_end     = ind % row(n+1) - 1
        !     do j = row_start, row_end
        !         sump = sump + mtx % elmn(j) * tmp(ind % col(j)) / mtx % elmn(ind % diag(n))
        !     end do
        !     x(n) = tmp(n) - sump
        ! enddo

    !     do n = 1, nnod
    !         x(n) = y(n) * mtx % elmn(ind % diag(n))
    !     end do
     
    ! end function
    
end module
    