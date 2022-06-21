module gpu

    use sdata, only: dp
    implicit none
    save

    interface gpu_allocate
        procedure allocate_real_1d
        procedure allocate_real_2d
    end interface

    interface gpu_initialize
        procedure initialize_real_1d
        procedure initialize_real_scalar
    end interface

    contains

    subroutine allocate_real_1d(x, len1)

        real(dp), allocatable :: x(:)
        integer, intent(in)   :: len1
  
        allocate(x(len1))
        !$acc enter data create(x)
    
    end subroutine
  
    subroutine allocate_real_2d(x, len1, len2)

  
        real(dp), allocatable :: x(:,:)
        integer, intent(in)   :: len1, len2
    
        allocate(x(len1, len2))
        !$acc enter data create(x)
    
    end subroutine

    subroutine initialize_real_1d(x, value)

        real(dp)              :: x(:)
        real(dp), intent(in)  :: value
  
        !$acc kernels present(x)
        x(:) = value
        !$acc end kernels
    
    end subroutine

    subroutine initialize_real_scalar(x, value)

        real(dp)              :: x
        real(dp), intent(in)  :: value
  
        !$acc kernels present(x)
        x = value
        !$acc end kernels
    
    end subroutine
  
  
  end module
  