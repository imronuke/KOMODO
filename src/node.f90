module node

    use iso_fortran_env, only: real64

    implicit none

    private

    save

    integer, parameter :: dp = real64

    integer, parameter  :: n_rect_surf = 6         ! order is => east, west, north, top, bottom

    ! Staggered mesh
    type, public :: staggered
        integer :: smax, smin                      ! imax and imin along x and y direction for staggered nodes
    end type

    ! Surface averaged values
    type, public :: surface_type
        real(dp), allocatable     :: dt(:)           ! D tilde : FDM nodal coupling coefficients
        real(dp), allocatable     :: dh(:)           ! D hat   : Corrected (higher order) nodal coupling coefficients (X+,X-,Y+, Y-, Z+, Z-)
    end type
    
    ! Base Node
    type :: node_base
        integer  :: n                 ! node index
    end type

    ! Rectangular node
    type, extends(node_base), public :: node_rect          

        type(surface_type)  :: sf(n_rect_surf)     ! surface
    
        type(node_rect), pointer  :: west   => null()   ! Linked list to west node
        type(node_rect), pointer  :: east   => null()   ! Linked list to east node
        type(node_rect), pointer  :: north  => null()   ! Linked list to north node
        type(node_rect), pointer  :: south  => null()   ! Linked list to south node
        type(node_rect), pointer  :: top    => null()   ! Linked list to top node
        type(node_rect), pointer  :: bottom => null()   ! Linked list to bottom node
    end type

end module
    