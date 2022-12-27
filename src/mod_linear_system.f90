module linear_system

    use constant

    implicit none

    save
    
    ! Linear system Matrix (Stored in Compressed Sparse Row aka CSR)
    type :: matrix_type
        real(dp), allocatable :: elmn(:)               ! Non-zero elements of FDM matrix for a row
    end type
    type(matrix_type), allocatable :: A(:)             ! FDM matrix for each energy group

    ! Matrix indexing
    type :: matrix_index
        integer, dimension(:), allocatable :: row      ! Row pointer
        integer, dimension(:), allocatable :: col      ! Column index for the non-zero element of the FDM matrix
    end type
    type(matrix_index) :: ind                          ! Index of the FDM matrix

    
end module
    