module geometry

    use constant
  
    implicit none
  
    save

    integer              :: nx, ny, nz                   ! Number of assemblies in x, y, and z directions
    integer              :: nxx, nyy, nzz                ! Number of nodes in x, y, and z directions
    integer, allocatable :: ix(:)                        ! index of x for given node index
    integer, allocatable :: iy(:)                        ! index of y for given node index
    integer, allocatable :: iz(:)                        ! index of z for given node index
    integer, allocatable :: xyz(:,:,:)                   ! index of node for given x, y, z index
    integer, allocatable :: xdiv(:)                      ! Assembly division in x direction
    integer, allocatable :: ydiv(:)                      ! Assembly division in y direction
    integer, allocatable :: zdiv(:)                      ! Assembly division in z direction
    integer              :: xwest, xeast                 ! West and East Boundary conditions
    integer              :: ysouth, ynorth               ! South and North Boundary conditions
    integer              :: zbott, ztop                  ! Bottom and Top Boundary conditions
   
    type :: STAGGERED
      integer :: smax, smin                              ! imax and imin along x and y direction for staggered nodes
    end type
    type(STAGGERED), dimension(:), allocatable :: ystag, xstag

    ! Below variables only used when reading input
    integer :: np                                               ! Number of planars
    integer, dimension(:), allocatable :: zpln                  ! Planar assignment to z direction
    real(dp), dimension(:), allocatable :: xsize, ysize, zsize  ! Assembly size for each direction
    
    type :: MAT_ASGN                                            ! Material assignment
        integer, dimension(:,:), allocatable :: asm             ! Material assignment into assembly
        integer, dimension(:,:), allocatable :: node            ! Material assignment into nodes
    end type
    type(MAT_ASGN), dimension(:), allocatable :: planar         ! planar
    integer, dimension(:,:,:), allocatable :: mnum

end module