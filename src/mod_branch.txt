module branch

    use constant

    implicit none

    save
    
    ! Data type for branch xsec data used if XTAB file present
    type :: XBRANCH
        real(dp), dimension(:), allocatable :: sigtr, siga, nuf, sigf   !XSEC
        real(dp), dimension(:,:), allocatable :: sigs
        real(dp), dimension(:,:), allocatable :: dc        !ASSEMBLY DISCONTINUITY FACTOR
    end type
    ! Data type to store data in XTAB file
    type :: MBRANCH
        integer :: nd, nb, nf, nm  ! BRANCH parameter DIMENSION (Coolant dens., boron conc., fuel and moderator temp.)
        real(dp), dimension(:), allocatable :: pd, pb, pf, pm         !Branch paramaters (Coolant dens., boron conc., fuel and moderator temp.)
        type(XBRANCH), dimension(:,:,:,:), allocatable :: xsec        !Unrodded XSEC
        type(XBRANCH), dimension(:,:,:,:), allocatable :: rxsec       !Rodded XSEC
        real(dp), dimension(:), allocatable :: velo   ! Neutron velocity
        real(dp), dimension(nf) :: ibeta, lamb          ! beta and decay constant
        integer :: tadf            !Control input: adf
        integer :: trod            !Control input: control rod
    end type
    type(MBRANCH), dimension(:), allocatable :: m

    
end module constant
    