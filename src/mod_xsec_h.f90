module xsec_h

    use constant

    implicit none

    save
    
    type :: xs_type
        real(dp), allocatable :: sigtr(:)          ! Transport macroscopic XSEC
        real(dp), allocatable :: siga (:)          ! Absorption macroscopic XSEC
        real(dp), allocatable :: nuf  (:)          ! nu*fission macroscopic XSEC
        real(dp), allocatable :: sigf (:)          ! fission macroscopic XSEC
        real(dp), allocatable :: D    (:)          ! Diffusion coefficient
        real(dp), allocatable :: sigr (:)          ! Removal macroscopic XSEC
        real(dp), allocatable :: sigs (:,:)        ! Scattering macroscopic XSEC
    end type

    type :: xs_type_mat
        real(dp), allocatable :: chi  (:)          ! neutron fission spectrum
        real(dp), allocatable :: beta_total  (:)   ! total delayed neutron for each precusor family
    end type

    ! Change of XS per change of other parameters
    type, private :: xs_changes
        real(dp), allocatable :: sigtr(:)          ! Transport macroscopic XSEC
        real(dp), allocatable :: siga (:)          ! Absorption macroscopic XSEC
        real(dp), allocatable :: nuf  (:)          ! nu*fission macroscopic XSEC
        real(dp), allocatable :: sigf (:)          ! fission macroscopic XSEC
    end type

    type, extends(xs_changes), private :: xs_changes_tab  ! if xtab file is used
        real(dp), allocatable :: adf (:,:)
    end type

    type :: xs_change_type
        type(xs_changes)  :: fuel_temp
        type(xs_changes)  :: mod_temp
        type(xs_changes)  :: mod_dens
        type(xs_changes)  :: rods
        type(xs_changes)  :: boron
    end type

    type :: rod_xs_change_type
        type(xs_changes_tab)  :: rod
    end type


    
end module
    