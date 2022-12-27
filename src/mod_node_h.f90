module node_h

    use constant
    use xsec_h,   only : xs_type

    implicit none

    save

    type, private :: surface_type
        real(dp), allocatable     :: dt(:)           ! D tilde : FDM nodal coupling coefficients
        real(dp), allocatable     :: dh(:)           ! D hat   : Corrected (higher order) nodal coupling coefficients (X+,X-,Y+, Y-, Z+, Z-)
        real(dp), allocatable     :: adf(:)          ! Assembly Discontinuity Factors
    end type
    
    type :: node_type
        real(dp), pointer     :: flux(:)           ! node-averaged flux
        real(dp), pointer     :: fsrc(:)           ! fission source

        real(dp), allocatable :: xdel              ! Delta x in cm
        real(dp), allocatable :: ydel              ! Delta y in cm
        real(dp), allocatable :: zdel              ! Delta z in cm
        real(dp), allocatable :: vdel              ! Delta nodes' volume in cm3

        real(dp), allocatable :: fuel_temp         ! Fuel temperature
        real(dp), allocatable :: mod_temp          ! Moderator temperature
        real(dp), allocatable :: mod_dens          ! Moderator density
        real(dp), allocatable :: node_power        ! nodes power (watt)
        real(dp), allocatable :: tfm(:)            ! Fuel pin mesh temperature for each nodes
        real(dp), allocatable :: ent               ! Coolant Enthalpy (J/Kg)
        real(dp), allocatable :: heatf             ! Heat flux (W/m2
        real(dp), allocatable :: frate             ! coolant mass flow rate

        real(dp), allocatable :: prec_dens(:)      ! precusor_density
        real(dp), allocatable :: trans_src(:)      ! source for transient mode
        real(dp), allocatable :: ext_src(:)        ! external source
        real(dp), allocatable :: omega(:)          ! Exponential transformation constant
        real(dp), allocatable :: sigrp(:)          ! Initial removal cross sections before added by parameters required for transient
        real(dp), allocatable :: L    (:)          ! Total leakages for node n and group g
        real(dp), allocatable :: dfis               

        integer               :: i_mat               ! material index

        type(surface_type), allocatable :: sw(:)      ! surface west
        type(surface_type), allocatable :: se(:)      ! surface east
        type(surface_type), allocatable :: sn(:)      ! surface north
        type(surface_type), allocatable :: ss(:)      ! surface south
        type(surface_type), allocatable :: st(:)      ! surface top
        type(surface_type), allocatable :: sb(:)      ! surface bottom

        type(xs_type), pointer :: xs                  ! XSEC data

        type(node_type), pointer  :: west   => null()   ! Linked list to west node
        type(node_type), pointer  :: east   => null()   ! Linked list to east node
        type(node_type), pointer  :: north  => null()   ! Linked list to north node
        type(node_type), pointer  :: south  => null()   ! Linked list to south node
        type(node_type), pointer  :: top    => null()   ! Linked list to top node
        type(node_type), pointer  :: bottom => null()   ! Linked list to bottom node 
    end type

    
end module node_h
    