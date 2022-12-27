module data

    use constant
    use node_h,   only: node_type
  
    implicit none
  
    save
    
    character(len=length_line) :: mode
    
    integer :: ng                                        ! number of groups
    integer :: nmat                                      ! number of materials
    integer :: nnod                                      ! Number of nodes
  
    type(node_type), allocatable  :: node(:)             ! Node data
  
    ! Flux and fission source
    real(dp), allocatable, target :: flux(:,:)           ! node-averaged flux
    real(dp), allocatable         :: flux_prev(:,:)      ! previousnode-averaged flux
    real(dp), allocatable, target :: fsrc(:)             ! Fission source
    real(dp), allocatable         :: fsrc_prev(:)        ! Previous fission source
  
    ! Iteration Control
    real(dp) :: max_flux_error = 1.e-5    ! Flux Error Criteria
    real(dp) :: max_fsrc_error = 1.e-5    ! Fission source Error CRITERIA
    real(dp) :: flux_error                ! Flux and Fission source error in BCSEARCH calcs.
    real(dp) :: fsrc_error                ! Flux and Fission source error in BCSEARCH calcs.
    integer  :: n_inner = 2               ! Maximum inner iteration
    integer  :: n_outer = 500             ! Maximum outer iteration
    integer  :: extrp_interval = 5        ! Fission source extrapolation interval
    integer  :: n_th_iter = 30            ! Maximum number of thermal-hydraulics iteration
    integer  :: n_outer_th = 20           ! Maximum number of outer iterations per thermal-hydraulics iteration
    integer  :: upd_interval              ! Nodal update interval
    
    ! Output print option
    integer :: print_rad_pow=YES, print_axi_pow=YES, print_flux=YES
    
    ! Parameter references for TH feedback
    real(dp) :: fuel_temp_ref      ! Fuel temperature Reference in Kelvin
    real(dp) :: mod_temp_ref       ! Moderator temperature Reference in Kelvin
    real(dp) :: mod_dens_ref       ! Coolant Density Reference in g/cm3
  
    ! Boron Concentration
    real(dp) :: bcon       ! Boron concentration in ppm
    real(dp) :: rbcon      ! Boron concentration in ppm Reference
    
    ! Transient parameters
    integer, parameter    :: nf = 6                ! Number of delayed neutron precusor family
    real(dp)              :: small_theta = 1._dp   ! Small theta and big theta for transient using theta method
    real(dp)              :: big_theta = 0._dp     ! Small theta and big theta for transient using theta method
    real(dp)              :: beta(nf)              ! beta (delayed neutron fraction)
    real(dp)              :: lambda(nf)            ! lambda (precusor decay constant)
    real(dp)              :: core_beta             ! Core averaged beta
    real(dp), allocatable :: neutron_velo(:)       ! Neutron velocity
    real(dp)              :: total_time            ! TOTAL SIMULATION TIME
    real(dp)              :: time_step_1           ! FIRST TIME STEP
    real(dp)              :: time_step_2           ! SECOND TIME STEP
    real(dp)              :: time_mid              ! WHEN SECOND TIME STEP START
    logical :: transient_warning = .false.         ! To activate unconverged  outer iteration warning
    
    ! Thermal-hydraulics parameters
    real(dp)              :: total_power                   ! sum of node power
    real(dp)              :: pow                           ! Reactor power for given geometry (watt)
    real(dp)              :: ppow                          ! Reactor percent power in percent
    real(dp)              :: tpow                          ! Total reactor power
    real(dp)              :: tin                           ! coolant inlet temperature (kelvin)
    real(dp)              :: cflow                         ! Sub-channel mass flow rate (kg/s)
    real(dp)              :: rf, tg, tc, ppitch            ! Fuel meat radius, gap thickness, clad thickness, and pin picth (m)
    real(dp)              :: rg, rc                        ! Outer radius of gap and cladding
    real(dp)              :: dia, dh, farea                ! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
    real(dp)              :: cf                            ! heat fraction deposited into coolant
    real(dp), allocatable :: node_nf(:,:)                  ! Number of fuel pin per node
    integer, parameter    :: nm = 10                       ! number of Fuel meat mesh
    integer, parameter    :: nt = nm+2                     ! Number Total mesh (+2 mesh for gap and clad)
    real(dp), allocatable :: rdel(:)                       ! mesh delta
    real(dp), allocatable :: rpos(:)                       ! mesh position
    real(dp)              :: th_err                        ! Doppler error
    integer, parameter    :: thunit = 300                  ! Unit number to open steam table file
    
    ! Steam Table data
    integer, parameter:: ntem = 9   ! Number of temperature in steam table
    real(dp), dimension(ntem,6) :: stab  ! Steam table matrixs
    
    ! CMFD
    real(dp), dimension(:), allocatable   :: a1n, a2n, a3n, a4n  ! Nodal expansion coefficients for current node
    real(dp), dimension(:), allocatable   :: a1p, a2p, a3p, a4p  ! Nodal expansion coefficients for following node
    real(dp), dimension(:), allocatable   :: Ln1, Lp1            ! First Transverse leakages moments
    real(dp)                              :: ndmax               ! Maximum change in nodal coupling coefficients
    character(len=4)                      :: kern = 'SANM'
    integer                               :: im, jm, km
    real(dp), allocatable                 :: r(:), rs(:), v(:), p(:), s(:), t(:), tmp(:)
    
    !Timing
    real(dp) :: fdm_time = 0., nod_time = 0., xs_time = 0., &
    inp_time = 0., th_time = 0.
    
    contains
    
    function get_time() result (time)
    
        implicit none
    
        real(dp) :: time
    
        call cpu_time(time)
    
    end function

end module