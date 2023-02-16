module data


    use iso_fortran_env, only: real64, output_unit
    use stagger

    implicit none
  
    public
    save

    ! Some constants
    integer, parameter, private  :: dp = real64
    real(dp), parameter          :: pi = acos(-1.0)
    integer , parameter          :: length_word = 30
    integer , parameter          :: length_line = 200
    integer , parameter          :: length_card = 4
    integer, parameter           :: YES = 1
    integer, parameter           :: NO = 0
    integer, parameter           :: ZERO_FLUX = 0
    integer, parameter           :: ZERO_INCOMING = 1
    integer, parameter           :: REFLECTIVE = 2

    integer :: ounit                      ! output file unit number

    ! Card active/inactive indicator (active = 1, inactive = 0)
    integer :: bmode = NO, bxsec = NO, bgeom = NO, bcase = NO, besrc = NO
    integer :: biter = NO, bprnt = NO, badf  = NO, bcrod = NO, bbcon = NO
    integer :: bftem = NO, bmtem = NO, bcden = NO, bcbcs = NO, bejct = NO
    integer :: bther = NO, bxtab = NO, bkern = NO, bextr = NO, bthet = NO
    integer :: boutp = NO
    
    character(len=length_word) :: mode    ! calculation mode

    ! Flux and fission source
    real(dp), allocatable :: flux(:,:)           ! node-averaged flux
    real(dp), allocatable :: fsrc(:)             ! node-averaged fission source
    real(dp), allocatable :: prev_flux(:,:)      ! previous time-step node-averaged flux
    real(dp), allocatable :: prev_fsrc(:)        ! previous time-step node-averaged fission source
    real(dp), allocatable :: exsrc(:,:)          ! node-averaged external source

    ! Material wise XS [dimension(nmat, ng)]
    integer               :: ng                  ! number of energy group
    integer               :: nmat                ! number of materials
    real(dp), allocatable :: sigtr(:,:)          ! Transport macroscopic XSEC
    real(dp), allocatable :: siga (:,:)          ! Absorption macroscopic XSEC
    real(dp), allocatable :: nuf  (:,:)          ! nu*fission macroscopic XSEC
    real(dp), allocatable :: sigf (:,:)          ! fission macroscopic XSEC
    real(dp), allocatable :: sigs (:,:,:)        ! Scattering macroscopic XSEC
    real(dp), allocatable :: chi  (:,:)          ! neutron fission spectrum
    real(dp), allocatable :: dc   (:,:,:)        ! ADF

    ! Geometry
    integer                       :: nnod                ! Total number of nodes
    integer                       :: nx, ny, nz          ! Number of assemblies in x, y, and z directions
    integer                       :: nxx, nyy, nzz       ! Number of nodes in x, y, and z directions
    real(dp), pointer             :: xdel(:)             ! Delta x in cm
    real(dp), pointer             :: ydel(:)             ! Delta y in cm
    real(dp), pointer             :: zdel(:)             ! Delta z in cm
    real(dp), pointer             :: vdel(:)             ! node volume cm3
    integer, allocatable          :: xdiv(:)             ! Assembly division in x direction
    integer, allocatable          :: ydiv(:)             ! Assembly division in y direction
    integer, allocatable          :: zdiv(:)             ! Assembly division in z direction
    real(dp)                      :: core_height         ! core height
    type(staggered), allocatable  :: xstag(:)            ! Staggered mesh data
    type(staggered), allocatable  :: ystag(:)            ! Staggered mesh data
    integer, allocatable          :: ix(:)               ! index in x direction for given node index
    integer, allocatable          :: iy(:)               ! index in y direction for given node index
    integer, allocatable          :: iz(:)               ! index in z direction for given node index
    integer, allocatable          :: xyz(:,:,:)          ! node index for given index in x, y, and z directions
    integer, allocatable          :: ind_mat(:)          ! node-wise material index (material map)
    integer                       :: east, west          ! BC
    integer                       :: north, south        ! 0 = zero flux, 1 = zero incoming current
    integer                       :: bottom, top         ! 2 = Reflective

    ! Thermal-hydraulics parameters
    real(dp), allocatable :: ftem(:)                ! Fuel temperature in Kelvin for each nodes
    real(dp), allocatable :: mtem(:)                ! Moderator temperature in Kelvin for each nodes
    real(dp), allocatable :: cden(:)                ! Coolant Density in g/cm3 for each nodes
    real(dp)              :: bcon                   ! Boron concentration in ppm
    real(dp)              :: total_power            ! sum of node power
    real(dp)              :: power                  ! Reactor power for given geometry (watt)
    real(dp)              :: percent_pow            ! Reactor percent power in percent
    real(dp)              :: tpow                   ! Total reactor power
    real(dp), allocatable :: npow(:)                ! nodes power (watt)
    real(dp)              :: t_inlet                ! coolant inlet temperature (kelvin)
    real(dp)              :: rf, tg, tc             ! Fuel meat radius, gap thickness, clad thickness (m)
    real(dp)              :: ppitch                 ! pin picth (m)
    real(dp)              :: cf                     ! heat fraction deposited into coolant
    real(dp)              :: cflow                  ! Sub-channel mass flow rate (kg/s)
    real(dp), allocatable :: node_nf(:,:)           ! Number of fuel pin per node
    real(dp), allocatable :: tfm(:,:)               ! Fuel pin mesh/ring temperature for each nodes
    real(dp)              :: th_err                 ! Doppler error
    real(dp), allocatable :: enthalpy  (:)          ! Coolant Enthalpy (J/Kg)
    real(dp), allocatable :: heat_flux (:)          ! sub-channel heat flux (W/m2
    real(dp), allocatable :: flow_rate(:)           ! sub-channel coolant mass flow rate(Kg/s)
    
    ! Crod changes
    integer               :: nbank                  ! Number of CR banks
    real(dp), allocatable :: bank_pos(:)            ! CR bank position
    real(dp), allocatable :: bank_pos_final(:)      ! Final CR bank position
    real(dp), allocatable :: t_move (:)             ! Time when CR bank starts moving
    real(dp), allocatable :: bank_speed(:)          ! CR bank movement speed
    integer,  allocatable :: direction(:)           ! To indicate CR movement direction (0=do not move, 1=down, 2 = up)
    real(dp)              :: nstep                  ! Number of steps
    real(dp)              :: coreh                  ! Core Height
    integer, allocatable  :: fbmap(:,:)             ! Radial control rod bank map (node wise)
    real(dp)              :: zero_pos, step_size    ! Zero step position and step size
  
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
    
    ! parameter references for TH feedback
    ! real(dp) :: fuel_temp_ref      ! Fuel temperature Reference in Kelvin
    ! real(dp) :: mod_temp_ref       ! Moderator temperature Reference in Kelvin
    ! real(dp) :: mod_dens_ref       ! Coolant Density Reference in g/cm3

    
    ! Transient parameters
    integer, parameter    :: nf = 6                ! Number of delayed neutron precusor family
    real(dp)              :: small_theta = 1._dp   ! Theta transient using theta method
    real(dp)              :: beta(nf)              ! beta (delayed neutron fraction)
    real(dp)              :: lambda(nf)            ! lambda (precusor decay constant)
    real(dp), allocatable :: total_beta(:)         ! beta total for each material
    real(dp), allocatable :: dfis(:)               
    real(dp)              :: core_beta             ! Core averaged beta
    real(dp), allocatable :: neutron_velo(:)       ! Neutron velocity
    real(dp)              :: total_time            ! TOTAL SIMULATION TIME
    real(dp)              :: time_step_1           ! FIRST TIME STEP
    real(dp)              :: time_step_2           ! SECOND TIME STEP
    real(dp)              :: time_mid              ! WHEN SECOND TIME STEP START
    logical :: transient_warning = .false.         ! To activate unconverged  outer iteration warning
    
    ! CMFD
    character(len=4)                    :: kern = 'SANM'   ! Nodal kernel

end module