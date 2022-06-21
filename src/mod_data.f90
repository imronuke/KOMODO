MODULE sdata

save

integer, parameter :: dp = selected_real_kind(10, 15)

CHARACTER(LEN=100) :: mode

INTEGER :: ng     ! number of groups
INTEGER :: nmat   ! number of materials

! XSECs Assigned to Nodes
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sigtr          ! Transport macroscopic XSEC
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: siga           ! Absorption macroscopic XSEC
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: nuf            ! nu* fission macroscopic XSEC
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sigf           ! fission macroscopic XSEC
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: chi            ! neutron fission spectrum
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: sigs         ! Scattering macroscopic XSEC
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: dc           ! ADF
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: D              ! Diffusion coefficient
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sigr           ! Removal macroscopic XSEC

! XSECs Assigned to Materials
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xsigtr          ! Transport macroscopic XSEC
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xsiga           ! Absorption macroscopic XSEC
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xnuf            ! nu* fission macroscopic XSEC
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xsigf           ! fission macroscopic XSEC
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xD              ! Diffusion coefficient
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xsigr           ! Removal macroscopic XSEC
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: xsigs         ! Scattering macroscopic XSEC
LOGICAL :: ccnuf = .TRUE.                            ! Logical variable to check the presence of fissile material
LOGICAL :: ccsigf = .TRUE.                           ! Logical variable to check the presence of fissile material

! Geometry
INTEGER :: nx, ny, nz                                ! Number of assemblies in x, y, and z directions
INTEGER :: nxx, nyy, nzz                             ! Number of nodes in x, y, and z directions
INTEGER :: nnod                                      ! Number of nodes
INTEGER, DIMENSION(:), ALLOCATABLE :: ix, iy, iz
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: xyz
INTEGER, DIMENSION(:), ALLOCATABLE :: xdiv, ydiv, zdiv     ! Assembly division
REAL(DP), DIMENSION(:), ALLOCATABLE :: xdel, ydel, zdel, vdel  ! Delta x, y and z and nodes' volume in cm3
INTEGER :: xwest, xeast, ysouth, ynorth, zbott, ztop       ! Boundary conditions
INTEGER, DIMENSION(:), ALLOCATABLE :: mat                  ! Material assignment to nodes

! FDM Matrix (Stored in Compressed Sparse Row aka CSR)
TYPE :: FDM_MATR
  REAL(DP), DIMENSION(:), ALLOCATABLE :: elmn               ! Non-zero elements of FDM matrix for a row
END TYPE
TYPE(FDM_MATR), DIMENSION(:), ALLOCATABLE :: A            ! FDM matrix
TYPE :: FDM_IND
  INTEGER, DIMENSION(:), ALLOCATABLE :: row                ! Row Pointer
  INTEGER, DIMENSION(:), ALLOCATABLE :: col                ! Column index for the non-zero element of the FDM matrix
END TYPE
TYPE(FDM_IND) :: ind            ! Index of the FDM matrix

! Keff, flux and currents
REAL(DP) :: Ke
TYPE :: NODE_DATA
  REAL(DP), DIMENSION(6) :: df             ! FDM nodal coupling coefficients (X+,X-,Y+, Y-, Z+, Z-)
  REAL(DP), DIMENSION(6) :: dn             ! Corrected (higher order) nodal coupling coefficients (X+,X-,Y+, Y-, Z+, Z-)
END TYPE
TYPE(NODE_DATA), DIMENSION(:,:), ALLOCATABLE :: nod

REAL(DP), DIMENSION(:,:), ALLOCATABLE :: f0, ft      ! current and previous Fluxes
REAL(DP), DIMENSION(:), ALLOCATABLE :: fs0, fst     ! Fission source
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: c0 ! neutron precusor density
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: s0 ! neutron precusor density

TYPE :: STAGGERED
    INTEGER :: smax, smin                             ! imax and imin along x and y direction for staggered nodes
END TYPE
TYPE(STAGGERED), DIMENSION(:), ALLOCATABLE :: ystag, xstag

! Extra Sources
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: exsrc
REAL(DP) :: powtot

! Iteration Control
REAL(DP) :: ferc = 1.e-5    ! Flux Error Criteria
REAL(DP) :: serc = 1.e-5    ! Fission source Error CRITERIA
REAL(DP) :: fer, ser        ! Flux and Fission source error in BCSEARCH calcs.
INTEGER :: nin = 2      ! Maximum inner iteration
INTEGER :: nout = 500   ! Maximum outer iteration
INTEGER :: nac = 5      ! Fission source extrapolation interval
INTEGER :: th_niter = 30   ! Maximum number of thermal-hydraulics iteration
INTEGER :: nth = 20     ! Maximum number of outer iterations per thermal-hydraulics iteration
integer :: nupd         ! Nodal update interval

! OUTPUT PRINT OPTION
INTEGER :: aprad=1, apaxi=1, afrad=1

! FUEL TEMPERATURE
REAL(DP), DIMENSION(:), ALLOCATABLE :: ftem       ! Fuel temperature in Kelvin for each nodes
REAL(DP) :: rftem      ! Fuel temperature Reference in Kelvin
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: fsigtr, fsiga, fnuf, fsigf   ! XSEC changes per fuel temp changes
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: fsigs

! MODERATOR TEMPERATURE
REAL(DP), DIMENSION(:), ALLOCATABLE :: mtem       ! Moderator temperature in Kelvin for each nodes
REAL(DP) :: rmtem      ! Moderator temperature Reference in Kelvin
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: msigtr, msiga, mnuf, msigf   ! XSEC changes per Moderator temp changes
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: msigs

! COOLANT DENSITY
REAL(DP), DIMENSION(:), ALLOCATABLE :: cden       ! Coolant Density in g/cm3 for each nodes
REAL(DP) :: rcden      ! Coolant Density Reference in g/cm3
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: lsigtr, lsiga, lnuf, lsigf   ! XSEC changes per Coolant density changes
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: lsigs

! Crod changes
INTEGER :: nb                                                     ! Number of CR banks
REAL(DP), DIMENSION(:), ALLOCATABLE :: bpos  ! CR bank position
REAL(DP), DIMENSION(:), ALLOCATABLE :: fbpos    ! Final CR bank position
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dsigtr, dsiga, dnuf, dsigf   ! XSEC incerement or decrement due to CR insertion
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: dsigs
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: ddc   ! increment or decreent for ADF
REAL(DP), DIMENSION(:), ALLOCATABLE :: tmove    ! Time when CR bank starts moving
REAL(DP), DIMENSION(:), ALLOCATABLE :: bspeed   ! CR bank movement speed
INTEGER, DIMENSION(:), ALLOCATABLE :: mdir  ! To indicate CR movement direction (0=do not move, 1=down, 2 = up)
REAL(DP) :: nstep                                         ! Number of steps
REAL(DP)    :: coreh                                      ! Core Height
INTEGER, DIMENSION(:,:), ALLOCATABLE :: fbmap             ! Radial control rod bank map (node wise)
REAL(DP) :: pos0, ssize                                   ! Zero step position and step size

! Boron Concentration
REAL(DP) :: bcon       ! Boron concentration in ppm
REAL(DP) :: rbcon      ! Boron concentration in ppm Reference
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: csigtr, csiga, cnuf, csigf   ! XSEC changes due to boron concentration
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: csigs                      ! Used only for CBCS card

! Transient parameters
INTEGER, PARAMETER :: nf = 6                           ! Number of delayed neutron precusor family
REAL(DP)           :: sth = 1._dp, bth = 0._dp         ! Small theta and big theta for transient using theta method
REAL(DP), DIMENSION(nf) :: ibeta, lamb                 ! beta (delayed neutron fraction) and precusor decay constant
REAL(DP), DIMENSION(:), ALLOCATABLE :: tbeta           ! total beta
REAL(DP) :: ctbeta                                     ! Core averaged
REAL(DP), DIMENSION(:), ALLOCATABLE :: velo            ! Neutron velocity
REAL(DP) :: ttot                                       ! TOTAL SIMULATION TIME
REAL(DP) :: tstep1                                     ! FIRST TIME STEP
REAL(DP) :: tstep2                                     ! SECOND TIME STEP
REAL(DP) :: tdiv                                       ! WHEN SECOND TIME STEP APPLY
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: omeg          ! Exponential transformation constant
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sigrp         ! Initial removal cross sections before added by parameters required for transient
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: L             ! Total leakages for node n and group g
REAL(DP), DIMENSION(:), ALLOCATABLE :: dfis
LOGICAL :: tranw = .FALSE.                             ! To activate unconverged  outer iteration warning

! Thermal-hydraulics parameters
REAL(DP) :: pow                                        ! Reactor power for given geometry (watt)
REAL(DP) :: ppow                                       ! Reactor percent power in percent
REAL(DP) :: tpow                                       ! Total reactor power
REAL(DP), DIMENSION(:), ALLOCATABLE :: npow            ! nodes power (watt)
REAL(DP) :: tin                                        ! coolant inlet temperature (kelvin)
REAL(DP) :: cflow                                      ! Sub-channel mass flow rate (kg/s)
REAL(DP) :: rf, tg, tc, ppitch                         ! Fuel meat radius, gap thickness, clad thickness, and pin picth (m)
REAL(DP) :: rg, rc                                     ! Outer radius of gap and cladding
REAL(DP) :: dia, dh, farea                             ! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
REAL(DP) :: cf                                         ! heat fraction deposited into coolant
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: node_nf       ! Number of fuel pin per node
INTEGER, PARAMETER :: nm = 10                          ! number of Fuel meat mesh
INTEGER, PARAMETER :: nt = nm+2                          ! Number Total mesh (+2 mesh for gap and clad)
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: tfm           ! Fuel pin mesh temperature for each nodes
REAL(DP), DIMENSION(:), ALLOCATABLE :: rdel            ! mesh delta
REAL(DP), DIMENSION(:), ALLOCATABLE :: rpos            ! mesh position
REAL(DP) :: th_err                                     ! Doppler error
REAL(DP), DIMENSION(:), ALLOCATABLE :: ent             ! Coolant Enthalpy (J/Kg)
REAL(DP), DIMENSION(:), ALLOCATABLE :: heatf           ! Heat flux (W/m2
REAL(DP), DIMENSION(:), ALLOCATABLE :: frate           ! coolant mass flow rate
INTEGER, PARAMETER :: thunit = 300                 ! Unit number to open steam table file
REAL(DP), PARAMETER :: pi = 3.14159265

! Steam Table data
INTEGER, PARAMETER:: ntem = 9   ! Number of temperature in steam table
REAL(DP), DIMENSION(ntem,6) :: stab  ! Steam table matrixs

! Data type for branch xsec data used if XTAB file present
TYPE :: XBRANCH
  REAL(DP), DIMENSION(:), ALLOCATABLE :: sigtr, siga, nuf, sigf   !XSEC
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sigs
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dc        !ASSEMBLY DISCONTINUITY FACTOR
END TYPE
! Data Type to store data in XTAB file
TYPE :: MBRANCH
  INTEGER :: nd, nb, nf, nm  ! BRANCH PARAMETER DIMENSION (Coolant dens., boron conc., fuel and moderator temp.)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: pd, pb, pf, pm         !Branch paramaters (Coolant dens., boron conc., fuel and moderator temp.)
  TYPE(XBRANCH), DIMENSION(:,:,:,:), ALLOCATABLE :: xsec        !Unrodded XSEC
  TYPE(XBRANCH), DIMENSION(:,:,:,:), ALLOCATABLE :: rxsec       !Rodded XSEC
  REAL(DP), DIMENSION(:), ALLOCATABLE :: velo   ! Neutron velocity
  REAL(DP), DIMENSION(nf) :: ibeta, lamb          ! beta and decay constant
  INTEGER :: tadf            !Control input: adf
  INTEGER :: trod            !Control input: control rod
END TYPE
TYPE(MBRANCH), DIMENSION(:), ALLOCATABLE :: m


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

END MODULE sdata
