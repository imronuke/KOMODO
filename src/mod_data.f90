module sdata

implicit none

save

integer, parameter :: dp = selected_real_kind(10, 15)

character(LEN=100) :: mode

integer :: ng     ! number of groups
integer :: nmat   ! number of materials

! XSECs Assigned to Nodes
real(dp), dimension(:,:), allocatable :: sigtr          ! Transport macroscopic XSEC
real(dp), dimension(:,:), allocatable :: siga           ! Absorption macroscopic XSEC
real(dp), dimension(:,:), allocatable :: nuf            ! nu* fission macroscopic XSEC
real(dp), dimension(:,:), allocatable :: sigf           ! fission macroscopic XSEC
real(dp), dimension(:,:), allocatable :: chi            ! neutron fission spectrum
real(dp), dimension(:,:,:), allocatable :: sigs         ! Scattering macroscopic XSEC
real(dp), dimension(:,:,:), allocatable :: dc           ! ADF
real(dp), dimension(:,:), allocatable :: D              ! Diffusion coefficient
real(dp), dimension(:,:), allocatable :: sigr           ! Removal macroscopic XSEC

! XSECs Assigned to Materials
real(dp), dimension(:,:), allocatable :: xsigtr          ! Transport macroscopic XSEC
real(dp), dimension(:,:), allocatable :: xsiga           ! Absorption macroscopic XSEC
real(dp), dimension(:,:), allocatable :: xnuf            ! nu* fission macroscopic XSEC
real(dp), dimension(:,:), allocatable :: xsigf           ! fission macroscopic XSEC
real(dp), dimension(:,:), allocatable :: xD              ! Diffusion coefficient
real(dp), dimension(:,:), allocatable :: xsigr           ! Removal macroscopic XSEC
real(dp), dimension(:,:,:), allocatable :: xsigs         ! Scattering macroscopic XSEC
logical :: ccnuf = .true.                            ! Logical variable to check the presence of fissile material
logical :: ccsigf = .true.                           ! Logical variable to check the presence of fissile material

! Geometry
integer :: nx, ny, nz                                ! Number of assemblies in x, y, and z directions
integer :: nxx, nyy, nzz                             ! Number of nodes in x, y, and z directions
integer :: nnod                                      ! Number of nodes
integer, dimension(:), allocatable :: ix, iy, iz
integer, dimension(:,:,:), allocatable :: xyz
integer, dimension(:), allocatable :: xdiv, ydiv, zdiv     ! Assembly division
real(dp), dimension(:), allocatable :: xdel, ydel, zdel, vdel  ! Delta x, y and z and nodes' volume in cm3
integer :: xwest, xeast, ysouth, ynorth, zbott, ztop       ! Boundary conditions
integer, dimension(:), allocatable :: mat                  ! Material assignment to nodes

! FDM Matrix (Stored in Compressed Sparse Row aka CSR)
type :: fdm_matr
  real(dp), dimension(:), allocatable :: elmn               ! Non-zero elements of FDM matrix for a row
end type
type(fdm_matr), dimension(:), allocatable :: A            ! FDM matrix
type :: fdm_ind
  integer, dimension(:), allocatable :: row                ! Row Pointer
  integer, dimension(:), allocatable :: col                ! Column index for the non-zero element of the FDM matrix
end type
type(fdm_ind) :: ind            ! Index of the FDM matrix

! Keff, flux and currents
real(dp) :: Ke
type :: node_data
  real(dp), dimension(6) :: df             ! FDM nodal coupling coefficients (X+,X-,Y+, Y-, Z+, Z-)
  real(dp), dimension(6) :: dn             ! Corrected (higher order) nodal coupling coefficients (X+,X-,Y+, Y-, Z+, Z-)
end type
type(node_data), dimension(:,:), allocatable :: nod

real(dp), dimension(:,:), allocatable :: f0, ft      ! current and previous Fluxes
real(dp), dimension(:), allocatable :: fs0, fst     ! Fission source
real(dp), dimension(:,:), allocatable :: c0 ! neutron precusor density
real(dp), dimension(:,:), allocatable :: s0 ! neutron precusor density

type :: staggered
    integer :: smax, smin                             ! imax and imin along x and y direction for staggered nodes
end type
type(staggered), dimension(:), allocatable :: ystag, xstag

! Extra Sources
real(dp), dimension(:,:), allocatable :: exsrc
real(dp) :: powtot

! Iteration Control
real(dp) :: ferc = 1.e-5    ! Flux Error Criteria
real(dp) :: serc = 1.e-5    ! Fission source Error CRITERIA
real(dp) :: fer, ser        ! Flux and Fission source error in BCSEARCH calcs.
integer :: nin = 2      ! Maximum inner iteration
integer :: nout = 500   ! Maximum outer iteration
integer :: nac = 5      ! Fission source extrapolation interval
integer :: th_niter = 30   ! Maximum number of thermal-hydraulics iteration
integer :: nth = 20     ! Maximum number of outer iterations per thermal-hydraulics iteration
integer :: nupd         ! Nodal update interval

character(16) :: matrix_solver = 'prec_cg' ! Linear system solver ('prec_cg', 'cg', 'bicg')
real(dp) :: inner_atol = 1d-8 ! Linear system solver absolute tolerance
real(dp) :: inner_rtol = 1d-5 ! Linear system solver relative tolerance

! OUTPUT PRINT OPTION
integer :: aprad=1, apaxi=1, afrad=1

! FUEL TEMPERATURE
real(dp), dimension(:), allocatable :: ftem       ! Fuel temperature in Kelvin for each nodes
real(dp) :: rftem      ! Fuel temperature Reference in Kelvin
real(dp), dimension(:,:), allocatable :: fsigtr, fsiga, fnuf, fsigf   ! XSEC changes per fuel temp changes
real(dp), dimension(:,:,:), allocatable :: fsigs

! MODERATOR TEMPERATURE
real(dp), dimension(:), allocatable :: mtem       ! Moderator temperature in Kelvin for each nodes
real(dp) :: rmtem      ! Moderator temperature Reference in Kelvin
real(dp), dimension(:,:), allocatable :: msigtr, msiga, mnuf, msigf   ! XSEC changes per Moderator temp changes
real(dp), dimension(:,:,:), allocatable :: msigs

! COOLANT DENSITY
real(dp), dimension(:), allocatable :: cden       ! Coolant Density in g/cm3 for each nodes
real(dp) :: rcden      ! Coolant Density Reference in g/cm3
real(dp), dimension(:,:), allocatable :: lsigtr, lsiga, lnuf, lsigf   ! XSEC changes per Coolant density changes
real(dp), dimension(:,:,:), allocatable :: lsigs

! Crod changes
integer :: nb                                                     ! Number of CR banks
real(dp), dimension(:), allocatable :: bpos  ! CR bank position
real(dp), dimension(:), allocatable :: fbpos    ! Final CR bank position
real(dp), dimension(:,:), allocatable :: dsigtr, dsiga, dnuf, dsigf   ! XSEC incerement or decrement due to CR insertion
real(dp), dimension(:,:,:), allocatable :: dsigs
real(dp), dimension(:,:,:), allocatable :: ddc   ! increment or decreent for ADF
real(dp), dimension(:), allocatable :: tmove    ! Time when CR bank starts moving
real(dp), dimension(:), allocatable :: bspeed   ! CR bank movement speed
integer, dimension(:), allocatable :: mdir  ! To indicate CR movement direction (0=do not move, 1=down, 2 = up)
real(dp) :: nstep                                         ! Number of steps
real(dp)    :: coreh                                      ! Core Height
integer, dimension(:,:), allocatable :: fbmap             ! Radial control rod bank map (node wise)
real(dp) :: pos0, ssize                                   ! Zero step position and step size

! Boron Concentration
real(dp) :: bcon       ! Boron concentration in ppm
real(dp) :: rbcon      ! Boron concentration in ppm Reference
real(dp), dimension(:,:), allocatable :: csigtr, csiga, cnuf, csigf   ! XSEC changes due to boron concentration
real(dp), dimension(:,:,:), allocatable :: csigs                      ! Used only for CBCS card

! Transient parameters
integer, parameter :: nf = 6                           ! Number of delayed neutron precusor family
real(dp)           :: sth = 1._dp, bth = 0._dp         ! Small theta and big theta for transient using theta method
real(dp), dimension(nf) :: ibeta, lamb                 ! beta (delayed neutron fraction) and precusor decay constant
real(dp), dimension(:), allocatable :: tbeta           ! total beta
real(dp) :: ctbeta                                     ! Core averaged
real(dp), dimension(:), allocatable :: velo            ! Neutron velocity
real(dp) :: ttot                                       ! TOTAL SIMULATION TIME
real(dp) :: tstep1                                     ! FIRST TIME STEP
real(dp) :: tstep2                                     ! SECOND TIME STEP
real(dp) :: tdiv                                       ! WHEN SECOND TIME STEP APPLY
real(dp), dimension(:,:), allocatable :: omeg          ! Exponential transformation constant
real(dp), dimension(:,:), allocatable :: sigrp         ! Initial removal cross sections before added by parameters required for transient
real(dp), dimension(:,:), allocatable :: L             ! Total leakages for node n and group g
real(dp), dimension(:), allocatable :: dfis
logical :: tranw = .false.                             ! To activate unconverged  outer iteration warning

! Thermal-hydraulics parameters
real(dp) :: pow                                        ! Reactor power for given geometry (watt)
real(dp) :: ppow                                       ! Reactor percent power in percent
real(dp) :: tpow                                       ! Total reactor power
real(dp), dimension(:), allocatable :: npow            ! nodes power (watt)
real(dp) :: tin                                        ! coolant inlet temperature (kelvin)
real(dp) :: cflow                                      ! Sub-channel mass flow rate (kg/s)
real(dp) :: rf, tg, tc, ppitch                         ! Fuel meat radius, gap thickness, clad thickness, and pin picth (m)
real(dp) :: rg, rc                                     ! Outer radius of gap and cladding
real(dp) :: dia, dh, farea                             ! Pin diameter, Hydraulic diameter (m) and sub-channel area (m2)
real(dp) :: cf                                         ! heat fraction deposited into coolant
real(dp), dimension(:,:), allocatable :: node_nf       ! Number of fuel pin per node
integer, parameter :: nm = 10                          ! number of Fuel meat mesh
integer, parameter :: nt = nm+2                        ! Number Total mesh (+2 mesh for gap and clad)
real(dp), dimension(:,:), allocatable :: tfm           ! Fuel pin mesh temperature for each nodes
real(dp), dimension(:), allocatable :: rdel            ! mesh delta
real(dp), dimension(:), allocatable :: rpos            ! mesh position
real(dp) :: th_err                                     ! Doppler error
real(dp), dimension(:), allocatable :: ent             ! Coolant Enthalpy (J/Kg)
real(dp), dimension(:), allocatable :: heatf           ! Heat flux (W/m2
real(dp), dimension(:), allocatable :: frate           ! coolant mass flow rate
integer, parameter :: thunit = 300                     ! Unit number to open steam table file
real(dp), parameter :: pi = 3.14159265

! Steam Table data
integer, parameter:: ntem = 9   ! Number of temperature in steam table
real(dp), dimension(ntem,6) :: stab  ! Steam table matrixs

! Data type for branch xsec data used if XTAB file present
type :: xbranch
  real(dp), dimension(:), allocatable :: sigtr, siga, nuf, sigf   !XSEC
  real(dp), dimension(:,:), allocatable :: sigs
  real(dp), dimension(:,:), allocatable :: dc        !ASSEMBLY DISCONTINUITY FACTOR
end type
! Data Type to store data in XTAB file
type :: mbranch
  integer :: nd, nb, nf, nm  ! BRANCH parameter dimension (Coolant dens., boron conc., fuel and moderator temp.)
  real(dp), dimension(:), allocatable :: pd, pb, pf, pm         !Branch paramaters (Coolant dens., boron conc., fuel and moderator temp.)
  type(xbranch), dimension(:,:,:,:), allocatable :: xsec        !Unrodded XSEC
  type(xbranch), dimension(:,:,:,:), allocatable :: rxsec       !Rodded XSEC
  real(dp), dimension(:), allocatable :: velo   ! Neutron velocity
  real(dp), dimension(nf) :: ibeta, lamb          ! beta and decay constant
  integer :: tadf            !Control input: adf
  integer :: trod            !Control input: control rod
end type
type(mbranch), dimension(:), allocatable :: m


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

! VTK format output
integer :: nver        !< Number of vertices

type :: vertice_type
    real(dp) :: x, y, z
    integer :: index
end type 

type :: cell_type
    type(vertice_type) :: vertice(8)
end type 

type(cell_type), allocatable :: cell(:)

!****************************************************************************!

contains

function get_time() result (time)

    real(dp) :: time

    call cpu_time(time)

end function

end module sdata
