module th

    use iso_fortran_env, only: real64, error_unit
    use utilities

    implicit none

    private

    save

    integer, parameter    :: dp = real64

    real(dp), parameter   :: pi = acos(-1.0)
    integer               :: fid                         ! File output ID
    real(dp), parameter   :: fdens = 10.412e3            ! UO2 density (kg/m3)
    real(dp), parameter   :: cdens = 6.6e3               ! Cladding density (kg/m3)
    real(dp), parameter   :: hg    = 1.0e4               ! Gap heat transfer coef

    ! Steam Table data for water at 15.5 MPa
    integer, parameter    :: ntem = 9                                 ! Number of temperature in steam table
    real(dp), allocatable :: stab(:,:)
    
    ! Thermal-hydraulics data
    type, public :: th_type
        integer :: nzz
        integer :: n_fuel     = 10                      ! number of Fuel meat mesh
        integer :: n_gap_clad = 2                       ! number of ring in gap and cladding (one each)
        integer :: n_ring                               ! Number Total mesh (+2 mesh for gap and clad)
        
        real(dp)              :: alpha = 0.7            ! parameter used to determine lumped fuel temp
        real(dp)              :: cflow                  ! Sub-channel mass flow rate (kg/s)
        real(dp)              :: s_area                 ! sub-channel area (m2)
        real(dp)              :: dia                    ! Pin diameter (m)
        real(dp)              :: cf                     ! heat fraction deposited into coolant
        real(dp)              :: dh                     ! sub-channel hydraulic diameter (m)
        real(dp)              :: t_inlet                ! sub-channel inlet temperature (K)
        real(dp)              :: rf                     ! Fuel meat radius (m)
        real(dp)              :: rg, rc                 ! Outer radius of gap and cladding (m)
        real(dp), allocatable :: rdel(:)                ! mesh delta (cm)
        real(dp), allocatable :: rpos(:)                ! mesh position (m)
    
        real(dp), pointer :: zdel(:)
    end type

    public :: set_th_data, th_update, th_transient

    contains

    !===============================================================================================!
    ! Set th data
    !===============================================================================================!

    subroutine set_th_data(t, fid_inp, nzz, zdel, rf, tg, tc, ppitch, cflow, cf, t_inlet)

        class(th_type) :: t
        integer, intent(in)          :: fid_inp
        integer, intent(in)          :: nzz
        real(dp), intent(in), target :: zdel(:)
        real(dp), intent(in)         :: rf, tg, tc, ppitch
        real(dp), intent(in)         :: cflow
        real(dp), intent(in)         :: cf, t_inlet

        fid = fid_inp
        t % nzz  = nzz
        t % zdel => zdel
        if (size(t % zdel) .ne. nzz) stop 'zdel size in th module does not match'

        t % rf      = rf
        t % rg      = rf + tg
        t % rc      = t % rg + tc
        t % dia     = 2. * t % rc
        t % dh      = t % dia * ((4./pi) * (ppitch/t % dia)**2 - 1.)
        t % s_area  = ppitch**2 - 0.25*pi*t % dia**2
        t % cflow   = cflow
        t % t_inlet = t_inlet
        t % cf      = cf

        call th_setup(t, tg, tc)
        call set_steam_table()
        
    end subroutine

    !===============================================================================================!
    ! Do final TH setup
    !===============================================================================================!

    subroutine th_setup(t, tg, tc)

        class(th_type)       :: t
        real(dp), intent(in) :: tg, tc

        real(dp) :: dummy
        integer  :: i
        
        ! Total number of rings
        t % n_ring = t % n_fuel + t % n_gap_clad

        ! Calculate fuel pin mesh delta and position
        dummy = t % rf / real(t % n_fuel)
        
        allocate(t % rdel(t % n_ring), t % rpos(t % n_ring))
        do i = 1, t % n_fuel
            t % rdel(i) = dummy
        end do
        
        ! Fuel pin mesh size
        t % rdel(t % n_fuel+1) = tg
        t % rdel(t % n_fuel+2) = tc
        
        ! Fuel pin mesh position
        t % rpos(1) = 0.5 * t % rdel(1)
        do i = 2, t % n_ring
            t % rpos(i) = t % rpos(i-1) + 0.5 * (t % rdel(i-1) + t % rdel(i))
        end do
        
    end subroutine

    !===============================================================================================!
    ! set steam table
    !===============================================================================================!

    subroutine set_steam_table()

        allocate(stab(ntem, 6))
        stab = transpose(reshape( [  &  ! Steam table matrixs
        !Temp     Density    Enthalpy   Prandtl    Kinematic Visc.   Thermal Conduct.
        ! oC      (g/cm3)     (J/kg)                  1e-6 m2/s          W/(mK)
        543.15, 0.78106745,  1182595.0, 0.820773,     0.128988,         0.60720, &
        553.15, 0.76428125,  1232620.0, 0.829727,     0.126380,         0.59285, &
        563.15, 0.74619716,  1284170.0, 0.846662,     0.124118,         0.57710, &
        573.15, 0.72650785,  1337630.0, 0.863597,     0.121856,         0.55970, &
        583.15, 0.70475081,  1393570.0, 0.915035,     0.120105,         0.54045, &
        593.15, 0.68018488,  1452895.0, 0.966472,     0.118354,         0.51880, &
        603.15, 0.65150307,  1517175.0, 1.166745,     0.143630,         0.49420, &
        613.15, 0.61590149,  1589770.0, 1.515852,     0.195931,         0.46550, &
        617.91, 0.59896404,  1624307.1, 1.681940,     0.220813,         0.45185  &
        ],[6, ntem]))
        
    end subroutine

    !===============================================================================================!
    ! To get enthalpy from steam table for given coolant temp.
    !===============================================================================================!

    subroutine getent(temperature, enthalpy)
        
        real(dp), intent(in) :: temperature
        real(dp), intent(out) :: enthalpy
        real(dp) :: t1, ent1
        real(dp) :: t2, ent2
        integer :: i
        
        if ((temperature < stab(1,1)) .OR. (temperature > stab(ntem,1))) then
            call fatal_error(fid, &
            'Coolant temp. : ' // n2c(temperature, 3) // ' K' // new_line('a') // &
            'ERROR : MODERATOR TEMP. IS out OF THE RANGE OF DATA in THE STEAM TABLE' // new_line('a') // &
            'CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER' )
        end if
        
        t2 = stab(1,1); ent2 = stab(1,3)
        do i = 2, ntem
            t1 = t2
            ent1 = ent2
            t2 = stab(i,1); ent2 = stab(i,3)
            if ((temperature >= t1) .and. (temperature <= t2)) then
                enthalpy = ent1 + (temperature - t1) / (t2 - t1) * (ent2 - ent1)
                exit
            end if
        end do
        
    end subroutine getent

    !===============================================================================================!
    ! To get coolant properties from steam table for given enthalpy.
    !===============================================================================================!

    subroutine gettd(ent,temp,rho,prx,kvx,tcx,Rx)
        
        real(dp), intent(in) :: ent
        real(dp), intent(out) :: temp   ! Temperature
        real(dp), intent(out) :: rho    ! coolant densisty
        real(dp), intent(out) :: prx    ! Prandtl number
        real(dp), intent(out) :: kvx    ! kinematic viscosity
        real(dp), intent(out) :: tcx    ! coolant thermal conductivity
        real(dp), intent(out), OPTIONAL :: Rx
        real(dp) :: ratx
        
        integer :: i, i1, i2
        
        ! Get two closest interpolation points
        if (ent >= stab(1,3) .and. ent <= stab(ntem,3)) then  !if enthalpy inside data range
            do i = 2, ntem
                if (ent >= stab(i-1,3) .and. ent <= stab(i,3)) then
                    i1 = i-1
                    i2 = i
                    exit
                end if
            end do
        else if (ent < stab(1,3)  .and. (stab(1,3) - ent) / stab(1,3) < 0.1) then !if 10% lower than min. steam table data
            i1 = 1
            i2 = 2
        else if (ent > stab(ntem,3) .and. (ent - stab(ntem,3)) / stab(ntem,3) < 0.1) then !if 10% higher than max. steam table data
            i1 = ntem-1
            i2 = ntem
        else
            call fatal_error(fid, &
            'CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER' // new_line('a') // &
            'ERROR: ENTHALPY '// n2c(ent/1000., 1) //' KJ/Kg IS out OF THE RANGE in THE STEAM TABLE')
        end if
        
        ! Interpolate
        ratx = (ent - stab(i1,3)) / (stab(i2,3) - stab(i1,3))
        temp = stab(i1,1) + ratx * (stab(i2,1) - stab(i1,1))
        rho  = stab(i1,2) + ratx * (stab(i2,2) - stab(i1,2))
        prx  = stab(i1,4) + ratx * (stab(i2,4) - stab(i1,4))
        kvx  = stab(i1,5) + ratx * (stab(i2,5) - stab(i1,5))
        tcx  = stab(i1,6) + ratx * (stab(i2,6) - stab(i1,6))
        if (present(Rx)) then
          Rx = 1000._DP * (stab(i2,2) - stab(i1,2)) / (stab(i2,3) - stab(i1,3))
        end if
        
    end subroutine gettd

    !===============================================================================================!
    ! To calculate thermal conductivity of cladding
    !===============================================================================================!

    pure function getkc(t) result(res)
    
        real(dp), intent(in) :: t
        real(dp)             :: res
        
        res = 7.51_DP + 2.09e-2_DP*t - 1.45e-5_DP*t**2 + 7.67e-9_DP*t**3
    
    end function getkc
    
    !===============================================================================================!
    ! To calculate thermal conductivity of fuel
    !===============================================================================================!
    
    pure function getkf(t) result(res)
    
        real(dp), intent(in) :: t
        real(dp)             :: res
        
        res = 1.05_DP + 2150.0_DP / (t - 73.15_DP)
    
    end function getkf
    
    !===============================================================================================!
    ! To calculate specific heat capacity of cladding
    !===============================================================================================!
    
    pure function getcpc(t) result(res)
    
        real(dp), intent(in) :: t
        real(dp)             :: res
        
        res = 252.54_DP + 0.11474_DP*t
    
    end function getcpc
    
    !===============================================================================================!
    ! To calculate specific heat capacity of fuel
    !===============================================================================================!
    
    pure function getcpf(t) result(res)
    
        real(dp), intent(in) :: t
        real(dp)             :: res
        
        res = 162.3_DP + 0.3038_DP*t - 2.391e-4_DP*t**2 + 6.404e-8_DP*t**3
    
    end function getcpf

    !===============================================================================================!
    ! To calculate heat transfer coef.
    !===============================================================================================!

    real(dp) function geths(t, xden, tc, kv, Pr)
    
        class(th_type)        :: t
        real(dp), intent(in)  :: xden  ! coolant densisty
        real(dp), intent(in)  :: tc    ! coolant thermal conductivity
        real(dp), intent(in)  :: kv    ! kinematic viscosity
        real(dp), intent(in)  :: Pr    ! Prandtl Number
        
        real(dp) :: cvelo, Nu, Re
        
        cvelo = t % cflow / (t % s_area * xden * 1000._DP)    ! Calculate flow velocity (m/s)
        Re = cvelo * t % dh / (kv * 1.e-6_DP)                 ! Calculate Reynolds Number
        Nu = 0.023_DP*(Pr**0.4_DP)*(Re**0.8_DP)               ! Calculate Nusselt Number
        geths = (tc / t % dh) * Nu                            ! Calculate heat transfer coefficient
    
    end function geths

    !===============================================================================================!
    ! To update thermal parameters
    !===============================================================================================!

    subroutine th_update(t, linear_pow, tfm, heat_flux, enthalpy, ftem, mtem, cden)
    
        class(th_type)          :: t
        real(dp), intent(in)    :: linear_pow(:) ! Linear Power Density (W/cm)
        real(dp), intent(inout) :: heat_flux(:)  ! heat flux
        real(dp), intent(inout) :: tfm(:,:)      ! Fuel pin mesh temperature (nzz, n_ring)
        real(dp), intent(out)   :: enthalpy(:)   ! coolant enthalpy
        real(dp), intent(out)   :: ftem(:)       ! fuel temperature
        real(dp), intent(out)   :: mtem(:)       ! moderator temperature
        real(dp), intent(out)   :: cden(:)       ! moderator density
        
        integer   :: i, k
        real(dp)  :: a(t % n_ring + 1)
        real(dp)  :: b(t % n_ring + 1)
        real(dp)  :: c(t % n_ring + 1)
        real(dp)  :: d(t % n_ring + 1)
        real(dp)  :: hs, kt, kt1, kt2
        real(dp)  :: alp = 0.7_DP
        real(dp)  :: xa, xc
        real(dp)  :: pdens           ! power densisty  (W/m3)
        real(dp)  :: ent_inlet       ! Coolant inlet enthalpy
        real(dp)  :: cpline          ! Coolant Linear power densisty (W/m)
        real(dp)  :: Pr, kv, tcon    ! Coolant Prandtl Number, Kinematic viscosity, and thermal conductivity
        real(dp)  :: zd              ! delta z in meter
        real(dp)  :: ent_bound       ! enthalpy at the node boundary
        integer   :: error_flag

        
        !set initial tridiagonal matrix element a, b
        a = 0._dp; b = 0._dp; c = 0._dp;
        
        ! Get enthalpy or the lowest node
        call getent(t % t_inlet, ent_inlet)

        ! sweep in the axial direction
        do k = 1, t % nzz
            cpline = heat_flux(k) * pi * t % dia  + t % cf * linear_pow(k) * 100._DP     ! Coolant Linear power densisty (W/m)
            zd = t % zdel(k) * 0.01_DP
            if (k == 1) then
                enthalpy(k) = ent_inlet + 0.5_DP * cpline * zd / t % cflow               ! Calculate coolant enthalpy
                ent_bound   = 2. * enthalpy(k) - ent_inlet
            else
                enthalpy(k) = ent_bound + 0.5_DP * cpline * zd / t % cflow
                ent_bound   = 2. * enthalpy(k) - ent_bound
            end if
            call gettd(enthalpy(k), mtem(k), cden(k), Pr, kv, tcon)
        
            hs = geths(t, cden(k), Pr, kv, tcon)                                         ! calculate heat transfer coeff
            pdens = (1._DP - t % cf) * 100._DP * linear_pow(k) / (pi * t % rf**2)        ! Fuel pin Power Density (W/m3)
        
            ! Calculate tridiagonal matrix: a, b, c and source: d
            ! For nt=1 [FUEL CENTERLINE]
            kt1  = getkf(tfm(k,1))                                                        ! Get thermal conductivity
            kt2  = getkf(tfm(k,2))
            kt   = 2._DP * kt1 * kt2 / (kt1 + kt2)
            xc   = kt * t % rpos(1) / t % rdel(1)
            b(1) =  xc
            c(1) = -xc
            d(1) = pdens * 0.5_DP * t % rpos(1)**2
        
            do i = 2, t % n_ring-2
                kt1  = kt2
                kt2  = getkf(tfm(k, i+1))
                kt   = 2._DP * kt1 * kt2 / (kt1 + kt2)
                xa   = xc
                xc   = kt * t % rpos(i) / t % rdel(i)
                a(i) = -xa
                b(i) =  xa + xc
                c(i) = -xc
                d(i) = pdens * 0.5_DP * (t % rpos(i)**2 - t % rpos(i-1)**2)
            end do
          
            ! For nt-1 [FUEL-GAP INTERFACE]
            xa = xc
            xc = t % rg * Hg
            a(t % n_ring-1) = -xa
            b(t % n_ring-1) =  xa + xc
            c(t % n_ring-1) = -xc
            d(t % n_ring-1) = pdens * 0.5_DP * (t % rf**2 - t % rpos(t % n_ring-2)**2)
          
            ! For nt [GAP-CLADDING INTERFACE]
            kt1           = getkc(tfm(k, t % n_ring))
            kt2           = getkc(tfm(k, t % n_ring+1))
            kt            = 2._DP * kt1 * kt2 / (kt1 + kt2)     ! For cladding
            xa            = xc
            xc            = kt * t % rpos(t % n_ring) / t % rdel(t % n_ring)
            a(t % n_ring) = -xa
            b(t % n_ring) =  xa + xc
            c(t % n_ring) = -xc
            d(t % n_ring) = 0.
          
            ! For nt+1  [CLADDING-COOLANT INTERFACE]
            xa              = xc
            a(t % n_ring+1) = -xa
            b(t % n_ring+1) = xa + hs * t % rc
            d(t % n_ring+1) = t % rc * hs * mtem(k)
          
            ! Solve tridiagonal matrix
            call TridiaSolve(a, b, c, d, tfm(k, :), error_flag)
            if(error_flag .ne. 0) then
                write(error_unit,*) "Error: TridiaSolve routine failed"
                stop
            end if
          
            ! Get lumped fuel temp
            ftem(k) = (1. - t % alpha) * tfm(k, 1) + alp * tfm(k, t % n_ring-1)
          
            ! Calculate heat flux
            heat_flux(k) = hs * (tfm(k, t % n_ring+1) - mtem(k))
        end do
    
    end subroutine

    !===============================================================================================!
    ! To update thermal parameters (transient case)
    !===============================================================================================!

    subroutine th_transient(t, h, lin_power, flow_rate, enthalpy, &
        tfm, heat_flux, ftem, mtem, cden)
        
        class(th_type)          :: t
        real(dp), intent(in)    :: h             ! Time step
        real(dp), intent(in)    :: lin_power(:)  ! Linear Power Density (W/cm)
        real(dp), intent(inout) :: flow_rate(:)  ! coolant mass flow rate (kg/s)
        real(dp), intent(inout) :: enthalpy(:)   ! updated enthalpy
        real(dp), intent(inout) :: heat_flux(:)  ! heat flux
        real(dp), intent(inout) :: tfm(:,:)      ! Fuel pin mesh temperature (nzz, n_ring)
        real(dp), intent(out)   :: ftem(:)       ! fuel temperature
        real(dp), intent(out)   :: mtem(:)       ! moderator temperature
        real(dp), intent(out)   :: cden(:)       ! moderator density
        
        integer   :: i, k
        real(dp)  :: a(t % n_ring + 1)
        real(dp)  :: b(t % n_ring + 1)
        real(dp)  :: c(t % n_ring + 1)
        real(dp)  :: d(t % n_ring + 1)
        real(dp)  :: hs                                    ! coolant heat transfer coef
        real(dp)  :: kt , kt1, kt2                         ! thermal conductivity
        real(dp)  :: xa, xc
        real(dp)  :: cp                                    ! Specific heat capacity
        real(dp)  :: eps, eta
        real(dp)  :: mdens, vol                            ! Coolant density and channel volume
        real(dp)  :: ent_prev(t % nzz)                     ! previous enthalpy
        real(dp)  :: pdens                                 ! power densisty  (W/m3)
        real(dp)  :: ent_inlet                             ! Coolant inlet enthalpy
        real(dp)  :: cpline                                ! Coolant Linear power densisty (W/m)
        real(dp)  :: Pr, kv, tcon                          ! Coolant Prandtl Number, Kinematic viscosity, and thermal conductivity
        real(dp)  :: R
        real(dp)  :: sq_rdel
        real(dp)  :: ent_bound                             ! enthalpy at the node boundary
        integer   :: error_flag
        
        !set initial tridiagonal matrix element a, b
        a = 0._dp; b = 0._dp; c = 0._dp;

        ! Get enthalpy or the lowest node
        call getent(t % t_inlet, ent_inlet)
        ent_prev = enthalpy
        
        do k = 1, t % nzz
            mdens  = cden(k) * 1000._DP                                                 ! Coolant density (kg/m3)
            cpline = heat_flux(k) * pi * t % dia + t % cf * lin_power(k) * 100._DP      ! Coolant Linear power densisty (W/m)
            vol    = t % s_area * t % zdel(k) * 0.01_DP                                 ! channel node volume
            eps    = mdens * vol / h
        
            ! Calculate coolant enthalpy
            if (k == 1) then
                enthalpy(k)   = (cpline * t % zdel(k) * 0.01_DP + 2._DP * flow_rate(k) &
                * ent_inlet + eps * ent_prev(k)) / (eps + 2._DP * flow_rate(k))          ! Calculate enthalpy
                ent_bound     = 2. * enthalpy(k) - ent_inlet
                flow_rate(k)  = t % cflow - 0.5_dp * vol / h * R * (enthalpy(k) - ent_prev(k))
            else
                enthalpy(k)  = (cpline * t % zdel(k) * 0.01_DP + 2._DP * flow_rate(k) &
                * ent_bound + eps * ent_prev(k)) / (eps + 2._DP * flow_rate(k))
                ent_bound    = 2. * enthalpy(k) - ent_bound
                flow_rate(k) = flow_rate(k-1) - 0.5_dp * vol / h * R * (enthalpy(k) - ent_prev(k))
            end if

            call gettd(enthalpy(k), mtem(k), cden(k), Pr, kv, tcon, R)                   ! Get corresponding temp and density
        
            hs = geths(t, cden(k), Pr, kv, tcon)                                         ! Calculate heat transfer coef
            pdens = (1._DP - t % cf) * 100._DP * lin_power(k) / (pi * t % rf**2)         ! Fuel pin Power Density (W/m3)
        
            ! Calculate tridiagonal matrix: a, b, c and source: d
            ! For nt=1 [FUEL CENTERLINE]
            kt1  = getkf(tfm(k,1))                                                     ! Get thermal conductivity
            kt2  = getkf(tfm(k,2))
            kt   = 2._DP * kt1 * kt2 / (kt1 + kt2)
            cp   = getcpf(tfm(k,1))                                                     ! Get specific heat capacity
            eta  = fdens * cp * t % rpos(1)**2 / (2._DP * h)
            xc   = kt * t % rpos(1) / t % rdel(1)
            b(1) = xc + eta
            c(1) = -xc
            d(1) = pdens * 0.5_DP * t % rpos(1)**2 + eta * tfm(k,1)
        
            do i = 2, t % n_ring-2
                kt1     = kt2
                kt2     = getkf(tfm(k,i+1))
                kt      = 2._DP * kt1 * kt2 / (kt1 + kt2)
                cp      = getcpf(tfm(k,i))
                sq_rdel = t % rpos(i)**2 - t % rpos(i-1)**2
                eta     = fdens * cp * sq_rdel / (2. * h)
                xa      = xc
                xc      = kt * t % rpos(i) / t % rdel(i)
                a(i)    = -xa
                b(i)    =  xa + xc + eta
                c(i)    = -xc
                d(i)    = pdens * 0.5_DP * sq_rdel + eta * tfm(k,i)
            end do
        
            ! For nt-1 [FUEL-GAP INTERFACE]
            cp              = getcpf(tfm(k,t % n_ring-1))
            sq_rdel         = t % rf**2 - t % rpos(t % n_ring-2)**2
            eta             = fdens * cp * sq_rdel / (2. * h)
            xa              = xc
            xc              = t % rg * hg
            a(t % n_ring-1) = -xa
            b(t % n_ring-1) =  xa + xc + eta
            c(t % n_ring-1) = -xc
            d(t % n_ring-1) = pdens * 0.5_DP * sq_rdel + eta * tfm(k,t % n_ring-1)
        
            ! For nt [GAP-CLADDING INTERFACE]
            kt1           = getkc(tfm(k,t % n_ring))
            kt2           = getkc(tfm(k,t % n_ring+1))
            kt            = 2._DP * kt1 * kt2 / (kt1 + kt2)     ! For cladding
            cp            = getcpc(tfm(k,t % n_ring))
            eta           = cdens * cp * (t % rpos(t % n_ring)**2 - t % rg**2) / (2. * h)
            xa            = xc
            xc            = kt * t % rpos(t % n_ring) / t % rdel(t % n_ring)
            a(t % n_ring) = -xa
            b(t % n_ring) =  xa + xc + eta
            c(t % n_ring) = -xc
            d(t % n_ring) = eta * tfm(k,t % n_ring)
        
            ! For nt+1  [CLADDING-COOLANT INTERFACE]
            cp              = getcpc(tfm(k,t % n_ring+1))
            eta             = cdens * cp * (t % rc**2 - t % rpos(t % n_ring)**2) / (2. * h)
            xa              = xc
            xc              = t % rc * hs
            a(t % n_ring+1) = -xa
            b(t % n_ring+1) =  xa + xc + eta
            d(t % n_ring+1) = t % rc * hs * mtem(k) + eta * tfm(k,t % n_ring+1)
        
            ! Solve tridiagonal matrix
            call TridiaSolve(a, b, c, d, tfm(k, :), error_flag)
            if(error_flag .ne. 0) then
                write(error_unit,*) "Error: TridiaSolve routine failed"
                stop
            end if
        
            ! Get lumped fuel temp
            ftem(k) = (1.- t % alpha) * tfm(k, 1) + t % alpha * tfm(k, t % n_ring-1)
        
            ! Calculate heat flux
            heat_flux(k) = hs * (tfm(k, t % n_ring+1) - mtem(k))
        end do
        
    end subroutine
    
end module
    