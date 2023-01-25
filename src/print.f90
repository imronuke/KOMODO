module print
    use iso_fortran_env, only: output_unit
    use data
    use utilities

    implicit none
  
    save

    public

    contains

    !===============================================================================================!
    ! Print header
    !===============================================================================================!

    SUBROUTINE print_head()
    
        IMPLICIT NONE
        if (mode == 'FORWARD' .and. bther == 1) then
            WRITE(ounit,*); WRITE(ounit,*)
            WRITE(ounit,3245); WRITE(ounit,3247) mode; WRITE(ounit,3245)
            WRITE(ounit,*)
            WRITE(ounit,3251); WRITE(ounit,1179)
            write(output_unit,*); write(output_unit,*)
            write(output_unit,3245); write(output_unit,3247) mode; write(output_unit,3245)
            write(output_unit,*)
            write(output_unit,3251); write(output_unit,1179)
        else if (mode == 'FORWARD' .or. mode == 'ADJOINT') then
            WRITE(ounit,*); WRITE(ounit,*)
            WRITE(ounit,3245); WRITE(ounit,3247) mode; WRITE(ounit,3245)
            WRITE(ounit,*)
            WRITE(ounit,3248); WRITE(ounit,3249)
            write(output_unit,*); write(output_unit,*)
            write(output_unit,3245); write(output_unit,3247) mode; write(output_unit,3245)
            write(output_unit,*)
            write(output_unit,3248); write(output_unit,3249)
        else if (mode == 'FIXEDSRC') then
            WRITE(ounit,*); WRITE(ounit,*)
            WRITE(ounit,3245); WRITE(ounit,3247) mode; WRITE(ounit,3245)
            WRITE(ounit,*)
            WRITE(ounit,3248); WRITE(ounit,3249)
            write(output_unit,*); write(output_unit,*)
            write(output_unit,3245); write(output_unit,3247) mode; write(output_unit,3245)
            write(output_unit,*)
            write(output_unit,3250); write(output_unit,3249)
        else
            if (bther == 0) then
                ! File Output
                WRITE(ounit,*); WRITE(ounit,*)
                WRITE(ounit,2176); WRITE(ounit,2177); WRITE(ounit,2176)
                WRITE(ounit,*); WRITE(ounit,2178); WRITE(ounit,2179)
                ! Terminal Output
                write(output_unit,*); write(output_unit,*)
                write(output_unit,2176); write(output_unit,2177); write(output_unit,2176)
                write(output_unit,*); write(output_unit,2178); write(output_unit,2179)
            else
                ! File Output
                WRITE(ounit,*); WRITE(ounit,*)
                WRITE(ounit,1176); WRITE(ounit,1177); WRITE(ounit,1176)
                WRITE(ounit,*); WRITE(ounit,1178); WRITE(ounit,1179)
                ! Terminal Output
                write(output_unit,*); write(output_unit,*)
                write(output_unit,1176); write(output_unit,1177); write(output_unit,1176)
                write(output_unit,*); write(output_unit,1178); write(output_unit,1179)
            end if
        end if
    
        3245 format (1X, ' ==============================================', &
        '================================')
        3247 format(23X, A8, ' CALCULATION RESULTS')
        3248 format(2X,'Itr     k-eff     Fis.Src Error     Flux error')
        3249 format(1X,'----------------------------------------------------')
        3250 format(2X,'Itr   Fis.Src Error     Flux error')
        2176 format(' ============================================================')
        2177 format(12X,'CRITICAL BORON CONCENTRATION SEARCH')
        2178 format('Itr  Boron Conc.   K-eff     Flux Error   Fiss. Src. Error')
        2179 format(' -----------------------------------------------------------')
        1176 format  &
        (' =========================================================================')
        1177 format(19X,'CRITICAL BORON CONCENTRATION SEARCH')
        1178 format &
        ('Itr  Boron Conc.   K-eff     Flux Err.    Fiss. Src. Err.  Fuel Temp. Error.')
        1179 format &
        (' -----------------------------------------------------------------------')
        3251 format(2X,'Itr     k-eff     Fis.Src Error     Flux error    Fuel Temp. Error')
    
    END SUBROUTINE print_head

    !===============================================================================================!
    ! To print final th parameters
    !===============================================================================================!

    SUBROUTINE print_tail()

        REAL(DP) :: tf, tm, mtm, mtf, otm, cd, ocd
    
        CALL par_ave(ftem, tf)
        CALL par_ave(mtem, tm)
    
        mtf = maxval(tfm(:,1))
        mtm = maxval(mtem)
    
        CALL par_ave_out(mtem, otm)
        CALL par_ave(cden, cd)
        CALL par_ave_out(cden, ocd)
    
        ! Write Output
        WRITE(ounit,*)
        WRITE(ounit, 5001) tf, tf-273.15; WRITE(ounit, 5002)  mtf, mtf-273.15
        WRITE(ounit, 5003) tm, tm-273.15; WRITE(ounit, 5004) mtm, mtm-273.15
        WRITE(ounit, 5005) otm, otm-273.15; WRITE(ounit, 5006) cd * 1000., cd
        WRITE(ounit, 5007) ocd * 1000., ocd
        write(output_unit,*)
        write(output_unit, 5001) tf, tf-273.15; write(output_unit, 5002)  mtf, mtf-273.15
        write(output_unit, 5003) tm, tm-273.15; write(output_unit, 5004) mtm, mtm-273.15
        write(output_unit, 5005) otm, otm-273.15; write(output_unit, 5006) cd * 1000., cd
        write(output_unit, 5007) ocd * 1000., ocd
    
        5001 FORMAT(2X, 'AVERAGE FUEL TEMPERATURE        : ', F7.1, ' K (', F7.1, ' C)')
        5002 FORMAT(2X, 'MAX FUEL CENTERLINE TEMPERATURE : ', F7.1, ' K (', F7.1, ' C)')
        5003 FORMAT(2X, 'AVERAGE MODERATOR TEMPERATURE   : ', F7.1, ' K (', F7.1, ' C)')
        5004 FORMAT(2X, 'MAXIMUM MODERATOR TEMPERATURE   : ', F7.1, ' K (', F7.1, ' C)')
        5005 FORMAT(2X, 'OUTLET MODERATOR TEMPERATURE    : ', F7.1, ' K (', F7.1, ' C)')
        5006 FORMAT(2X, 'AVERAGE MODERATOR DENSITY       : ', F7.1, ' kg/m3 (', F7.3, ' g/cc)')
        5007 FORMAT(2X, 'OUTLET MODERATOR DENSITY        : ', F7.1, ' kg/m3 (', F7.3, ' g/cc)')
    
    END SUBROUTINE print_tail

    !===============================================================================================!
    ! To calculate average fuel temp (only for active core)
    !===============================================================================================!

    SUBROUTINE par_ave(par, ave)
        
        REAL(DP), DIMENSION(:), INTENT(IN) :: par
        REAL(DP), INTENT(OUT) :: ave
        REAL(DP) :: dum, dum2
        INTEGER :: n
        
        dum = 0.; dum2 = 0.
        DO n = 1, nnod
           IF (fdm % xs % nuf(n,ng) > 0.) THEN
              dum = dum + par(n) * vdel(n)
              dum2 = dum2 + vdel(n)
           END IF
        END DO
        
        ave = dum / dum2
        
    END SUBROUTINE par_ave
        
    !===============================================================================================!
    ! To calculate average outlet coolant temperature
    !===============================================================================================!
        
    SUBROUTINE par_ave_out(par, ave)
        
        REAL(DP), DIMENSION(:), INTENT(IN) :: par
        REAL(DP), INTENT(OUT) :: ave
        REAL(DP) :: dum, dum2
        INTEGER, DIMENSION(nxx,nyy) :: zmax
        INTEGER :: n, i, j, k
        
        ! get number of nodex in axial direction from bottom -> fuel
        DO j = 1, nyy
          DO i = ystag(j)%smin, ystag(j)%smax
            zmax(i,j) = 0
            DO k = 1, nzz/2
              if (fdm % xs % nuf(xyz(i,j,k),ng) < 1.e-5_dp) zmax(i,j) = zmax(i,j) + 1
            END DO
          END DO
        END DO
        ! get number of nodex in axial direction from fuel -> top reflectors
        DO j = 1, nyy
          DO i = ystag(j)%smin, ystag(j)%smax
            DO k = 1, nzz
              if (fdm % xs % nuf(xyz(i,j,k),ng) > 1.e-5) zmax(i,j) = zmax(i,j) + 1
            END DO
          END DO
        END DO
        
        dum = 0.; dum2 = 0.
        DO n = 1, nnod
           IF (iz(n) == zmax(ix(n),iy(n)) .AND. fdm % xs % nuf(n,ng) > 1.e-5) THEN
              dum = dum + par(n) * vdel(n)
              dum2 = dum2 + vdel(n)
           END IF
        END DO
        
        ave = dum / dum2
        
    END SUBROUTINE par_ave_out

    !===============================================================================================!
    ! To print axially averaged assembly-wise power distribution
    !===============================================================================================!

    subroutine print_power_map(fn)
    
        real(dp), intent(in) :: fn(:)
        
        real(dp) :: fx(nxx, nyy, nzz)
        integer  :: i, j, k, n
        integer  :: ly, lx, ys, xs, yf, xf
        real(dp) :: summ, vsumm
        real(dp) :: fnode(nxx, nyy)
        real(dp) :: fasm(nx, ny)
        real(dp) :: totp
        integer  :: nfuel
        real(dp) :: fmax
        integer  :: xmax, ymax
        character(LEN=6) :: cpow(nx, ny)
        
        integer, parameter :: xm = 12
        integer :: ip, ipr

        nnod = size(fn)
        
        fx = 0._DP
        do n = 1, nnod
            fx(ix(n), iy(n), iz(n)) = fn(n)
        end do
        
        !Calculate axially averaged node-wise distribution
        fnode = 0._DP
        do j = 1, nyy
            do i = ystag(j)%smin, ystag(j)%smax
                summ = 0._DP
                vsumm = 0._DP
                do k = 1, nzz
                    summ = summ + fx(i,j,k)*zdel(k)
                    vsumm = vsumm + zdel(k)
                end do
                fnode(i,j)= summ/vsumm
            end do
        end do
        
        !Calculate assembly power
        nfuel = 0
        totp  = 0._DP
        ys = 1
        yf = 0
        do j= 1, ny
            yf = yf + ydiv(j)
            xf = 0
            xs = 1
            do i= 1, nx
                xf = xf + xdiv(i)
                summ = 0._DP
                vsumm = 0._DP
                do ly= ys, yf
                    do lx= xs, xf
                        summ = summ + fnode(lx,ly)*xdel(lx)*ydel(ly)
                        vsumm = vsumm + xdel(lx)*ydel(ly)
                    end do
                end do
                fasm(i,j) = summ / vsumm
                xs = xs + xdiv(i)
                if (fasm(i,j) > 0._DP) nfuel = nfuel + 1
                if (fasm(i,j) > 0._DP) totp  = totp + fasm(i,j)
            end do
            ys = ys + ydiv(j)
        end do
        
        
        ! Normalize assembly power to 1._DP
        xmax = 1; ymax = 1
        fmax = 0._DP
        do j = 1, ny
            do i = 1, nx
                if (totp > 0.) fasm(i,j) = real(nfuel) / totp * fasm(i,j)
                if (fasm(i,j) > fmax) then     ! Get max position
                    xmax = i
                    ymax = j
                    fmax = fasm(i,j)
                end if
                ! Convert power to character (if power == 0 convert to blank spaces)
                if ((fasm(i,j) - 0.) < 1.e-5_DP) then
                    cpow(i,j) = '     '
                else
                    write (cpow(i,j),'(F6.3)') fasm(i,j)
                    cpow(i,j) = trim(cpow(i,j))
                end if
            end do
        end do
        
        
        ! Print assembly power distribution
        write(ounit,*)
        write(ounit,*)
        write(ounit,*) '    Radial Power Distribution'
        write(ounit,*) '  =============================='
        
        ip = nx/xm
        ipr = MOD(nx,xm) - 1
        xs = 1; xf = xm
        do k = 1, ip
            write(ounit,'(4X,100I8)') (i, i = xs, xf)
            do j= ny, 1, -1
                write(ounit,'(2X,I4,100A8)') j, (cpow(i,j), i=xs, xf)
            end do
            write(ounit,*)
            xs = xs + xm
            xf = xf + xm
        end do
        
        
        write(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
        if (xs+ipr > xs) then
            do j= ny, 1, -1
                write(ounit,'(2X,I4,100A8)') j, (cpow(i,j), i=xs, xs+ipr)
            end do
        end if
        
        
        
        write(ounit,*)
        
        write(ounit,*) '  MAX POS.       Maximum Value'
        write(ounit,1101) ymax, xmax, fasm(xmax, ymax)
        
        1101 format(2X, '(' , I3, ',', I3,')', F15.3)
    
    
    end subroutine print_power_map

    !===============================================================================================!
    ! calculate power distribution (normalize to number of nodes)
    !===============================================================================================!

    subroutine power_dist_print(f0, pdist)

        real(dp), intent(in) :: f0(:,:)  ! Flux
        real(dp), intent(out) :: pdist(:)

        integer :: g, n
        real(dp) :: pow, vtot, powtot
        
        pdist = 0._dp
        do g= 1, ng
            do n= 1, nnod
              pow = f0(n,g) * fdm % xs % sigf(n,g)
              if (pow < 0.) pow = 0.
              pdist(n) = pdist(n) + pow
            end do
        end do
        
        ! Normalize the nnod
        powtot = 0._dp
        vtot = 0.0
        do n = 1, nnod
            powtot = powtot + pdist(n) * vdel(n)
            vtot   = vtot + vdel(n)
        end do
        
        if (powtot <= 0 .and. (mode .ne. 'FIXEDSRC')) THEN
            call fatal_error(ounit, '   ERROR: TOTAL NODES POWER IS ZERO OR LESS' // &
            new_line('a') // '   STOP IN subroutine POWDIS')
        end if
        
        ! sum of pdist must be equal to nnod
        do n = 1, nnod
            pdist(n) = pdist(n) * vtot / powtot
        end do

    end subroutine

end module