module print

    use data

    implicit none
  
    save

    public

    contains

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

end module