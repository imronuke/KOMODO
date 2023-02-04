module transient_aux

    ! Module auxilaries for the transient module

    use iso_fortran_env, only: real64, error_unit, output_unit
    use data
    use nodal
    use xsec
    use fdm
    use print, only: print_fail_converge

    implicit none

    private

    save

    integer, parameter  :: dp = real64

    public :: adjust_keff, get_leakage, integrate, xs_update, move_crod, get_adjoint_flux

    contains

    !===============================================================================================!
    ! To adjuts the Keff to 1.0 if it is not equal to 1.0
    ! When starting transient, normally reactor at steady state condition (k = 1.0).
    ! In case the k is not equal to 1.0, then it is adjusted to 1.0
    !===============================================================================================!

    subroutine adjust_keff(c, xsc)

        class(fdm_rect)       :: c
        class(xs_change_type) :: xsc
        
        integer :: i
        logical :: converge
        
        write(output_unit,'(A)', advance='no') '  steady state k-eff not equal to one, force it to one ... '
        write(ounit, *)
        write(ounit, '(A46,F9.6)') '  INITIAL MULTIPLICATION EFFECTIVE (K-EFF) = ',  c % Keff
        write(ounit, *) '  WARNING: THE STEADY STATE K-EFF IS NOT EQUAL TO 1.0'
        write(ounit, *) '  AND NOW IT IS FORCED TO 1.0 BY MODIFYING THE nu*sigf CROSS SECTIONS '
        write(ounit, *)
        do i = 1, 10
            ! Adjust material-wise nu*fission and xs changes due to CR
            nuf = nuf / c % Keff
            if (allocated(xsc % crod)) xsc % crod % nuf = xsc % crod % nuf / c % Keff

            call xsec_update(c % xs, xsc, bcon, ftem, mtem, cden, bank_pos)
            call outer_iter(c, print_iter = .false., is_converge=converge)
            if (.not. converge) call print_fail_converge()
            if (abs(c % Keff-1.0) < 1.e-5) exit
        end do
        if (i == 10) then
            write(error_unit,*) "K-EFF STILL NOT EQUAL TO ONE. KOMODO IS STOPPING"
            stop
        end if

        write(output_unit,*) 'done'
        
    end subroutine

    !===============================================================================================!
    ! To calculate adjoint flux
    !===============================================================================================!

    subroutine get_adjoint_flux(xsc, adj_flux)
        
        class(xs_change_type), intent(in) :: xsc
        real(dp), intent(out)             :: adj_flux(:,:)

        type(fdm_rect), allocatable :: fdm_adj
        real(dp)                     :: flux_tmp(nnod, ng), fsrc_tmp(nnod)
        logical                      :: converge
        integer                      :: nodal_interval

        allocate(fdm_adj)

        nodal_interval = ceiling((nxx + nyy + nzz) / 2.5)
        call set_fdm_data(fdm_adj, ounit, kern, ng, nnod, nxx, nyy, nzz, nmat, &
        east, west, north, south, top, bottom, 40, 1.e-5_dp, 1.e-5_dp, 1000, 2, 10)

        call set_fdm_pointer(fdm_adj, flux_tmp, fsrc_tmp, exsrc, ind_mat, xdel, ydel, zdel, vdel, &
        xstag, ystag, ix, iy, iz, xyz)

        call xsec_setup(fdm_adj % xs, sigtr, siga, nuf, sigf, sigs, chi, dc)
        call xsec_update(fdm_adj % xs, xsc, bcon, ftem, mtem, cden, bank_pos)

        call outer_adjoint(fdm_adj, print_iter = .false., is_converge=converge)
        if (.not. converge) call print_fail_converge()

        adj_flux = flux_tmp

        call deallocate_fdm(fdm_adj)
 
    end subroutine

    !===============================================================================================!
    ! To calculate leakage in x, y and z directions
    !===============================================================================================!

    subroutine get_leakage(c, Lx, Ly, Lz)
        
        class(fdm_rect)                      :: c
        real(dp), dimension(:,:), intent(out) :: Lx, Ly, Lz

        call TLUpd0(c % d, Lx, Ly, Lz)
 
    end subroutine

    !===============================================================================================!
    ! To update xs
    !===============================================================================================!

    subroutine xs_update(xs, xsc)

        type(xs_rect)         :: xs
        class(xs_change_type) :: xsc
        
        call xsec_update(xs, xsc, bcon, ftem, mtem, cden, bank_pos)
        
    end subroutine

    !===============================================================================================!
    ! To update xs
    !===============================================================================================!

    subroutine move_crod(t_step, t_next)

        real(dp), intent(in) :: t_step, t_next

        integer :: n

        do n = 1, nbank
            if (direction(n) == 1) then   ! If CR moving down
                if (t_next-t_move(n) > 1.e-5_dp .AND. bank_pos_final(n)-bank_pos(n) < 1.e-5_dp) then
                    bank_pos(n) = bank_pos(n) - t_step*bank_speed(n)
                    if (bank_pos(n) < bank_pos_final(n)) bank_pos(n) = bank_pos_final(n)  ! If bpos exceed, set to fbpos
                end if
            else if (direction(n) == 2) then ! If CR moving up
                if (t_next-t_move(n) > 1.e-5_dp .AND. bank_pos_final(n)-bank_pos(n) > 1.e-5_dp) then
                    bank_pos(n) = bank_pos(n) + t_step*bank_speed(n)
                    if (bank_pos(n) > bank_pos_final(n)) bank_pos(n) = bank_pos_final(n)  ! If bpos exceed, set to fbpos
                end if
            end if
         end do
        
    end subroutine

    !===============================================================================================!
    ! do volume integration on variable a
    !===============================================================================================!

    pure function integrate(a) result(res)

        real(dp), intent(in)          :: a(:)
        real(dp)                      :: res

        integer                 :: n

        res = 0.
        do n = 1, nnod
            res = res + a(n) * vdel(n)
        end do
    
    end function

end module