module control

    use data
    use print
    use fdm
    use xsec
    use time

    implicit none

    save

    private :: calculation_init

    contains

    !===============================================================================================!
    ! forwrad calculation
    !===============================================================================================!

    subroutine forward()

        integer :: g, n
        real(dp) :: fn(nnod)

        call calculation_init()

        call outer_iter(fdm, 1.e-5_dp, 1.e-5_dp, max_outer=1000, max_inner=2, &
        extrp = 5, print_iter = .true.)

        fn = 0
        do g = 1, ng
            do n = 1, nnod
                fn(n) = fn(n) + fdm % xs % nuf(n,g) * flux(n,g)
            end do
        end do

        call print_power_map(fn)

    end subroutine

    !===============================================================================================!
    ! forwrad calculation
    !===============================================================================================!

    subroutine adjoint()

        integer :: g, n
        real(dp) :: fn(nnod)

        call calculation_init()

        call outer_iter(fdm, 1.e-5_dp, 1.e-5_dp, max_outer=1000, max_inner=2, &
        extrp = 5, print_iter = .true.)

        fn = 0
        do g = 1, ng
            do n = 1, nnod
                fn(n) = fn(n) + fdm % xs % nuf(n,g) * flux(n,g)
            end do
        end do

        call print_power_map(fn)

    end subroutine

    !===============================================================================================!
    ! Initialize calculation
    !===============================================================================================!

    subroutine calculation_init()

        call fdm_time % on

        allocate(flux(nnod, ng))
        allocate(fsrc(nnod))
        allocate(exsrc(nnod, ng))
        exsrc = 0.0

        if (.not. allocated(dc)) then
            allocate(dc(nnod, ng, 6))
            dc = 1.0
        end if

        call set_rect_data(fdm, kern, ng, nnod, nxx, nyy, nzz, nmat, &
        east, west, north, south, top, bottom, ounit)
        
        call set_rect_pointer(fdm, flux, fsrc, exsrc, mat_map, xdel, ydel, zdel, vdel, &
        xstag, ystag, ix, iy, iz, xyz)

        call fdm_time % on

        call xs_time % on

        call alloc_xsec(fdm % xs, nnod, ng)
        call set_xsec(fdm % xs, nnod, ix, iy, iz, mat_map, &
        sigtr, siga, nuf, sigf, sigs, chi, dc)

        call xs_time % off

    end subroutine

end module
    