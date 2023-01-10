module stagger

    use iso_fortran_env, only: real64

    implicit none

    private

    save

    integer, parameter :: dp = real64

    ! Staggered mesh
    type, public :: staggered
        integer :: smax, smin                      ! imax and imin along x and y direction for staggered nodes
    end type

end module
    