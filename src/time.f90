module time

    use iso_fortran_env, only: real64, int64

    implicit none

    private

    save

    integer, parameter :: dp = real64
    integer, parameter :: int8 = int64

    type, public :: timer
    integer(int8)  :: start_time
    ! real(dp)  :: start_time
    real(dp)  :: stop_time
    real(dp)    :: elapsed_time = 0.0
    contains
    procedure :: on => start_timer
    procedure :: off => stop_timer
    end type

    !Timing
    type(timer), public :: inp_time
    type(timer), public :: fdm_time
    type(timer), public :: nodal_time
    type(timer), public :: xs_time
    type(timer), public :: th_time

    contains

    subroutine start_timer(this)

        class(timer) :: this

        call system_clock(this % start_time)
        ! call cpu_time(this % start_time)

    end subroutine

    subroutine stop_timer(this)

        class(timer) :: this
        integer(int8) :: end_counts, count_rate

        call system_clock(end_counts, count_rate)
        ! call cpu_time(this % stop_time)
        this % elapsed_time = this % elapsed_time + real(end_counts - this % start_time, dp)/real(count_rate, dp)
        ! this % elapsed_time = this % elapsed_time + (this % stop_time - this % start_time)

    end subroutine

end module