module time

    use iso_fortran_env, only: real64, int64

    implicit none

    save

    integer, parameter :: dp = real64
    integer, parameter :: i8 = int64

    type :: timer
    ! integer(int8)  :: start_time
    real(8)  :: start_time
    real(8)  :: stop_time
    real(dp)    :: elapsed_time = 0.0
    contains
    procedure :: on => start_timer
    procedure :: off => stop_timer
    end type

    contains

    subroutine start_timer(this)

        class(timer) :: this

        ! call system_clock(this % start_time)
        call cpu_time(this % start_time)

    end subroutine

    subroutine stop_timer(this)

        class(timer) :: this
        ! integer(int8) :: end_counts, count_rate

        ! call system_clock(end_counts, count_rate)
        call cpu_time(this % stop_time)
        ! this % elapsed_time = this % elapsed_time + real(end_counts - this % start_time, dp)/real(count_rate, dp)
        this % elapsed_time = this % elapsed_time + (this % stop_time - this % start_time)

    end subroutine

end module