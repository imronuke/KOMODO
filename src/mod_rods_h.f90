module rods_h

    use constant

    implicit none

    save

    ! Crod changes
    integer              :: nbank                  ! Number of CR banks
    real(dp)             :: nstep                  ! Number of steps
    real(dp)             :: coreh                  ! Core Height
    integer, allocatable :: fbmap(:,:)             ! Radial control rod bank map (node wise)
    real(dp)             :: zero_pos               ! Zero step position
    real(dp)             :: step_size              ! Step size
    
    type :: cr_bank_type
        real(dp), allocatable     :: bpos  (:)    ! CR bank position
        real(dp), allocatable     :: fbpos (:)    ! Final CR bank position
        real(dp), allocatable     :: tmove (:)    ! Time when CR bank starts moving
        real(dp), allocatable     :: bspeed(:)    ! CR bank movement speed
        integer, allocatable      :: mdir  (:)    ! To indicate CR movement direction (0=do not move, 1=down, 2 = up)
    end type


    
end module
    