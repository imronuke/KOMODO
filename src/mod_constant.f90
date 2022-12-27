module constant

    implicit none

    save
    
    integer, parameter :: dp = selected_real_kind(10, 15)

    real(dp), parameter   :: pi = acos(-1.0)
    integer , parameter   :: length_word = 30
    integer , parameter   :: length_line = 200
    integer , parameter   :: length_card = 4

    integer, parameter :: YES = 1
    integer, parameter :: NO = 0

    
end module constant
    