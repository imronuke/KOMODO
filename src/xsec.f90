module xsec

    use iso_fortran_env, only: real64

    implicit none

    private

    save

    integer, parameter :: dp = real64
    
    ! Node-wise xs base
    type :: xs_base
        real(dp), allocatable :: D(:,:)                ! Diffusion coefficient
        real(dp), allocatable :: sigtr(:,:)          ! Transport macroscopic XSEC
        real(dp), allocatable :: siga (:,:)          ! Absorption macroscopic XSEC
        real(dp), allocatable :: nuf  (:,:)          ! nu*fission macroscopic XSEC
        real(dp), allocatable :: sigf (:,:)          ! fission macroscopic XSEC
        real(dp), allocatable :: sigr (:,:)          ! Removal macroscopic XSEC
        real(dp), allocatable :: sigs (:,:,:)        ! Scattering macroscopic XSEC
        real(dp), allocatable :: chi(:,:)            ! MG fission spectrum
    end type

    type, extends(xs_base), public :: xs_rect
        real(dp), allocatable :: dc(:,:,:)
    end type
    
    ! Change of XS per change of other parameters
    type :: xs_changes
        real(dp), allocatable :: sigtr(:)          ! Transport macroscopic XSEC
        real(dp), allocatable :: siga (:)          ! Absorption macroscopic XSEC
        real(dp), allocatable :: nuf  (:)          ! nu*fission macroscopic XSEC
        real(dp), allocatable :: sigf (:)          ! fission macroscopic XSEC
    end type

    type, extends(xs_changes) :: xs_changes_tab  ! if xtab file is used
        real(dp), allocatable :: adf (:,:)
    end type

    type :: xs_change_type
        type(xs_changes)  :: fuel_temp
        type(xs_changes)  :: mod_temp
        type(xs_changes)  :: mod_dens
        type(xs_changes)  :: rods
        type(xs_changes)  :: boron
    end type

    type :: rod_xs_change_type
        type(xs_changes_tab)  :: rod
    end type

    public :: alloc_xsec, set_xsec

    contains

    !===============================================================================================!
    ! Allocate MG XSEC (node-wise)                                                                  !
    !===============================================================================================!

    subroutine alloc_xsec(xs, nnod, ng, n_surf)

        class(xs_base), intent(inout)       :: xs
        integer, intent(in)                :: ng, nnod
        integer, intent(in), optional      :: n_surf

        allocate(xs % sigtr(nnod, ng), xs % siga(nnod, ng), xs % nuf(nnod, ng))
        allocate(xs % sigf(nnod, ng), xs % sigs(nnod, ng, ng), xs % chi(nnod, ng))
        allocate(xs % D(nnod, ng), xs % sigr(nnod, ng))

        select type (xs)
        type is (xs_rect)
            if (present(n_surf)) then
                allocate(xs % dc(nnod, ng, n_surf))
            else
                stop "ERROR: number of surgface is necessary in alloc_xsec routine"
            endif
        end select

    end subroutine

    !===============================================================================================!
    ! Set MG XSEC (node-wise)                                                                       !
    !===============================================================================================!

    subroutine set_xsec(xs, nnod, ix, iy, iz, mat_map, sigtr, siga, nuf, sigf, sigs, chi, dc)

        class(xs_base), intent(inout)   :: xs
        integer, intent(in)            :: nnod
        integer, intent(in)            :: ix(:), iy(:), iz(:)
        integer, intent(in)            :: mat_map(:,:,:)
        real(dp), intent(in)           :: sigtr(:,:)          ! Transport macroscopic XSEC
        real(dp), intent(in)           :: siga (:,:)          ! Absorption macroscopic XSEC
        real(dp), intent(in)           :: nuf  (:,:)          ! nu*fission macroscopic XSEC
        real(dp), intent(in)           :: sigf (:,:)          ! fission macroscopic XSEC
        real(dp), intent(in)           :: sigs (:,:,:)        ! Scattering macroscopic XSEC
        real(dp), intent(in)           :: chi(:,:)            ! MG fission spectrum
        real(dp), intent(in), optional :: dc(:,:,:)           ! ADF

        integer   :: n, i, j, k, nn

        ! Allocate node wise xs
        do n = 1, nnod
            i = ix(n); j = iy(n); k = iz(n)
            nn = mat_map(i, j, k)

            xs % sigtr(n,:)      = sigtr(nn,:)
            xs % siga (n,:)      = siga (nn,:)
            xs % nuf  (n,:)      = nuf  (nn,:)
            xs % sigf (n,:)      = sigf (nn,:)
            xs % sigs (n,:,:)    = sigs (nn,:,:)
            xs % chi (n,:)       = chi  (nn,:)

            select type (xs)
            type is (xs_rect)
                if (present(dc)) then
                    xs % dc(n, :, :) = dc(nn, :, :)
                else
                    xs % dc(n, :, :) = 1.
                endif
            end select
            
        end do
        
        call upd_diff_sigr(xs)

    end subroutine

    !===============================================================================================!
    ! ! Calculate diff coeff and removal XS                                                         !
    !===============================================================================================!

    pure subroutine upd_diff_sigr(xs)

        type(xs_base), intent(inout) :: xs

        integer :: nnod, ng
        integer :: g, h, n

        xs % sigr = xs % siga

        nnod = size(xs % siga, dim=1)
        ng   = size(xs % siga, dim=2)

        do h = 1, ng
            do g = 1, ng
                do n = 1, nnod
                    if (g /= h) xs % sigr(n, g) = xs % sigr(n, g) + xs % sigs(n, g, h)
                end do
            end do
        end do

        do g = 1, ng
            do n = 1, nnod
                xs % D(n, g) = 1 / (3. * xs % sigtr(n, g))
            end do
        end do
        
    end subroutine

    
end module


    