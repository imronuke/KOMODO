module CMFD

  use sdata, only: dp
  implicit none
  save

contains

  !****************************************************************************!

  subroutine coup_coef()

    !Purpose: to calculate FDM nodal coupling coefficients

    use sdata, only: ng, nnod, ix, iy, iz, xyz, D, xdel, ydel, zdel, &
    ystag, xstag, nod, nzz, xeast, xwest, ysouth, ynorth, zbott, ztop

    implicit none

    real(dp) :: alb = 1.e30_dp
    real(dp) :: d1, d2
    integer  :: n, g, i, j, k
    logical  :: first = .true.

    if (first) then
      allocate(nod(nnod,ng))
      do g = 1, ng
          do n = 1, nnod
              nod(n,g)%dn = 0._dp    ! initial nodal coupling coefficients
          end do
      end do
      first = .false.
    end if

    ! Calculate FDM coupling coefficients
    do g = 1, ng
      do n = 1, nnod

        ! Set i, j, k
        i = ix(n); j = iy(n); k = iz(n)

        ! Set FDM coupling coefficients in x direction
        if (i == ystag(j)%smax) then
          if (xeast == 0) then
            nod(n,g)%df(1) = 2._dp * alb * D(n,g) / &
                         (2._dp * D(n,g) + alb * xdel(i))
          else if (xeast == 1) then
            nod(n,g)%df(1) =  D(n,g) / (2._dp * D(n,g) + 0.5_dp * xdel(i))
          else
            nod(n,g)%df(1) = 0._dp
          end if
        else
          d2 = D(xyz(i+1, j, k), g)
          nod(n,g)%df(1) = 2._dp * D(n,g) * d2 / &
                       (D(n,g) * xdel(i+1) +  d2 * xdel(i))
        end if

        if (i == ystag(j)%smin) then
          if (xwest == 0) then
            nod(n,g)%df(2) = 2._dp * alb * D(n,g) / &
                         (2._dp * D(n,g) + alb * xdel(i))
          else if (xwest == 1) then
            nod(n,g)%df(2) =  D(n,g) / (2._dp * D(n,g) + 0.5_dp * xdel(i))
          else
            nod(n,g)%df(2) = 0._dp
          end if
        else
          d1 = D(xyz(i-1, j, k), g)
          nod(n,g)%df(2) = 2._dp * D(n,g) * d1 / &
                       (D(n,g) * xdel(i-1) +  d1 * xdel(i))
        end if

        ! Set nodal coupling coefficients in y direction
        if (j == xstag(i)%smax) then
          if (ynorth == 0) then
            nod(n,g)%df(3) = 2._dp * alb * D(n,g) / &
                         (2._dp * D(n,g) + alb * ydel(j))
          else if (ynorth == 1) then
            nod(n,g)%df(3) =  D(n,g) / (2._dp * D(n,g) + 0.5_dp * ydel(j))
          else
            nod(n,g)%df(3) = 0._dp
          end if
        else
          d2 = D(xyz(i, j+1, k), g)
          nod(n,g)%df(3) = 2._dp * D(n,g) * d2 / &
                       (D(n,g) * ydel(j+1) +  d2 * ydel(j))
        end if

        if (j == xstag(i)%smin) then
          if (ysouth == 0) then
            nod(n,g)%df(4) = 2._dp * alb * D(n,g) / &
                         (2._dp * D(n,g) + alb * ydel(j))
          else if (ysouth == 1) then
            nod(n,g)%df(4) =  D(n,g) / (2._dp * D(n,g) + 0.5_dp * ydel(j))
          else
            nod(n,g)%df(4) = 0._dp
          end if
        else
          d1 = D(xyz(i, j-1, k), g)
          nod(n,g)%df(4) = 2._dp * D(n,g) * d1 / &
                       (D(n,g) * ydel(j-1) +  d1 * ydel(j))
        end if

        ! Set nodal coupling coefficients in z direction
        if (k == nzz) then
          if (ztop == 0) then
            nod(n,g)%df(5) = 2._dp * alb * D(n,g) / &
                         (2._dp * D(n,g) + alb * zdel(k))
          else if (ztop == 1) then
            nod(n,g)%df(5) =  D(n,g) / (2._dp * D(n,g) + 0.5_dp * zdel(k))
          else
            nod(n,g)%df(5) = 0._dp
          end if
        else
          d2 = D(xyz(i, j, k+1), g)
          nod(n,g)%df(5) = 2._dp * D(n,g) * d2 / &
                       (D(n,g) * zdel(k+1) +  d2 * zdel(k))
        end if

        if (k == 1) then
          if (zbott == 0) then
            nod(n,g)%df(6) = 2._dp * alb * D(n,g) / &
                         (2._dp * D(n,g) + alb * zdel(k))
          else if (zbott == 1) then
            nod(n,g)%df(6) =  D(n,g) / (2._dp * D(n,g) + 0.5_dp * zdel(k))
          else
            nod(n,g)%df(6) = 0._dp
          end if
        else
          d1 = D(xyz(i, j, k-1), g)
          nod(n,g)%df(6) = 2._dp * D(n,g) * d1 / &
                       (D(n,g) * zdel(k-1) +  d1 * zdel(k))
        end if
      end do
    end do

  end subroutine coup_coef

  !****************************************************************************!

  subroutine set_ind()

    !Purpose: to set indexes for CMFD matrix

    use sdata, only: ix, iy, iz, ystag, xstag, nnod, nxx, nyy, nzz, ind

    implicit none

    integer, dimension(nxx,nyy) :: nodp  !radial node position
    integer :: n, np, rec, i, j, k

    ! setup radial node position nodp
    nodp = 0
    rec = 0
    do j = 1, nyy
      do i = ystag(j)%smin, ystag(j)%smax
        rec = rec + 1
        nodp(i,j) = rec
      end do
    end do

    ! Calculate number of nodes for one planar
    np = rec

    do n = 1, nnod

      ! Set i, j, k
      i = ix(n); j = iy(n); k = iz(n)
      rec = 0

      ! Lower diagonal matrix element for z-direction
      if (k /= 1) then
        rec = rec + 1
        ind(n)%col(rec) = n - np
      end if

      ! Lower diagonal matrix element for y-direction
      if (j /= xstag(i)%smin) then
        rec = rec + 1
        ind(n)%col(rec) = n - (nodp(i,j) - nodp(i,j-1))
      end if

      ! Lower diagonal matrix element for x-direction
      if (i /= ystag(j)%smin) then
        rec = rec + 1
        ind(n)%col(rec) = n - 1
      end if

      ! Diagonal matrix elementss
      rec = rec + 1
      ind(n)%col(rec) = n

       ! Upper diagonal matrix element for x-direction
      if (i /= ystag(j)%smax) then
        rec = rec + 1
        ind(n)%col(rec) = n + 1
      end if

      ! Upper diagonal matrix element for y-direction
      if (j /= xstag(i)%smax) then
        rec = rec + 1
        ind(n)%col(rec) = n + (nodp(i,j+1) - nodp(i,j))
      end if

      ! Upper diagonal matrix element for z-direction
      if (k /= nzz) then
        rec = rec + 1
        ind(n)%col(rec) = n + np
      end if

    end do

  end subroutine set_ind

  !****************************************************************************!

  subroutine matrix_setup(opt)

    !Purpose: to setup sparse penta-diagonal matrix. Elements are indexed in
    ! two-dimensional vector ind and non-zero elements strored in A

    use sdata, only: nod, ix, iy, iz, xdel, ydel, zdel, &
    ystag, xstag, nnod, sigr, nzz, ng, A, ind

    implicit none

    integer, intent(in)         :: opt
    integer :: n, g, rec, i, j, k
    logical :: first = .true.

    ! Allocate FDM matrix for first time
    if (first) then
      allocate(A(nnod,ng))
      do n = 1, nnod
        do g = 1, ng
          allocate(A(n,g)%elmn(ind(n)%ncol))
        end do
      end do
      call set_ind()
      first = .false.
    end if

    ! If need to recalculate FDM coupling coefficients
    if (opt > 0) call coup_coef()


    ! Setup CMFD linear system
    do g = 1, ng
      do n = 1, nnod

        ! Set i, j, k
        i = ix(n); j = iy(n); k = iz(n)
        rec = 0

        ! Lower diagonal matrix element for z-direction
        if (k /= 1) then
          rec = rec + 1
          A(n,g)%elmn(rec) = -(nod(n,g)%df(6) - nod(n,g)%dn(6)) / zdel(k)
        end if

        ! Lower diagonal matrix element for y-direction
        if (j /= xstag(i)%smin) then
          rec = rec + 1
          A(n,g)%elmn(rec) = -(nod(n,g)%df(4) - nod(n,g)%dn(4)) / ydel(j)
        end if

        ! Lower diagonal matrix element for x-direction
        if (i /= ystag(j)%smin) then
          rec = rec + 1
          A(n,g)%elmn(rec) = -(nod(n,g)%df(2) - nod(n,g)%dn(2)) / xdel(i)
        end if

        ! Diagonal matrix elementss
        rec = rec + 1
        A(n,g)%elmn(rec) = (nod(n,g)%df(1) + nod(n,g)%df(2) - &
                           nod(n,g)%dn(1) + nod(n,g)%dn(2)) / xdel(i) + &
                           (nod(n,g)%df(3) + nod(n,g)%df(4) - &
                           nod(n,g)%dn(3) + nod(n,g)%dn(4)) / ydel(j) + &
                           (nod(n,g)%df(5) + nod(n,g)%df(6) - &
                           nod(n,g)%dn(5) + nod(n,g)%dn(6)) / zdel(k) + &
                           sigr(n,g)

         ! Upper diagonal matrix element for x-direction
        if (i /= ystag(j)%smax) then
          rec = rec + 1
          A(n,g)%elmn(rec) = -(nod(n,g)%df(1) + nod(n,g)%dn(1)) / xdel(i)
        end if

        ! Upper diagonal matrix element for y-direction
        if (j /= xstag(i)%smax) then
          rec = rec + 1
          A(n,g)%elmn(rec) = -(nod(n,g)%df(3) + nod(n,g)%dn(3)) / ydel(j)
        end if

        ! Upper diagonal matrix element for z-direction
        if (k /= nzz) then
          rec = rec + 1
          A(n,g)%elmn(rec) = -(nod(n,g)%df(5) + nod(n,g)%dn(5)) / zdel(k)
        end if

      end do
    end do

  end subroutine matrix_setup

  !****************************************************************************!

  subroutine fiss_extrp(popt,e1,e2,erro,errn,fs)

  !
  ! Purpose:
  !    To perform fission source extrapolation

  USE io,    ONLY: ounit, scr

  implicit none

  integer, intent(in)                   :: popt
  real(dp), intent(in)                  :: e1, e2
  real(dp), dimension(:), intent(in)    :: erro, errn
  real(dp), dimension(:), intent(inout) :: fs


  real(dp) :: domiR, mval

  domiR = e2 / e1
  mval = MAXVAL(ABS(erro))
  if (mval * mval < 0.0) domiR = -domiR
  fs = fs + domiR / (1._DP - domiR) * errn
  if (popt > 0) then
    write(ounit,*) '    ...FISSION SOURCE EXTRAPOLATED...'
    if (scr) write(*,*) '    ...FISSION SOURCE EXTRAPOLATED...'
  end if

  end subroutine fiss_extrp

  !****************************************************************************!

  subroutine nodal_upd(popt, nmode)

  !
  ! Purpose:
  !    To update nodal coupling coefficients

  USE sdata, ONLY: ndmax, im, jm ,km, kern, get_time, nod_time
  USE nodal, ONLY: nodal_update, nodal_update_pnm
  USE io,    ONLY: ounit, scr


  implicit none

  integer, intent(in) :: popt
  integer, intent(in) :: nmode      ! Nodal update mode
  real(dp) :: st, fn

  st = get_time()

  ndmax = 0._dp
  !Update nodal coupling coefficients
  if (kern == 'SANM') then
    call nodal_update(nmode)
  else
    call nodal_update_pnm(nmode)
  end if
  !Update CMFD matrix
  call matrix_setup(0)
  if (popt > 0) then
    write(ounit,*) '    .....NODAL COUPLING UPDATED..... '
    write(ounit,1145) ndmax, im, jm, km
    if (scr) then
      write(*,*) '    .....NODAL COUPLING UPDATED..... '
      write(*,1145) ndmax, im, jm, km
    end if
  end if

  fn = get_time()

  nod_time = nod_time + (fn-st)

  1145 FORMAT ('MAX. CHANGE IN NODAL COUPLING COEF.= ', ES12.5, &
               ' AT NODE I = ', I2, ', J = ', I2, ', K = ', I2)

end subroutine nodal_upd

!****************************************************************************!

subroutine print_keff(popt)

!
! Purpose:
!    To update nodal coupling coefficients

USE sdata, ONLY: Ke
USE io,    ONLY: ounit, scr

implicit none

integer, intent(in) :: popt

if (popt > 0) then
  write(ounit,*)
  write(ounit,1146) Ke
  if (scr) then
    write(*,*)
    write(*,1146) Ke
  end if
end if

1146 format(2X,'MULTIPLICATION EFFECTIVE (K-EFF) = ', F9.6)

end subroutine print_keff

  !****************************************************************************!

  subroutine outer(popt)

  !
  ! Purpose:
  !    To perform forward outer iteration


  USE sdata, ONLY: ng, nnod, nout, nin, serc, ferc, fer, ser, f0, nupd, &
                   Ke, nac, fs0, s0, ndmax, kern, get_time, fdm_time
  USE io,    ONLY: ounit, scr, bther
  USE nodal, ONLY: nodal_update

  implicit none

  integer, optional, intent(in) :: popt

  real(dp)                     :: Keo                !Old Multiplication factor (Keff)
  real(dp), dimension(nnod)    :: fs0c               !old fission source
  real(dp), dimension(nnod)    :: bs                 !total source
  real(dp), dimension(nnod,ng) :: f0c                !Old flux
  real(dp) :: f, fc                                  ! new and old integrated fission sources
  integer  :: p, g
  real(dp) :: e1, e2
  real(dp), dimension(nnod) :: errn, erro            ! current and past error vectors
  logical  :: first = .true.
  real(dp) :: st, fn

  st = get_time()

  !Setup CMFD matrix
  CALL matrix_setup(1)

  !Allocate flux and fission source for first time
  if (first .and. bther == 0) then
    allocate (f0(nnod,ng), fs0(nnod), s0(nnod,ng))
    Ke = 1._dp
    f0 = 1._dp
    call FSrc(fs0)
    first = .false.
  end if

  ! Initialize keff and fission source
  f = Integrate(fs0)
  errn = 1._DP
  e1 = Integrate(errn)

  fn = get_time()
  fdm_time = fdm_time + (fn-st) ! Get FDM time

  !Start outer iteration
  do p=1, nout
    st = get_time()
    fc   = f         ! Save old integrated fission source
    fs0c = fs0       ! Save old fission source
    f0c  = f0        ! Save old flux
    Keo  = Ke        ! Save old multiplication factor
    erro = errn      ! Save old fission source error/difference
    do g = 1, ng
      !!!Calculate total source
      call TSrc(g, Ke, bs)

      !!!Inner Iteration
      call bicg(nin, g, bs, f0(:,g))
    end do
    CALL FSrc (fs0)               !Update fission source
    errn = fs0 - fs0c
    e2 = l2norm(errn)
    if (MOD(p,nac) == 0) call fiss_extrp(popt,e1,e2,erro,errn,fs0)     ! Fission source extrapolation
    e1 = e2                       ! Save l2 norm of the fission source error
    f = Integrate(fs0)            ! Integrate fission source
    Ke = Keo * f / fc             ! Update Keff
    CALL RelE(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
    CALL RelEg(f0, f0c, fer)      ! Search maximum point wise flux error
    fn = get_time()
    fdm_time = fdm_time + (fn-st) ! Get FDM time
    if (MOD(p,nupd) == 0 .and. kern /= ' FDM') call nodal_upd(popt, 1)  ! Nodal coefficients update
    if (popt > 0) then
      write(ounit,'(I5,F13.6,2ES15.5)') p, Ke, ser, fer                 ! Write outer iteration evolution
      if (scr) write(*,'(I5,F13.6,2ES15.5)') p, Ke, ser, fer            ! Write outer iteration evolution
    end if
    if ((ser < serc) .AND. (fer < ferc) .AND. (ndmax < 1.e-2)) exit
  end do

  if (p-1 == nout) THEN
    write(*,*)
    write(*,*) '  MAXIMUM NUMBER OF OUTER ITERATION IS REACHED IN FORWARD CALCULATION.'
    write(*,*) '  CHECK PROBLEM SPECIFICATION OR CHANGE ITERATION CONTROL (%ITER).'
    write(*,*) '  PERHAPS BY MAKING FISSION SOURCE INTERPOLATION MORE FREQUENT'
    write(*,*) '  ADPRES IS STOPING...'
    STOP
  end if

  call print_keff(popt)

  end subroutine outer

  !****************************************************************************!

  subroutine outer_fs(popt)

  !
  ! Purpose:
  !    To perform fixed-source outer iteration


  USE sdata, ONLY: ng, nnod, nout, nin, serc, ferc, fer, ser, f0, nupd, &
                   Ke, nac, fs0, s0, ndmax, kern, get_time, fdm_time
  USE io,    ONLY: ounit, scr
  USE nodal, ONLY: nodal_update

  implicit none

  integer, optional, intent(in) :: popt

  real(dp), dimension(nnod)    :: fs0c               !old fission source
  real(dp), dimension(nnod,ng) :: f0c                !Old flux
  real(dp), dimension(nnod)    :: bs                 !total source
  integer  :: p, g
  real(dp) :: e1, e2
  real(dp), dimension(nnod) :: errn, erro            ! current and past error vectors
  logical  :: first = .true.
  real(dp) :: st, fn

  st = get_time()

  !Setup CMFD matrix
  CALL matrix_setup(1)

  !Allocate flux and fission source for first time
  if (first) then
    allocate (f0(nnod,ng), fs0(nnod), s0(nnod,ng))
    Ke = 1._dp
    f0 = 1._dp
    call FSrc(fs0)
    first = .false.
  end if

  ! Initialize keff and fission source
  errn = 1._DP
  e1 = Integrate(errn)

  fn = get_time()
  fdm_time = fdm_time + (fn-st) ! Get FDM time

  !Start outer iteration
  do p=1, nout
    st = get_time()
    fs0c = fs0       ! Save old fission source
    f0c  = f0        ! Save old flux
    erro = errn       ! Save old fission source error/difference
    do g = 1, ng
      !!!Calculate total source
      call TSrc(g, Ke, bs)

      !!!Inner Iteration
      call bicg(nin, g, bs, f0(:,g))
    end do
    CALL FSrc (fs0)               !Update fission source
    errn = fs0 - fs0c
    e2 = l2norm(errn)
    if (MOD(p,nac) == 0) call fiss_extrp(popt,e1,e2,erro,errn,fs0)     ! Fission source extrapolation
    e1 = e2                       ! Save l2 norm of the fission source error
    CALL RelE(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
    CALL RelEg(f0, f0c, fer)      ! Search maximum point wise flux error
    fn = get_time()
    fdm_time = fdm_time + (fn-st) ! Get FDM time
    if (MOD(p,nupd) == 0 .and. kern /= ' FDM') call nodal_upd(popt, 1)  ! Nodal coefficients update
    if (popt > 0) then
      write(ounit,'(I5,2ES15.5)') p, ser, fer            ! Write outer iteration evolution
      if (scr) write(*,'(I5,2ES15.5)') p, ser, fer      ! Write outer iteration evolution
    end if
    if ((ser < serc) .AND. (fer < ferc) .AND. (ndmax < 1.e-2)) exit
  end do

  if (p-1 == nout) THEN
    write(*,*)
    write(*,*) '  MAXIMUM NUMBER OF OUTER ITERATION IS REACHED IN FIXED-SOURCE CALCULATION.'
    write(*,*) '  CHECK PROBLEM SPECIFICATION OR CHANGE ITERATION CONTROL (%ITER).'
    write(*,*) '  PERHAPS BY MAKING FISSION SOURCE INTERPOLATION MORE FREQUENT'
    write(*,*) '  ADPRES IS STOPING...'
    STOP
  end if

  end subroutine outer_fs

  !****************************************************************************!

  subroutine outer_ad(popt)

  !
  ! Purpose:
  !    To perform adjoint outer iteration


  USE sdata, ONLY: ng, nnod, nout, nin, serc, ferc, fer, ser, f0, nupd, &
                   Ke, nac, fs0, s0, ndmax, kern, get_time, fdm_time
  USE io,    ONLY: ounit, scr
  USE nodal, ONLY: nodal_update

  implicit none

  integer, optional, intent(in) :: popt

  real(dp)                     :: Keo                !Old Multiplication factor (Keff)
  real(dp), dimension(nnod)    :: fs0c               !old fission source
  real(dp), dimension(nnod,ng) :: f0c                !Old flux
  real(dp), dimension(nnod)    :: bs                 !total source
  real(dp) :: f, fc                                  ! new and old integrated fission sources
  integer  :: p, g
  real(dp) :: e1, e2
  real(dp), dimension(nnod) :: errn, erro            ! current and past error vectors
  logical  :: first = .true.
  real(dp) :: st, fn

  st = get_time()

  !Setup CMFD matrix
  CALL matrix_setup(1)

  !Allocate flux and fission source for first time
  if (first .and. popt > 0) then
    allocate (f0(nnod,ng), fs0(nnod), s0(nnod,ng))
    Ke = 1._dp
    f0 = 1._dp
    call FSrcAd(fs0)
    first = .false.
  end if

  ! Initialize keff and fission source
  f = Integrate(fs0)
  errn = 1._DP
  e1 = Integrate(errn)

  fn = get_time()
  fdm_time = fdm_time + (fn-st) ! Get FDM time

  !Start outer iteration
  do p=1, nout
    st = get_time()
    fc   = f         ! Save old integrated fission source
    fs0c = fs0       ! Save old fission source
    f0c  = f0        ! Save old flux
    Keo  = Ke        ! Save old multiplication factor
    erro = errn       ! Save old fission source error/difference
    do g = ng,1,-1
      !!!Calculate total source
      call TSrcAd(g, Ke, bs)

      !!!Inner Iteration
      call bicg(nin, g, bs, f0(:,g))
    end do
    CALL FSrcAd (fs0)               !Update fission source
    errn = fs0 - fs0c
    e2 = l2norm(errn)
    if (MOD(p,nac) == 0) call fiss_extrp(popt,e1,e2,erro,errn,fs0)     ! Fission source extrapolation
    e1 = e2                       ! Save l2 norm of the fission source error
    f = Integrate(fs0)            ! Integrate fission source
    Ke = Keo * f / fc             ! Update Keff
    CALL RelE(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
    CALL RelEg(f0, f0c, fer)      ! Search maximum point wise flux error
    fn = get_time()
    fdm_time = fdm_time + (fn-st) ! Get FDM time
    ! For RODEJECT mode, adjoint calculation is approximated using
    ! nodal coefficients from forward calculation.
    if (MOD(p,nupd) == 0 .and. kern /= ' FDM' &
    .and. popt > 0) call nodal_upd(popt, 0)                       ! Nodal coefficients update
    if (popt > 0) then
      write(ounit,'(I5,F13.6,2ES15.5)') p, Ke, ser, fer           ! Write outer iteration evolution
      if (scr) write(*,'(I5,F13.6,2ES15.5)') p, Ke, ser, fer      ! Write outer iteration evolution
    end if
    if ((ser < serc) .AND. (fer < ferc) .AND. (ndmax < 1.e-2)) exit
  end do

  if (p-1 == nout) THEN
    write(*,*)
    write(*,*) '  MAXIMUM NUMBER OF OUTER ITERATION IS REACHED IN ADJOINT CALCULATION.'
    write(*,*) '  CHECK PROBLEM SPECIFICATION OR CHANGE ITERATION CONTROL (%ITER).'
    write(*,*) '  PERHAPS BY MAKING FISSION SOURCE INTERPOLATION MORE FREQUENT'
    write(*,*) '  ADPRES IS STOPING...'
    STOP
  end if

  call print_keff(popt)

  end subroutine outer_ad

  !****************************************************************************!

  subroutine outer_th(maxn)

  !
  ! Purpose:
  !    To perform TH outer iteration


  USE sdata, ONLY: ng, nnod, nin, serc, ferc, fer, ser, f0, nupd, &
                   Ke, nac, s0, fs0, ndmax, nth, kern, get_time, fdm_time
  USE io,    ONLY: ounit, biter
  USE nodal, ONLY: nodal_update

  implicit none

  integer, intent(in) :: maxn                        ! Max. number of outer iteration

  real(dp)                     :: Keo                !Old Multiplication factor (Keff)
  real(dp), dimension(nnod)    :: fs0c               !old fission source
  real(dp), dimension(nnod,ng) :: f0c                !Old flux
  real(dp), dimension(nnod)    :: bs                 !total source
  real(dp) :: f, fc                                  ! new and old integrated fission sources
  integer  :: p, g
  real(dp) :: e1, e2
  real(dp), dimension(nnod) :: errn, erro            ! current and past error vectors
  logical  :: first = .true.
  real(dp) :: st, fn
  logical  :: lnupd

  st = get_time()

  !Setup CMFD matrix
  CALL matrix_setup(1)

  !Allocate flux and fission source for first time
  if (first) then
    allocate (f0(nnod,ng), fs0(nnod), s0(nnod,ng))
    Ke = 1._dp
    f0 = 1._dp
    call FSrc(fs0)
    first = .false.
  end if

  ! Initialize keff and fission source
  f = Integrate(fs0)
  errn = 1._DP
  e1 = Integrate(errn)

  fn = get_time()
  fdm_time = fdm_time + (fn-st) ! Get FDM time

  if (biter == 0) nupd = int(nth/2)

  !Start outer iteration
  do p=1, maxn
    st = get_time()
    fc   = f         ! Save old integrated fission source
    fs0c = fs0       ! Save old fission source
    f0c  = f0        ! Save old flux
    Keo  = Ke        ! Save old multiplication factor
    erro = errn       ! Save old fission source error/difference
    do g = 1, ng
      !!!Calculate total source
      call TSrc(g, Ke, bs)

      !!!Inner Iteration
      call bicg(nin, g, bs, f0(:,g))
    end do
    CALL FSrc (fs0)               !Update fission source
    errn = fs0 - fs0c
    e2 = l2norm(errn)
    if (MOD(p,nac) == 0) call fiss_extrp(0,e1,e2,erro,errn,fs0)     ! Fission source extrapolation
    e1 = e2                       ! Save l2 norm of the fission source error
    f = Integrate(fs0)            ! Integrate fission source
    Ke = Keo * f / fc             ! Update Keff
    CALL RelE(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
    CALL RelEg(f0, f0c, fer)      ! Search maximum point wise flux error
    fn = get_time()
    fdm_time = fdm_time + (fn-st) ! Get FDM time
    if (MOD(p,nupd) == 0 .and. kern /= ' FDM') then
      lnupd = .false.
      call nodal_upd(0, 1)        ! Nodal coefficients update
    end if
    if ((ser < serc) .AND. (fer < ferc) .AND. (ndmax < 1.e-2)) exit
  end do

  if (lnupd .and. kern /= ' FDM') then
    write(*,*) 'ERROR: OUTER ITERATION WITHIN T-H ITERATION FINISHED WITHOUT NODAL UPDATE'
    write(*,*) 'CHANGE ITERATION CONTROL USING %ITER CARD'
    write(ounit,*) 'ERROR: OUTER ITERATION WITHIN T-H ITERATION FINISHED WITHOUT NODAL UPDATE'
    write(ounit,*) 'CHANGE ITERATION CONTROL USING %ITER CARD'
    stop
  end if

  end subroutine outer_th

  !****************************************************************************!

  subroutine outer_tr(ht,maxi)

  !
  ! Purpose:
  !    To perform transient fixed source outer iteration


  USE sdata, ONLY: ng, nnod, nout, nin, serc, ferc, fer, ser, f0, nupd, &
                   nac, fs0, ndmax, exsrc, kern, get_time, fdm_time
  USE nodal, ONLY: nodal_update

  implicit none

  real(dp), intent(in)        :: ht                ! Time step
  logical, intent(out)         :: maxi              ! Does it reach max number of iteration?

  real(dp), dimension(nnod)    :: fs0c               !old fission source
  real(dp), dimension(nnod,ng) :: f0c                !Old flux
  real(dp), dimension(nnod)    :: bs                 !total source
  integer  :: p, g
  real(dp) :: e1, e2
  real(dp), dimension(nnod) :: errn, erro            ! current and past error vectors
  real(dp) :: st, fn

  st = get_time()

  !Setup CMFD matrix
  CALL matrix_setup(1)

  ! Get terms that do not appear in static calculation
  call get_exsrc(ht, exsrc)

  fn = get_time()
  fdm_time = fdm_time + (fn-st) ! Get FDM time

  !Start outer iteration
  do p=1, nout
    st = get_time()
    fs0c = fs0       ! Save old fission source
    f0c  = f0        ! Save old flux
    erro = errn      ! Save old fission source error/difference
    do g = 1, ng
      !!!Calculate total source
      call TSrcTr(g, bs)

      !!!Inner Iteration
      call bicg(nin, g, bs, f0(:,g))
    end do
    CALL FSrc (fs0)               !Update fission source
    errn = fs0 - fs0c
    e2 = l2norm(errn)
    if (MOD(p,nac) == 0) call fiss_extrp(0,e1,e2,erro,errn,fs0)     ! Fission source extrapolation
    e1 = e2                       ! Save l2 norm of the fission source error
    CALL RelE(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
    CALL RelEg(f0, f0c, fer)      ! Search maximum point wise flux error
    fn = get_time()
    fdm_time = fdm_time + (fn-st) ! Get FDM time
    if (MOD(p,nupd) == 0 .and. kern /= ' FDM') call nodal_upd(0, 2)  ! Nodal coefficients update
    if ((ser < serc) .AND. (fer < ferc) .AND. (ndmax < 1.e-2)) exit
  end do

  if (p==nout+1) THEN
      maxi = .TRUE.
  ELSE
      maxi = .FALSE.
  end if


  end subroutine outer_tr

  !****************************************************************************!

   subroutine get_exsrc(ht, exsrc)

  !
  ! Purpose:
  !    To calculate paramaters from previous time step
  !

  USE sdata, ONLY: lamb, c0, iBeta, chi, mat, velo, &
  fst, ft, nf, m, omeg, ng, nnod, dfis, nuf, &
  bth, sth, s0, tbeta, sigrp, L
  USE io, ONLY: bxtab

  ! Purpose:
     ! To update get external source for transient fixed source problem
     ! This external source is the terms that do not appear in ..
     ! static calculation

  implicit none

  real(dp), intent(in)                      :: ht
  real(dp), dimension(:,:), intent(out)     :: exsrc

  real(dp) :: dt, dtp
  integer  :: n, i, g
  real(dp) :: a1, a2, pxe, pthet

  if (bxtab == 1) then
    dfis = 0._dp
    do n = 1, nnod
      dt = 0._dp; dtp = 0._dp
      do i = 1, nf
        pxe  = exp(-m(mat(n))%lamb(i)*ht)
        if (nuf(n,ng) > 0.) then
          a1   = (1._dp - pxe) / (m(mat(n))%lamb(i)*ht)
        else
          a1 = 0._dp
        end if
        a2   = 1._dp - a1
        a1   = a1 - pxe
        dfis(n) = dfis(n) + m(mat(n))%iBeta(i) * a2
        dt   = dt   + m(mat(n))%lamb(i) * c0(n,i) * pxe &
             + m(mat(n))%iBeta(i) * a1 * fst(n)
        dtp  = dtp + m(mat(n))%lamb(i) * c0(n,i)
      end do

      do g = 1, ng
        pthet = -L(n,g) - sigrp(n,g) * ft(n,g) + s0(n,g) &
              + (1._dp - tbeta(mat(n))) * chi(mat(n),g) * fst(n) &
              +  chi(mat(n),g) * dtp
        exsrc(n,g) = chi(mat(n),g) * dt &
                   + exp(omeg(n,g) * ht) * ft(n,g) &
                   / (sth * m(mat(n))%velo(g) * ht) + bth * pthet
      end do
    end do
  else
    dfis = 0._dp
    do n = 1, nnod
      dt = 0._dp; dtp = 0._dp
      do i = 1, nf
        pxe  = exp(-lamb(i)*ht)
        a1   = (1._dp - pxe) / (lamb(i)*ht)
        a2   = 1._dp - a1
        a1   = a1 - pxe
        dfis(n) = dfis(n) + iBeta(i) * a2
        dt   = dt   + lamb(i) * c0(n,i) * pxe &
             + iBeta(i) * a1 * fst(n)
        dtp  = dtp + lamb(i) * c0(n,i)
      end do

      do g = 1, ng
        pthet = -L(n,g) - sigrp(n,g) * ft(n,g) + s0(n,g) &
              + (1._dp - tbeta(mat(n))) * chi(mat(n),g) * fst(n) &
              +  chi(mat(n),g) * dtp
        exsrc(n,g) = chi(mat(n),g) * dt &
                   + exp(omeg(n,g) * ht) * ft(n,g) / (sth * velo(g) * ht) &
                   + bth * pthet
      end do
    end do
  end if

  end subroutine get_exsrc

 !****************************************************************************!

  subroutine FSrc(fs)
  !
  ! Purpose:
  !   To calculate fission source and fission source moments
  !

  USE sdata, ONLY: nnod, nuf, ng, f0

  implicit none

  real(dp), dimension(:), intent(out) :: fs

  integer :: n, g

  fs = 0._dp
  do g = 1, ng
      do n = 1, nnod
         fs(n)   = fs(n)   + f0 (n,g) * nuf(n,g)
      end do
  end do

  END subroutine FSrc

  !******************************************************************************!

  subroutine FSrcAd(fs)
    !
    ! Purpose:
    !   To calculate fission source (adjoint)
    !

    USE sdata, ONLY: nnod, chi, mat, ng, f0

    implicit none

    real(dp), dimension(:), intent(out) :: fs

    integer :: n, g

    fs = 0._dp
    do g = 1, ng
        do n = 1, nnod
           fs(n)   = fs(n)   + f0 (n,g) * chi(mat(n),g)
        end do
    end do

  end subroutine FSrcAd

  !****************************************************************************!

  subroutine TSrc(g, Keff, bs)
  !
  ! Purpose:
  !   To update total source
  !

  USE sdata, ONLY: chi, mat, nnod, fs0, f0, ng, sigs, chi, s0, exsrc

  implicit none

  integer, intent(in) :: g
  real(dp), intent(in) :: Keff
  real(dp), dimension(:), intent(out) :: bs

  integer :: n, h

  s0 = 0._dp
  do h = 1, ng
      do n = 1, nnod
          if (g /= h) s0(n,g)   = s0(n,g) + sigs(n,h,g) * f0(n,h)
      end do
  end do

  do n = 1, nnod
    bs(n) = chi(mat(n),g) * fs0(n)/Keff  + s0(n,g) + exsrc(n,g)
  end do

  end subroutine TSrc

  !******************************************************************************!

  subroutine TSrcAd(g, Keff, bs)
    !
    ! Purpose:
    !   To update total source (adjoint)
    !

    USE sdata, ONLY: nnod, fs0, f0, ng, sigs, nuf, s0, exsrc

    implicit none

    integer, intent(in) :: g
    real(dp), intent(in) :: Keff
    real(dp), dimension(:), intent(out) :: bs

    integer :: n, h

    s0 = 0._dp
    do h = 1, ng
        do n = 1, nnod
            if (g /= h) s0(n,g)   = s0(n,g) + sigs(n,g,h) * f0(n,h)
        end do
    end do

    do n = 1, nnod
      bs(n) = nuf(n,g) * fs0(n)/Keff  + s0(n,g) + exsrc(n,g)
    end do

  end subroutine TSrcAd

  !****************************************************************************!

  subroutine TSrcTr(g, bs)
  !
  ! Purpose:
  !   To update total source for transient fixed source problem
  !

  USE sdata, ONLY: chi, mat, nnod, fs0, f0, ng, sigs, chi, exsrc, &
  tbeta, s0, dfis

  implicit none

  integer, intent(in)       :: g
  real(dp), dimension(:), intent(out) :: bs

  integer :: n, h

  s0 = 0._dp
  do h = 1, ng
      do n = 1, nnod
          if (g /= h) s0(n,g)   = s0(n,g) + sigs(n,h,g) * f0(n,h)
      end do
  end do

  do n = 1, nnod
    bs(n) =  (1._dp - tbeta(mat(n)) + dfis(n)) * chi(mat(n),g) * fs0(n) &
                + s0(n,g) + exsrc(n,g)
  end do

  end subroutine TSrcTr

  !****************************************************************************!

  function l2norm(a) result(x)

    !purpose: to perform dot product

    real(dp), dimension(:), intent(in)   :: a
    real(dp)                             :: x   ! resulting vector

    integer   :: n

    x = 0._dp
    do n = 1, size(a)
      x = x + a(n)**2
    end do

    x = sqrt(x)

  end function l2norm

  !****************************************************************************!

  FUNCTION Integrate(s) RESULT(intg)

    !
    ! Purpose:
    !    To perform volume integration

  USE sdata, ONLY: nnod, vdel

  implicit none

  real(dp), dimension (:), intent(in) :: s
  real(dp) :: intg
  integer :: n

  intg = 0.
  do n = 1, nnod
      intg = intg + vdel(n) * s(n)
  end do

  end FUNCTION Integrate

  !******************************************************************************!

  subroutine RelE(newF, oldF, rel)

    !
    ! Purpose:
    !    To calculate Max Relative error

  USE sdata, ONLY: nnod

  implicit none

  real(dp), dimension(:), intent(in) :: newF, oldF
  real(dp), intent(out) :: rel

  real(dp) :: error
  integer :: n

  rel = 0.

  do n= 1, nnod
      if (ABS(newF(n)) > 1.e-10_DP) THEN
          error = ABS(newF(n) - oldF(n)) / ABS(newF(n))
          if (error > rel) rel = error
      end if
  end do

  end subroutine RelE

  !******************************************************************************!

  subroutine RelEg(newF, oldF, rel)

    !
    ! Purpose:
    !    To calculate Max Relative error for flux

  USE sdata, ONLY: nnod, ng

  implicit none

  real(dp), dimension(:,:), intent(in) :: newF, oldF
  real(dp), intent(out) :: rel

  real(dp) :: error
  integer :: n, g

  rel = 0.

  do n= 1, nnod
     do g = 1, ng
        if (ABS(newF(n,g)) > 1.d-10) THEN
           error = ABS(newF(n,g) - oldF(n,g)) / ABS(newF(n,g))
           if (error > rel) rel = error
        end if
    end do
  end do

  end subroutine RelEg

  !****************************************************************************!

  subroutine bicg(imax,g,b,x)

    !Purpose: to solve linear of system equation with BiCGSTAB method
    ! (without preconditioner). Sparse matrix saved in a and indexed in rc.
    ! a dimension is (#non_zero_elements)
    ! rc dimension is (2,#non_zero_elements+1)
    ! adapted from:
    ! https://www.cfd-online.com/Wiki/Sample_code_for_BiCGSTAB_-_Fortran_90

    implicit none

    integer, intent(in)                  :: imax, g  ! Max. number of iteration and group number
    real(dp), dimension(:), intent(in)   :: b   ! source
    real(dp), dimension(:), intent(inout)  :: x

    real(dp), dimension(size(b, dim=1))  :: r, rs, v, p, s, t
    real(dp)                             :: rho      , rho_prev
    real(dp)                             :: alpha    , omega   , beta, theta
    integer                              :: i

    r  = b - sp_matvec(g,x)
    rs = r
    rho   = 1.0_dp; alpha = 1.0_dp; omega = 1.0_dp
    v  = 0.0_dp; p  = 0.0_dp

    do i = 1, imax
      rho_prev = rho
      rho      = dproduct(rs,r)
      beta     = (rho/rho_prev) * (alpha/omega)
      p        = r + beta * (p - omega*v)
      v        = sp_matvec(g,p)
      alpha    = rho/dproduct(rs,v)
      s        = r - alpha*v
      t        = sp_matvec(g,s)
      theta    = dproduct(t,t)
      omega    = dproduct(t,s)/theta
      x        = x + alpha*p + omega*s
      r        = s - omega*t
    end do

  end subroutine bicg

 !****************************************************************************!

 function sp_matvec(g,x) result(v)

   USE sdata, ONLY: A, ind, nnod

   !purpose: to perform matrix vector multiplication Axb. A is a square matrix.
   ! Sparse matrix saved in A and indexed in ind

   integer, intent(in) :: g                    ! group number of FDM matrix
   real(dp), dimension(:), intent(in)   :: x   ! vector
   real(dp), dimension(size(x))         :: v   ! resulting vector

   integer   :: n
   integer   :: i

   v = 0._dp
   do n = 1, nnod
     do i = 1, ind(n)%ncol
       v(n) = v(n) + A(n,g)%elmn(i)*x(ind(n)%col(i))
     end do
   end do

 end function sp_matvec

 !****************************************************************************!

 function dproduct(a,b) result(x)

   !purpose: to perform dot product

   real(dp), dimension(:), intent(in)   :: a, b   ! vector
   real(dp)                             :: x   ! resulting vector

   integer   :: n

   x = 0._dp
   do n = 1, size(a)
     x = x + a(n) * b(n)
   end do

 end function dproduct

  !****************************************************************************!

  subroutine PowDis (p)

  !
  ! Purpose:
  !    To calculate power distribution
  !


  USE sdata, ONLY: ng, nnod, sigf, f0, vdel, mode
  USE io,    ONLY: ounit

  implicit none

  real(dp), dimension(:), intent(out) :: p
  integer :: g, n
  real(dp) :: tpow, pow

  p = 0._DP
  do g= 1, ng
      do n= 1, nnod
        pow = f0(n,g) * sigf(n,g) * vdel(n)
        if (pow < 0.) pow = 0.
        p(n) = p(n) + pow
      end do
  end do

  ! Normalize to 1._DP
  tpow = 0._DP
  do n = 1, nnod
      tpow = tpow + p(n)
  end do

  if (tpow <= 0 .AND. mode /= 'FIXEDSRC') THEN
     write(ounit, *) '   ERROR: TOTAL NODES POWER IS ZERO OR LESS'
     write(ounit, *) '   STOP IN subroutine POWDIS'
     STOP
  end if

  do n = 1, nnod
      p(n) = p(n) / tpow
  end do


  end subroutine PowDis

  !****************************************************************************!

end module
