module cmfd

   use sdata, only: dp
   use gpu
   implicit none
   save

   contains

   !****************************************************************************!

   subroutine coup_coef()

      !Purpose: to calculate FDM nodal coupling coefficients

      use sdata, only: ng, nnod, ix, iy, iz, xyz, d, xdel, ydel, zdel, &
      ystag, xstag, nod, nzz, xeast, xwest, ysouth, ynorth, zbott, ztop

      implicit none

      real(dp) :: alb = 1.e30_dp
      real(dp) :: d1, d2
      integer :: n, g, i, j, k
      logical :: first = .true.

      if (first) then
         allocate(nod(nnod, ng))
         do g = 1, ng
            do n = 1, nnod
               nod(n, g)%dn = 0._dp    ! initial nodal coupling coefficients
            end do
         end do
         first = .false.
      end if

      ! Calculate FDM coupling coefficients
      do g = 1, ng
         do n = 1, nnod

            ! Set i, j, k
            i = ix(n) ;
            j = iy(n) ;
            k = iz(n)

            ! Set FDM coupling coefficients in x direction
            if (i == ystag(j)%smax) then
               if (xeast == 0) then
                  nod(n, g)%df(1) = 2._dp * alb * d(n, g) / &
                  (2._dp * d(n, g) + alb * xdel(i))
               else if (xeast == 1) then
                  nod(n, g)%df(1) = d(n, g) / (2._dp * d(n, g) + 0.5_dp * xdel(i))
               else
                  nod(n, g)%df(1) = 0._dp
               end if
            else
               d2 = d(xyz(i + 1, j, k), g)
               nod(n, g)%df(1) = 2._dp * d(n, g) * d2 / &
               (d(n, g) * xdel(i + 1) + d2 * xdel(i))
            end if

            if (i == ystag(j)%smin) then
               if (xwest == 0) then
                  nod(n, g)%df(2) = 2._dp * alb * d(n, g) / &
                  (2._dp * d(n, g) + alb * xdel(i))
               else if (xwest == 1) then
                  nod(n, g)%df(2) = d(n, g) / (2._dp * d(n, g) + 0.5_dp * xdel(i))
               else
                  nod(n, g)%df(2) = 0._dp
               end if
            else
               d1 = d(xyz(i - 1, j, k), g)
               nod(n, g)%df(2) = 2._dp * d(n, g) * d1 / &
               (d(n, g) * xdel(i - 1) + d1 * xdel(i))
            end if

            ! Set nodal coupling coefficients in y direction
            if (j == xstag(i)%smax) then
               if (ynorth == 0) then
                  nod(n, g)%df(3) = 2._dp * alb * d(n, g) / &
                  (2._dp * d(n, g) + alb * ydel(j))
               else if (ynorth == 1) then
                  nod(n, g)%df(3) = d(n, g) / (2._dp * d(n, g) + 0.5_dp * ydel(j))
               else
                  nod(n, g)%df(3) = 0._dp
               end if
            else
               d2 = d(xyz(i, j + 1, k), g)
               nod(n, g)%df(3) = 2._dp * d(n, g) * d2 / &
               (d(n, g) * ydel(j + 1) + d2 * ydel(j))
            end if

            if (j == xstag(i)%smin) then
               if (ysouth == 0) then
                  nod(n, g)%df(4) = 2._dp * alb * d(n, g) / &
                  (2._dp * d(n, g) + alb * ydel(j))
               else if (ysouth == 1) then
                  nod(n, g)%df(4) = d(n, g) / (2._dp * d(n, g) + 0.5_dp * ydel(j))
               else
                  nod(n, g)%df(4) = 0._dp
               end if
            else
               d1 = d(xyz(i, j - 1, k), g)
               nod(n, g)%df(4) = 2._dp * d(n, g) * d1 / &
               (d(n, g) * ydel(j - 1) + d1 * ydel(j))
            end if

            ! Set nodal coupling coefficients in z direction
            if (k == nzz) then
               if (ztop == 0) then
                  nod(n, g)%df(5) = 2._dp * alb * d(n, g) / &
                  (2._dp * d(n, g) + alb * zdel(k))
               else if (ztop == 1) then
                  nod(n, g)%df(5) = d(n, g) / (2._dp * d(n, g) + 0.5_dp * zdel(k))
               else
                  nod(n, g)%df(5) = 0._dp
               end if
            else
               d2 = d(xyz(i, j, k + 1), g)
               nod(n, g)%df(5) = 2._dp * d(n, g) * d2 / &
               (d(n, g) * zdel(k + 1) + d2 * zdel(k))
            end if

            if (k == 1) then
               if (zbott == 0) then
                  nod(n, g)%df(6) = 2._dp * alb * d(n, g) / &
                  (2._dp * d(n, g) + alb * zdel(k))
               else if (zbott == 1) then
                  nod(n, g)%df(6) = d(n, g) / (2._dp * d(n, g) + 0.5_dp * zdel(k))
               else
                  nod(n, g)%df(6) = 0._dp
               end if
            else
               d1 = d(xyz(i, j, k - 1), g)
               nod(n, g)%df(6) = 2._dp * d(n, g) * d1 / &
               (d(n, g) * zdel(k - 1) + d1 * zdel(k))
            end if
         end do
      end do

   end subroutine coup_coef

   !****************************************************************************!

   subroutine set_ind()

      !Purpose: to set indexes for CMFD matrix

      use sdata, only: ix, iy, iz, ystag, xstag, nnod, nxx, nyy, nzz, ind

      implicit none

      integer, dimension(nxx, nyy) :: nodp  !radial node position
      integer :: n, np, rec, i, j, k

      ! setup radial node position nodp
      nodp = 0
      rec = 0
      do j = 1, nyy
         do i = ystag(j)%smin, ystag(j)%smax
            rec = rec + 1
            nodp(i, j) = rec
         end do
      end do

      ! Calculate number of nodes for one planar
      np = rec

      rec = 1
      do n = 1, nnod

         ! Set i, j, k
         i = ix(n) ;
         j = iy(n) ;
         k = iz(n)

         ind%row(n) = rec

         ! Lower diagonal matrix element for z-direction
         if (k /= 1) then
            ind%col(rec) = n - np
            rec = rec + 1
         end if

         ! Lower diagonal matrix element for y-direction
         if (j /= xstag(i)%smin) then
            ind%col(rec) = n - (nodp(i, j) - nodp(i, j - 1))
            rec = rec + 1
         end if

         ! Lower diagonal matrix element for x-direction
         if (i /= ystag(j)%smin) then
            ind%col(rec) = n - 1
            rec = rec + 1
         end if

         ! Diagonal matrix elementss
         ind%col(rec) = n
         rec = rec + 1

         ! Upper diagonal matrix element for x-direction
         if (i /= ystag(j)%smax) then
            ind%col(rec) = n + 1
            rec = rec + 1
         end if

         ! Upper diagonal matrix element for y-direction
         if (j /= xstag(i)%smax) then
            ind%col(rec) = n + (nodp(i, j + 1) - nodp(i, j))
            rec = rec + 1
         end if

         ! Upper diagonal matrix element for z-direction
         if (k /= nzz) then
            ind%col(rec) = n + np
            rec = rec + 1
         end if

      end do
      ind%row(nnod + 1) = rec

      !$acc enter data copyin(ind)
      !$acc enter data copyin(ind%row, ind%col)

   end subroutine set_ind

   !****************************************************************************!

   subroutine matrix_setup(opt)

      !Purpose: to setup sparse penta-diagonal matrix. Elements are indexed in
      ! two-dimensional vector ind and non-zero elements strored in A

      use sdata, only: nod, ix, iy, iz, xdel, ydel, zdel, &
      ystag, xstag, nnod, sigr, nzz, ng, a

      implicit none

      integer, intent(in) :: opt
      integer :: n, g, rec, i, j, k
      logical :: first = .true.

      ! Allocate FDM matrix for first time
      if (first) then
         call set_ind()
         first = .false.
      end if

      ! If need to calculate FDM coupling coefficients
      if (opt > 0) call coup_coef()


      ! Setup CMFD linear system
      do g = 1, ng
         rec = 0
         do n = 1, nnod

            ! Set i, j, k
            i = ix(n) ;
            j = iy(n) ;
            k = iz(n)

            ! Lower diagonal matrix element for z-direction
            if (k /= 1) then
               rec = rec + 1
               a(g)%elmn(rec) = -(nod(n, g)%df(6) - nod(n, g)%dn(6)) / zdel(k)
            end if

            ! Lower diagonal matrix element for y-direction
            if (j /= xstag(i)%smin) then
               rec = rec + 1
               a(g)%elmn(rec) = -(nod(n, g)%df(4) - nod(n, g)%dn(4)) / ydel(j)
            end if

            ! Lower diagonal matrix element for x-direction
            if (i /= ystag(j)%smin) then
               rec = rec + 1
               a(g)%elmn(rec) = -(nod(n, g)%df(2) - nod(n, g)%dn(2)) / xdel(i)
            end if

            ! Diagonal matrix elementss
            rec = rec + 1
            a(g)%elmn(rec) = (nod(n, g)%df(1) + nod(n, g)%df(2) - &
            nod(n, g)%dn(1) + nod(n, g)%dn(2)) / xdel(i) + &
            (nod(n, g)%df(3) + nod(n, g)%df(4) - &
            nod(n, g)%dn(3) + nod(n, g)%dn(4)) / ydel(j) + &
            (nod(n, g)%df(5) + nod(n, g)%df(6) - &
            nod(n, g)%dn(5) + nod(n, g)%dn(6)) / zdel(k) + &
            sigr(n, g)

            ! Upper diagonal matrix element for x-direction
            if (i /= ystag(j)%smax) then
               rec = rec + 1
               a(g)%elmn(rec) = -(nod(n, g)%df(1) + nod(n, g)%dn(1)) / xdel(i)
            end if

            ! Upper diagonal matrix element for y-direction
            if (j /= xstag(i)%smax) then
               rec = rec + 1
               a(g)%elmn(rec) = -(nod(n, g)%df(3) + nod(n, g)%dn(3)) / ydel(j)
            end if

            ! Upper diagonal matrix element for z-direction
            if (k /= nzz) then
               rec = rec + 1
               a(g)%elmn(rec) = -(nod(n, g)%df(5) + nod(n, g)%dn(5)) / zdel(k)
            end if

         end do
      end do

      if (opt > 0) then
         !$acc enter data copyin(A)
         do g = 1, ng
            !$acc enter data copyin(A(g)%elmn)
         end do
      else
         do g = 1, ng
            !$acc update device(A(g)%elmn(:))
         end do
      end if



   end subroutine matrix_setup

   !****************************************************************************!

   subroutine fiss_extrp(popt, e1, e2, erro, errn, fs)

      !
      ! Purpose:
      !    To perform fission source extrapolation

      use io, only: ounit, scr

      implicit none

      integer, intent(in) :: popt
      real(dp), intent(in) :: e1, e2
      real(dp), dimension(:), intent(in) :: erro, errn
      real(dp), dimension(:), intent(inout) :: fs


      real(dp) :: domir, mval

      domir = e2 / e1
      mval = maxval(abs(erro))
      if (mval * mval < 0.0) domir = -domir
      fs = fs + domir / (1._dp - domir) * errn
      if (popt > 0) then
         write(ounit, *) '    ...FISSION SOURCE EXTRAPOLATED...'
         if (scr) write(*, *) '    ...FISSION SOURCE EXTRAPOLATED...'
      end if

   end subroutine fiss_extrp

   !****************************************************************************!

   subroutine nodal_upd(popt, nmode)

      !
      ! Purpose:
      !    To update nodal coupling coefficients

      use sdata, only: ndmax, im, jm, km, kern, get_time, nod_time
      use nodal, only: nodal_update, nodal_update_pnm
      use io, only: ounit, scr


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
         write(ounit, *) '    .....NODAL COUPLING UPDATED..... '
         write(ounit, 1145) ndmax, im, jm, km
         if (scr) then
            write(*, *) '    .....NODAL COUPLING UPDATED..... '
            write(*, 1145) ndmax, im, jm, km
         end if
      end if

      fn = get_time()

      nod_time = nod_time + (fn - st)

      1145 format ('MAX. CHANGE IN NODAL COUPLING COEF.= ', es12.5, &
      ' AT NODE I = ', i2, ', J = ', i2, ', K = ', i2)

   end subroutine nodal_upd

   !****************************************************************************!

   subroutine print_keff()

      !
      ! Purpose:
      !    To update nodal coupling coefficients

      use sdata, only: ke
      use io, only: ounit, scr

      implicit none

      write(ounit, *)
      write(ounit, 1146) ke
      if (scr) then
         write(*, *)
         write(*, 1146) ke
      end if

      1146 format(2x, 'MULTIPLICATION EFFECTIVE (K-EFF) = ', f9.6)

   end subroutine print_keff

   !****************************************************************************!

   subroutine outer(popt)

      !
      ! Purpose:
      !    To perform forward outer iteration


      use sdata, only: ng, nnod, nout, nin, serc, ferc, fer, ser, f0, nupd, &
      ke, nac, fs0, s0, ndmax, kern, get_time, fdm_time, matrix_solver, &
      use_wielandt_shift, wielandt_init_iter, wielandt_shift, a, ind, nuf, chi, mat
      use io, only: ounit, scr, bther
      use nodal, only: nodal_update

      implicit none

      integer, optional, intent(in) :: popt

      real(dp) :: keo                !Old Multiplication factor (Keff)
      real(dp), dimension(nnod) :: fs0c               !old fission source
      real(dp), dimension(nnod) :: bs                 !total source
      real(dp), dimension(nnod, ng) :: f0c                !Old flux
      real(dp) :: f, fc                                  ! new and old integrated fission sources
      integer :: p, g
      real(dp) :: e1, e2
      real(dp), dimension(nnod) :: errn, erro            ! current and past error vectors
      logical :: first = .true.
      real(dp) :: st, fn

      real(dp) :: shift
      real(dp) :: beta ! shifted eigenvalue
      real(dp), allocatable :: shift_diag(:,:) ! (nnod, ng)
      integer :: i, j, row_start, row_end
      real(dp) :: xsum

      st = get_time()

      !Setup CMFD matrix
      call matrix_setup(1)

      !Allocate flux and fission source for first time
      if (first .and. bther == 0) then
         allocate (f0(nnod, ng), fs0(nnod), s0(nnod, ng))
         ke = 1._dp
         f0 = 1._dp
         call fsrc(fs0)
         first = .false.
      end if

      ! Initialize keff and fission source
      f = integrate(fs0)
      errn = 1._dp
      e1 = integrate(errn)

      beta = ke
      if (use_wielandt_shift) then
         allocate(shift_diag(nnod, ng))
         shift_diag = 0.0_dp
      end if

      fn = get_time()
      fdm_time = fdm_time + (fn - st) ! Get FDM time

      !Start outer iteration
      do p = 1, nout

         st = get_time()
         fc = f         ! Save old integrated fission source
         fs0c = fs0       ! Save old fission source
         f0c = f0        ! Save old flux
         keo = ke        ! Save old multiplication factor

         ! update shifted eigenvalue
         if ((use_wielandt_shift) .and. (p > wielandt_init_iter)) then
            if (wielandt_shift < 0.0_dp) then
               shift = ke - wielandt_shift
            else
               shift = wielandt_shift
            end if
            beta = shift * ke / (shift - ke)
         else
            beta = ke
         end if

         if ((use_wielandt_shift) .and. (p > wielandt_init_iter)) then
            ! need to adjust the diagonal for the Wielandt shift
            ! importantly, I also have to store the adjustment so that I can
            ! undo it on the next iteration
            do i = 1,nnod
               row_start = ind%row(i)
               row_end = ind%row(i + 1) - 1
               do j = row_start, row_end
                  if (ind%col(j) == i) then
                     xsum = sum(nuf(i, :) * f0(i, :))
                     do g = 1, ng
                        a(g)%elmn(j) = a(g)%elmn(j) + shift_diag(i, g)
                        shift_diag(i, g) = chi(mat(i), g) / (shift * f0(i, g)) * xsum
                        a(g)%elmn(j) = a(g)%elmn(j) - shift_diag(i, g)
                     end do
                     exit
                  endif
               end do
            end do
         end if

         erro = errn      ! Save old fission source error/difference
         do g = 1, ng
            !!!Calculate total source
            call tsrc(g, beta, bs)

            !!!Inner Iteration
            select case (matrix_solver)
               case ('prec_cg')
                  call prec_cg(nin, g, bs)
               case ('cg')
                  call cg(nin, g, bs)
               case ('bicg')
                   call bicg(nin, g, bs)
            end select
         end do
         call fsrc (fs0)               !Update fission source
         errn = fs0 - fs0c
         e2 = l2norm(errn)
         if (mod(p, nac) == 0) call fiss_extrp(popt, e1, e2, erro, errn, fs0)     ! Fission source extrapolation
         e1 = e2                       ! Save l2 norm of the fission source error
         f = integrate(fs0)            ! Integrate fission source

         ! update shifted eigenvalue
         ! then, recompute keff
         beta = beta * f / fc
         if ((use_wielandt_shift) .and. (p > wielandt_init_iter)) then
            ke = shift * beta / (shift + beta)
         else
            ke = beta
         end if

         call rele(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
         call releg(f0, f0c, fer)      ! Search maximum point wise flux error
         fn = get_time()
         fdm_time = fdm_time + (fn - st) ! Get FDM time
         if (mod(p, nupd) == 0 .and. kern /= ' FDM') call nodal_upd(popt, 1)  ! Nodal coefficients update
         if (popt > 0) then
            write(ounit, '(I5,F13.6,2ES15.5)') p, ke, ser, fer                 ! Write outer iteration evolution
            if (scr) write(*, '(I5,F13.6,2ES15.5)') p, ke, ser, fer            ! Write outer iteration evolution
         end if
         if ((ser < serc) .and. (fer < ferc) .and. (ndmax < 1.e-2)) exit
      end do

      if (p - 1 == nout) then
         write(*, *)
         write(*, *) '  MAXIMUM NUMBER OF OUTER ITERATION IS REACHED IN FORWARD CALCULATION.'
         write(*, *) '  CHECK PROBLEM SPECIFICATION OR CHANGE ITERATION CONTROL (%ITER).'
         write(*, *) '  PERHAPS BY MAKING FISSION SOURCE INTERPOLATION MORE FREQUENT'
         write(*, *) '  KOMODO IS STOPING...'
         stop
      end if

      if (use_wielandt_shift) then
         deallocate(shift_diag)
      endif

   end subroutine outer

   !****************************************************************************!

   subroutine outer_fs(popt)

      !
      ! Purpose:
      !    To perform fixed-source outer iteration


      use sdata, only: ng, nnod, nout, nin, serc, ferc, fer, ser, f0, nupd, &
      ke, nac, fs0, s0, ndmax, kern, get_time, fdm_time, matrix_solver
      use io, only: ounit, scr
      use nodal, only: nodal_update

      implicit none

      integer, optional, intent(in) :: popt

      real(dp), dimension(nnod) :: fs0c               !old fission source
      real(dp), dimension(nnod, ng) :: f0c                !Old flux
      real(dp), dimension(nnod) :: bs                 !total source
      integer :: p, g
      real(dp) :: e1, e2
      real(dp), dimension(nnod) :: errn, erro            ! current and past error vectors
      logical :: first = .true.
      real(dp) :: st, fn

      st = get_time()

      !Setup CMFD matrix
      call matrix_setup(1)

      !Allocate flux and fission source for first time
      if (first) then
         allocate (f0(nnod, ng), fs0(nnod), s0(nnod, ng))
         ke = 1._dp
         f0 = 1._dp
         call fsrc(fs0)
         first = .false.
      end if

      ! Initialize keff and fission source
      errn = 1._dp
      e1 = integrate(errn)

      fn = get_time()
      fdm_time = fdm_time + (fn - st) ! Get FDM time

      !Start outer iteration
      do p = 1, nout
         st = get_time()
         fs0c = fs0       ! Save old fission source
         f0c = f0        ! Save old flux
         erro = errn       ! Save old fission source error/difference
         do g = 1, ng
            !!!Calculate total source
            call tsrc(g, ke, bs)

            !!!Inner Iteration
            select case (matrix_solver)
               case ('prec_cg')
                  call prec_cg(nin, g, bs)
               case ('cg')
                  call cg(nin, g, bs)
               case ('bicg')
                   call bicg(nin, g, bs)
            end select
         end do
         call fsrc (fs0)               !Update fission source
         errn = fs0 - fs0c
         e2 = l2norm(errn)
         if (mod(p, nac) == 0) call fiss_extrp(popt, e1, e2, erro, errn, fs0)     ! Fission source extrapolation
         e1 = e2                       ! Save l2 norm of the fission source error
         call rele(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
         call releg(f0, f0c, fer)      ! Search maximum point wise flux error
         fn = get_time()
         fdm_time = fdm_time + (fn - st) ! Get FDM time
         if (mod(p, nupd) == 0 .and. kern /= ' FDM') call nodal_upd(popt, 1)  ! Nodal coefficients update
         if (popt > 0) then
            write(ounit, '(I5,2ES15.5)') p, ser, fer            ! Write outer iteration evolution
            if (scr) write(*, '(I5,2ES15.5)') p, ser, fer      ! Write outer iteration evolution
         end if
         if ((ser < serc) .and. (fer < ferc) .and. (ndmax < 1.e-2)) exit
      end do

      if (p - 1 == nout) then
         write(*, *)
         write(*, *) '  MAXIMUM NUMBER OF OUTER ITERATION IS REACHED IN FIXED-SOURCE CALCULATION.'
         write(*, *) '  CHECK PROBLEM SPECIFICATION OR CHANGE ITERATION CONTROL (%ITER).'
         write(*, *) '  PERHAPS BY MAKING FISSION SOURCE INTERPOLATION MORE FREQUENT'
         write(*, *) '  KOMODO IS STOPING...'
         stop
      end if

   end subroutine outer_fs

   !****************************************************************************!

   subroutine outer_ad(popt)

      !
      ! Purpose:
      !    To perform adjoint outer iteration


      use sdata, only: ng, nnod, nout, nin, serc, ferc, fer, ser, f0, nupd, &
      ke, nac, fs0, s0, ndmax, kern, get_time, fdm_time, matrix_solver
      use io, only: ounit, scr
      use nodal, only: nodal_update

      implicit none

      integer, optional, intent(in) :: popt

      real(dp) :: keo                !Old Multiplication factor (Keff)
      real(dp), dimension(nnod) :: fs0c               !old fission source
      real(dp), dimension(nnod, ng) :: f0c                !Old flux
      real(dp), dimension(nnod) :: bs                 !total source
      real(dp) :: f, fc                                  ! new and old integrated fission sources
      integer :: p, g
      real(dp) :: e1, e2
      real(dp), dimension(nnod) :: errn, erro            ! current and past error vectors
      logical :: first = .true.
      real(dp) :: st, fn

      st = get_time()

      !Setup CMFD matrix
      call matrix_setup(1)

      !Allocate flux and fission source for first time
      if (first .and. popt > 0) then
         allocate (f0(nnod, ng), fs0(nnod), s0(nnod, ng))
         ke = 1._dp
         f0 = 1._dp
         call fsrcad(fs0)
         first = .false.
      end if

      ! Initialize keff and fission source
      f = integrate(fs0)
      errn = 1._dp
      e1 = integrate(errn)

      fn = get_time()
      fdm_time = fdm_time + (fn - st) ! Get FDM time

      !Start outer iteration
      do p = 1, nout
         st = get_time()
         fc = f         ! Save old integrated fission source
         fs0c = fs0       ! Save old fission source
         f0c = f0        ! Save old flux
         keo = ke        ! Save old multiplication factor
         erro = errn       ! Save old fission source error/difference
         do g = ng, 1, -1
            !!!Calculate total source
            call tsrcad(g, ke, bs)

            !!!Inner Iteration
            select case (matrix_solver)
               case ('prec_cg')
                  call prec_cg(nin, g, bs)
               case ('cg')
                  call cg(nin, g, bs)
               case ('bicg')
                   call bicg(nin, g, bs)
            end select
         end do
         call fsrcad (fs0)               !Update fission source
         errn = fs0 - fs0c
         e2 = l2norm(errn)
         if (mod(p, nac) == 0) call fiss_extrp(popt, e1, e2, erro, errn, fs0)     ! Fission source extrapolation
         e1 = e2                       ! Save l2 norm of the fission source error
         f = integrate(fs0)            ! Integrate fission source
         ke = keo * f / fc             ! Update Keff
         call rele(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
         call releg(f0, f0c, fer)      ! Search maximum point wise flux error
         fn = get_time()
         fdm_time = fdm_time + (fn - st) ! Get FDM time
         ! For RODEJECT mode, adjoint calculation is approximated using
         ! nodal coefficients from forward calculation.
         if (mod(p, nupd) == 0 .and. kern /= ' FDM' &
         .and. popt > 0) call nodal_upd(popt, 0)                       ! Nodal coefficients update
         if (popt > 0) then
            write(ounit, '(I5,F13.6,2ES15.5)') p, ke, ser, fer           ! Write outer iteration evolution
            if (scr) write(*, '(I5,F13.6,2ES15.5)') p, ke, ser, fer      ! Write outer iteration evolution
         end if
         if ((ser < serc) .and. (fer < ferc) .and. (ndmax < 1.e-2)) exit
      end do

      if (p - 1 == nout) then
         write(*, *)
         write(*, *) '  MAXIMUM NUMBER OF OUTER ITERATION IS REACHED IN ADJOINT CALCULATION.'
         write(*, *) '  CHECK PROBLEM SPECIFICATION OR CHANGE ITERATION CONTROL (%ITER).'
         write(*, *) '  PERHAPS BY MAKING FISSION SOURCE INTERPOLATION MORE FREQUENT'
         write(*, *) '  KOMODO IS STOPING...'
         stop
      end if

   end subroutine outer_ad

   !****************************************************************************!

   subroutine outer_th()

      !
      ! Purpose:
      !    To perform TH outer iteration


      use sdata, only: ng, nnod, nin, serc, ferc, fer, ser, f0, nupd, &
      ke, nac, s0, fs0, ndmax, nth, kern, get_time, fdm_time, matrix_solver, &
      use_wielandt_shift, wielandt_init_iter, wielandt_shift, a, ind, nuf, chi, mat
      use io, only: ounit, biter
      use nodal, only: nodal_update

      implicit none

      real(dp) :: keo                !Old Multiplication factor (Keff)
      real(dp), dimension(nnod) :: fs0c               !old fission source
      real(dp), dimension(nnod, ng) :: f0c                !Old flux
      real(dp), dimension(nnod) :: bs                 !total source
      real(dp) :: f, fc                                  ! new and old integrated fission sources
      integer :: p, g
      real(dp) :: e1, e2
      real(dp), dimension(nnod) :: errn, erro            ! current and past error vectors
      logical :: first = .true.
      real(dp) :: st, fn
      logical :: lnupd = .true.

      real(dp) :: shift
      real(dp) :: beta ! shifted eigenvalue
      real(dp), allocatable :: shift_diag(:,:) ! (nnod, ng)
      integer :: i, j, row_start, row_end
      real(dp) :: xsum

      st = get_time()

      !Setup CMFD matrix
      call matrix_setup(1)

      !Allocate flux and fission source for first time
      if (first) then
         allocate (f0(nnod, ng), fs0(nnod), s0(nnod, ng))
         ke = 1._dp
         f0 = 1._dp
         call fsrc(fs0)
         first = .false.
      end if

      ! Initialize keff and fission source
      f = integrate(fs0)
      errn = 1._dp
      e1 = integrate(errn)

      beta = ke
      if (use_wielandt_shift) then
         allocate(shift_diag(nnod, ng))
         shift_diag = 0.0_dp
      end if

      fn = get_time()
      fdm_time = fdm_time + (fn - st) ! Get FDM time

      if (biter == 0) nupd = int(nth / 2)

      !Start outer iteration
      do p = 1, nth

         st = get_time()
         fc = f         ! Save old integrated fission source
         fs0c = fs0       ! Save old fission source
         f0c = f0        ! Save old flux
         keo = ke        ! Save old multiplication factor

         ! update shifted eigenvalue
         if ((use_wielandt_shift) .and. (p > wielandt_init_iter)) then
            if (wielandt_shift < 0.0_dp) then
               shift = ke - wielandt_shift
            else
               shift = wielandt_shift
            end if
            beta = shift * ke / (shift - ke)
         else
            beta = ke
         end if

         if ((use_wielandt_shift) .and. (p > wielandt_init_iter)) then
            ! need to adjust the diagonal for the Wielandt shift
            ! importantly, I also have to store the adjustment so that I can
            ! undo it on the next iteration
            do i = 1,nnod
               row_start = ind%row(i)
               row_end = ind%row(i + 1) - 1
               do j = row_start, row_end
                  if (ind%col(j) == i) then
                     xsum = sum(nuf(i, :) * f0(i, :))
                     do g = 1, ng
                        a(g)%elmn(j) = a(g)%elmn(j) + shift_diag(i, g)
                        shift_diag(i, g) = chi(mat(i), g) / (shift * f0(i, g)) * xsum
                        a(g)%elmn(j) = a(g)%elmn(j) - shift_diag(i, g)
                     end do
                     exit
                  endif
               end do
            end do
         end if

         erro = errn      ! Save old fission source error/difference
         do g = 1, ng
            !!!Calculate total source
            call tsrc(g, beta, bs)

            !!!Inner Iteration
            select case (matrix_solver)
               case ('prec_cg')
                  call prec_cg(nin, g, bs)
               case ('cg')
                  call cg(nin, g, bs)
               case ('bicg')
                   call bicg(nin, g, bs)
            end select
         end do
         call fsrc (fs0)               !Update fission source
         errn = fs0 - fs0c
         e2 = l2norm(errn)
         if (mod(p, nac) == 0) call fiss_extrp(0, e1, e2, erro, errn, fs0)     ! Fission source extrapolation
         e1 = e2                       ! Save l2 norm of the fission source error
         f = integrate(fs0)            ! Integrate fission source

         ! update shifted eigenvalue
         ! then, recompute keff
         beta = beta * f / fc
         if ((use_wielandt_shift) .and. (p > wielandt_init_iter)) then
            ke = shift * beta / (shift + beta)
         else
            ke = beta
         end if

         call rele(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
         call releg(f0, f0c, fer)      ! Search maximum point wise flux error
         fn = get_time()
         fdm_time = fdm_time + (fn - st) ! Get FDM time
         if (mod(p, nupd) == 0 .and. kern /= ' FDM') then
            lnupd = .false.
            call nodal_upd(0, 1)        ! Nodal coefficients update
            shift_diag = 0.0_dp
         end if

         if ((ser < serc) .and. (fer < ferc) .and. (ndmax < 1.e-2)) exit
      end do

      if (lnupd .and. kern /= ' FDM') then
         write(*, *) 'ERROR: OUTER ITERATION WITHIN T-H ITERATION FINISHED WITHOUT NODAL UPDATE'
         write(*, *) 'CHANGE ITERATION CONTROL USING %ITER CARD'
         write(ounit, *) 'ERROR: OUTER ITERATION WITHIN T-H ITERATION FINISHED WITHOUT NODAL UPDATE'
         write(ounit, *) 'CHANGE ITERATION CONTROL USING %ITER CARD'
         stop
      end if

      if (use_wielandt_shift) then
         deallocate(shift_diag)
      endif

   end subroutine outer_th

   !****************************************************************************!

   subroutine outer_tr(ht, maxi)

      !
      ! Purpose:
      !    To perform transient fixed source outer iteration


      use sdata, only: ng, nnod, nout, nin, serc, ferc, fer, ser, f0, nupd, &
      nac, fs0, ndmax, exsrc, kern, get_time, fdm_time, matrix_solver
      use nodal, only: nodal_update

      implicit none

      real(dp), intent(in) :: ht                ! Time step
      logical, intent(out) :: maxi              ! Does it reach max number of iteration?

      real(dp), dimension(nnod) :: fs0c               !old fission source
      real(dp), dimension(nnod, ng) :: f0c                !Old flux
      real(dp), dimension(nnod) :: bs                 !total source
      integer :: p, g
      real(dp) :: e1, e2
      real(dp), dimension(nnod) :: errn, erro            ! current and past error vectors
      real(dp) :: st, fn

      st = get_time()

      !Setup CMFD matrix
      call matrix_setup(1)

      ! Get terms that do not appear in static calculation
      call get_exsrc(ht, exsrc)

      fn = get_time()
      fdm_time = fdm_time + (fn - st) ! Get FDM time

      !Start outer iteration
      do p = 1, nout
         st = get_time()
         fs0c = fs0       ! Save old fission source
         f0c = f0        ! Save old flux
         erro = errn      ! Save old fission source error/difference
         do g = 1, ng
            !!!Calculate total source
            call tsrctr(g, bs)

            !!!Inner Iteration
            select case (matrix_solver)
               case ('prec_cg')
                  call prec_cg(nin, g, bs)
               case ('cg')
                  call cg(nin, g, bs)
               case ('bicg')
                   call bicg(nin, g, bs)
            end select
         end do
         call fsrc (fs0)               !Update fission source
         errn = fs0 - fs0c
         e2 = l2norm(errn)
         if (mod(p, nac) == 0) call fiss_extrp(0, e1, e2, erro, errn, fs0)     ! Fission source extrapolation
         e1 = e2                       ! Save l2 norm of the fission source error
         call rele(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
         call releg(f0, f0c, fer)      ! Search maximum point wise flux error
         fn = get_time()
         fdm_time = fdm_time + (fn - st) ! Get FDM time
         if (mod(p, nupd) == 0 .and. kern /= ' FDM') call nodal_upd(0, 2)  ! Nodal coefficients update
         if ((ser < serc) .and. (fer < ferc) .and. (ndmax < 1.e-2)) exit
      end do

      if (p == nout + 1) then
         maxi = .true.
      else
         maxi = .false.
      end if


   end subroutine outer_tr

   !****************************************************************************!

   subroutine get_exsrc(ht, exsrc)

      use sdata, only: lamb, c0, ibeta, chi, mat, velo, &
      fst, ft, nf, m, omeg, ng, nnod, dfis, nuf, &
      bth, sth, s0, tbeta, sigrp, l
      use io, only: bxtab

      ! Purpose:
      ! To update get external source for transient fixed source problem
      ! This external source contain the terms that appear in the kinetic
      ! calculations but do not appear in static calculation

      implicit none

      real(dp), intent(in) :: ht
      real(dp), dimension(:, :), intent(out) :: exsrc

      real(dp) :: dt, dtp
      integer :: n, i, g
      real(dp) :: a1, a2, pxe, pthet

      if (bxtab == 1) then
         dfis = 0._dp
         do n = 1, nnod
            dt = 0._dp;
            dtp = 0._dp
            do i = 1, nf
               pxe = exp(-m(mat(n))%lamb(i) * ht)
               if (nuf(n, ng) > 0.) then
                  a1 = (1._dp - pxe) / (m(mat(n))%lamb(i) * ht)
               else
                  a1 = 0._dp
               end if
               a2 = 1._dp - a1
               a1 = a1 - pxe
               dfis(n) = dfis(n) + m(mat(n))%ibeta(i) * a2
               dt = dt + m(mat(n))%lamb(i) * c0(n, i) * pxe &
               + m(mat(n))%ibeta(i) * a1 * fst(n)
               dtp = dtp + m(mat(n))%lamb(i) * c0(n, i)
            end do

            do g = 1, ng
               pthet = -l(n, g) - sigrp(n, g) * ft(n, g) + s0(n, g) &
               + (1._dp - tbeta(mat(n))) * chi(mat(n), g) * fst(n) &
               + chi(mat(n), g) * dtp
               exsrc(n, g) = chi(mat(n), g) * dt &
               + exp(omeg(n, g) * ht) * ft(n, g) &
                / (sth * m(mat(n))%velo(g) * ht) + bth * pthet
            end do
         end do
      else
         dfis = 0._dp
         do n = 1, nnod
            dt = 0._dp;
            dtp = 0._dp
            do i = 1, nf
               pxe = exp(-lamb(i) * ht)
               a1 = (1._dp - pxe) / (lamb(i) * ht)
               a2 = 1._dp - a1
               a1 = a1 - pxe
               dfis(n) = dfis(n) + ibeta(i) * a2
               dt = dt + lamb(i) * c0(n, i) * pxe &
               + ibeta(i) * a1 * fst(n)
               dtp = dtp + lamb(i) * c0(n, i)
            end do

            do g = 1, ng
               pthet = -l(n, g) - sigrp(n, g) * ft(n, g) + s0(n, g) &
               + (1._dp - tbeta(mat(n))) * chi(mat(n), g) * fst(n) &
               + chi(mat(n), g) * dtp
               exsrc(n, g) = chi(mat(n), g) * dt &
               + exp(omeg(n, g) * ht) * ft(n, g) / (sth * velo(g) * ht) &
               + bth * pthet
            end do
         end do
      end if

   end subroutine get_exsrc

   !****************************************************************************!

   subroutine fsrc(fs)
      !
      ! Purpose:
      !   To calculate fission source and fission source moments
      !

      use sdata, only: nnod, nuf, ng, f0

      implicit none

      real(dp), dimension(:), intent(out) :: fs

      integer :: n, g

      fs = 0._dp
      do g = 1, ng
         do n = 1, nnod
            fs(n) = fs(n) + f0 (n, g) * nuf(n, g)
         end do
      end do

   end subroutine fsrc

   !******************************************************************************!

   subroutine fsrcad(fs)
      !
      ! Purpose:
      !   To calculate fission source (adjoint)
      !

      use sdata, only: nnod, chi, mat, ng, f0

      implicit none

      real(dp), dimension(:), intent(out) :: fs

      integer :: n, g

      fs = 0._dp
      do g = 1, ng
         do n = 1, nnod
            fs(n) = fs(n) + f0 (n, g) * chi(mat(n), g)
         end do
      end do

   end subroutine fsrcad

   !****************************************************************************!

   subroutine tsrc(g, keff, bs)
      !
      ! Purpose:
      !   To update total source
      !

      use sdata, only: chi, mat, nnod, fs0, f0, ng, sigs, chi, s0, exsrc

      implicit none

      integer, intent(in) :: g
      real(dp), intent(in) :: keff
      real(dp), dimension(:), intent(out) :: bs

      integer :: n, h

      s0 = 0._dp
      do h = 1, ng
         do n = 1, nnod
            if (g /= h) s0(n, g) = s0(n, g) + sigs(n, h, g) * f0(n, h)
         end do
      end do

      do n = 1, nnod
         bs(n) = chi(mat(n), g) * fs0(n) / keff + s0(n, g) + exsrc(n, g)
      end do

   end subroutine tsrc

   !******************************************************************************!

   subroutine tsrcad(g, keff, bs)
      !
      ! Purpose:
      !   To update total source (adjoint)
      !

      use sdata, only: nnod, fs0, f0, ng, sigs, nuf, s0, exsrc

      implicit none

      integer, intent(in) :: g
      real(dp), intent(in) :: keff
      real(dp), dimension(:), intent(out) :: bs

      integer :: n, h

      s0 = 0._dp
      do h = 1, ng
         do n = 1, nnod
            if (g /= h) s0(n, g) = s0(n, g) + sigs(n, g, h) * f0(n, h)
         end do
      end do

      do n = 1, nnod
         bs(n) = nuf(n, g) * fs0(n) / keff + s0(n, g) + exsrc(n, g)
      end do

   end subroutine tsrcad

   !****************************************************************************!

   subroutine tsrctr(g, bs)
      !
      ! Purpose:
      !   To update total source for transient fixed source problem
      !

      use sdata, only: chi, mat, nnod, fs0, f0, ng, sigs, chi, exsrc, &
      tbeta, s0, dfis

      implicit none

      integer, intent(in) :: g
      real(dp), dimension(:), intent(out) :: bs

      integer :: n, h

      s0 = 0._dp
      do h = 1, ng
         do n = 1, nnod
            if (g /= h) s0(n, g) = s0(n, g) + sigs(n, h, g) * f0(n, h)
         end do
      end do

      do n = 1, nnod
         bs(n) = (1._dp - tbeta(mat(n)) + dfis(n)) * chi(mat(n), g) * fs0(n) &
         + s0(n, g) + exsrc(n, g)
      end do

   end subroutine tsrctr

   !****************************************************************************!

   function l2norm(a) result(x)

      !purpose: to perform dot product

      real(dp), dimension(:), intent(in) :: a
      real(dp) :: x   ! resulting vector

      integer :: n

      x = 0._dp
      do n = 1, size(a)
         x = x + a(n)**2
      end do

      x = sqrt(x)

   end function l2norm

   !****************************************************************************!

   function integrate(s) result(intg)

      !
      ! Purpose:
      !    To perform volume integration

      use sdata, only: nnod, vdel

      implicit none

      real(dp), dimension (:), intent(in) :: s
      real(dp) :: intg
      integer :: n

      intg = 0.
      do n = 1, nnod
         intg = intg + vdel(n) * s(n)
      end do

   end function integrate

   !******************************************************************************!

   subroutine rele(newf, oldf, rel)

      !
      ! Purpose:
      !    To calculate Max Relative error

      use sdata, only: nnod

      implicit none

      real(dp), dimension(:), intent(in) :: newf, oldf
      real(dp), intent(out) :: rel

      real(dp) :: error
      integer :: n

      rel = 0.

      do n = 1, nnod
         if (abs(newf(n)) > 1.e-10_dp) then
            error = abs(newf(n) - oldf(n)) / abs(newf(n))
            if (error > rel) rel = error
         end if
      end do

   end subroutine rele

   !******************************************************************************!

   subroutine releg(newf, oldf, rel)

      !
      ! Purpose:
      !    To calculate Max Relative error for flux

      use sdata, only: nnod, ng

      implicit none

      real(dp), dimension(:, :), intent(in) :: newf, oldf
      real(dp), intent(out) :: rel

      real(dp) :: error
      integer :: n, g

      rel = 0.

      do n = 1, nnod
         do g = 1, ng
            if (abs(newf(n, g)) > 1.d-10) then
               error = abs(newf(n, g) - oldf(n, g)) / abs(newf(n, g))
               if (error > rel) rel = error
            end if
         end do
      end do

   end subroutine releg

   !****************************************************************************!

   subroutine prec_cg(imax, g, b)

      use sdata, only: r, rs, p, v, t, ng, nnod, f0, inner_atol, inner_rtol, &
      ind, a

      !Purpose: to solve linear of system equation with Conjugate Gradient (CG) method
      ! with Jacobi (diagonal) preconditioner. Sparse matrix saved in a and indexed in rc.
      ! a dimension is (#non_zero_elements)
      ! rc dimension is (2,#non_zero_elements+1)
      ! on first entry for each group, compute and store the inverse diagonal to use for preconditioner
      ! adapted from:
      ! "Iterative Methods for Linear and Nonlinear Equations" By: C. T. Kelley

      implicit none

      integer, intent(in) :: imax, g  ! Max. number of iteration and group number
      real(dp), dimension(:), intent(in) :: b   ! source

      real(dp) :: rho_k1, rho_k2
      real(dp) :: tau_k1, tau_k2
      real(dp) :: alpha, beta
      real(dp) :: converge, xsto
      integer :: i, j
      integer :: row_start, row_end

      real(dp), allocatable :: inv_diag(:)

      logical, save :: first = .true.
      logical, parameter :: verbose = .false.

      if (first) then
        call gpu_allocate(r, nnod)
        call gpu_allocate(rs, nnod)
        call gpu_allocate(p, nnod)
        call gpu_allocate(v, nnod)
        call gpu_allocate(t, nnod)
        first = .false.
      end if

      ! Construct the preconditioner. This is done every time this subroutine is visited.
      ! It may be possible to efficiently cache this information, but would require more knowledge of the algorithm.
      ! Be cautious if this is used in association with Sutton's method for Wielandt shift or the matrix is updated
      ! during CMFD iterations.
      call gpu_allocate(inv_diag, nnod)
      call gpu_initialize(inv_diag, 0.0_dp)
      do i = 1,nnod
        row_start = ind%row(i)
        row_end   = ind%row(i + 1) - 1
        do j = row_start, row_end
           if (ind%col(j) == i) then
              inv_diag(i) = 1.0_dp / a(g)%elmn(j)
              exit
           end if
        end do
      end do

      ! 1)
      call sp_matvec(g, f0(:, g), r)
      call xpby(b, -1._dp, r, rs)
      call xew(rs, r)

      rho_k2 = 0.0_dp
      rho_k1 = dproduct(r, r)
      xsto = sqrt(rho_k1)

      tau_k1 = 0.0_dp
      tau_k2 = 0.0_dp

      call gpu_initialize(v, 0.0_dp)
      call gpu_initialize(p, 0.0_dp)
      call gpu_initialize(t, 0.0_dp)

      ! 2)
      beta = 0.0_dp
      do i = 1,imax
         ! 2) (a) 
         call xmuly(inv_diag, r, t)
         ! 2) (b)
         tau_k1 = dproduct(t, r)
         ! 2) (c)
         if (i == 1) then
            beta = 0.0_dp
         else
            beta = tau_k1 / tau_k2
         end if
         call xpby(t, beta, p, rs)
         call xew(rs, p)
         ! 2) (d)
         call sp_matvec(g, p, v)
         ! 2) (e)
         alpha = tau_k1 / dproduct(p, v)
         ! 2) (f)
         call xpby(f0(:, g), alpha, p, rs)
         call xew(rs, f0(:, g))
         ! 2) (g)
         call xpby(r, -alpha, v, rs)
         call xew(rs, r)
         ! 2) (h)
         rho_k2 = rho_k1
         rho_k1 = dproduct(r, r)
         tau_k2 = tau_k1

         converge = sqrt(rho_k1)
         if (verbose) then
           write(*,'(a,i0,1x,a,es13.6)') 'it=', i, 'converge=', converge
         end if
         if ((converge < inner_atol) .or. (converge < inner_rtol*xsto)) exit
      end do

      ! NOTE: consider checking if the solution actually converged before
      ! the number of iterations reached imax.

   end subroutine prec_cg

   !****************************************************************************!

   subroutine cg(imax, g, b)

      use sdata, only: r, rs, p, v, nnod, f0, inner_atol, inner_rtol

      !Purpose: to solve linear of system equation with Conjugate Gradient (CG) method
      ! (without preconditioner). Sparse matrix saved in a and indexed in rc.
      ! a dimension is (#non_zero_elements)
      ! rc dimension is (2,#non_zero_elements+1)
      ! adapted from:
      ! "Iterative Methods for Linear and Nonlinear Equations" By: C. T. Kelley

      implicit none

      integer, intent(in) :: imax, g  ! Max. number of iteration and group number
      real(dp), dimension(:), intent(in) :: b   ! source

      real(dp) :: rho_k1, rho_k2
      real(dp) :: alpha, beta
      real(dp) :: converge, xsto
      integer :: i

      logical, save :: first = .true.
      logical, parameter :: verbose = .false.

      if (first) then
        call gpu_allocate(r, nnod)
        call gpu_allocate(rs, nnod)
        call gpu_allocate(p, nnod)
        call gpu_allocate(v, nnod)
        first = .false.
      end if

      ! 1)
      call sp_matvec(g, f0(:, g), r)
      call xpby(b, -1._dp, r, rs)
      call xew(rs, r)

      rho_k2 = 0.0_dp
      rho_k1 = dproduct(r, r)
      xsto = sqrt(rho_k1)

      call gpu_initialize(v, 0.0_dp)
      call gpu_initialize(p, 0.0_dp)

      ! 2)
      beta = 0.0_dp
      do i = 1,imax
         ! 2) (a)
         call xpby(r, beta, p, rs)
         call xew(rs, p)
         ! 2) (b)
         call sp_matvec(g, p, v)
         ! 2) (c)
         alpha = rho_k1 / dproduct(p, v)
         ! 2) (d)
         call xpby(f0(:, g), alpha, p, rs)
         call xew(rs, f0(:, g))
         ! 2) (e)
         call xpby(r, -alpha, v, rs)
         call xew(rs, r)
         ! 2) (f)
         ! store old values
         rho_k2 = rho_k1
         rho_k1 = dproduct(r, r)
         beta = rho_k1 / rho_k2

         converge = sqrt(rho_k1)
         if (verbose) then
           write(*,'(a,i0,1x,a,es13.6)') 'it=', i, 'converge=', converge
         end if
         if ((converge < inner_atol) .or. (converge < inner_rtol*xsto)) exit
      end do

      ! NOTE: consider checking if the solution actually converged before
      ! the number of iterations reached imax.

   end subroutine cg

   !****************************************************************************!

   subroutine bicg(imax, g, b)

      use sdata, only: r, rs, v, p, s, t, tmp, nnod, f0, inner_atol, inner_rtol

      !Purpose: to solve linear of system equation with BiCGSTAB method
      ! (without preconditioner). Sparse matrix saved in a and indexed in rc.
      ! a dimension is (#non_zero_elements)
      ! rc dimension is (2,#non_zero_elements+1)
      ! adapted from:
      ! https://www.cfd-online.com/Wiki/Sample_code_for_BiCGSTAB_-_Fortran_90

      implicit none

      integer, intent(in) :: imax, g  ! Max. number of iteration and group number
      real(dp), dimension(:), intent(in) :: b   ! source

      real(dp) :: rho, rho_prev
      real(dp) :: alpha, omega, beta, theta
      real(dp) :: converge, xsto
      integer :: i

      logical :: first = .true.
      logical, parameter :: verbose = .false.

      if (first) then
         call gpu_allocate(r, nnod)
         call gpu_allocate(rs, nnod)
         call gpu_allocate(v, nnod)
         call gpu_allocate(p, nnod)
         call gpu_allocate(t, nnod)
         call gpu_allocate(s, nnod)
         call gpu_allocate(tmp, nnod)
         !$acc enter data copyin(f0(:,:))
         first = .false.
      end if

      !$acc enter data copyin(b)

      ! write(*,*) r(1:3)
      ! call gpu_initialize(r, 3.0_dp)
      ! !$acc update self(r)
      ! write(*,*) r(1:3)

      call sp_matvec(g, f0(:, g), r)

      call xpby(b, -1._dp, r, rs)
      call xew(rs, r)

      rho = dproduct(rs, r)
      xsto = sqrt(rho)
      alpha = 1.0_dp
      omega = 1.0_dp
      call gpu_initialize(v, 0.0_dp)
      call gpu_initialize(p, 0.0_dp)

      ! NOTE: this algorithim can be modified to allow two "exit" points so that
      ! it would be possible to quit early after the first mat*vec product.
      ! It would require slightly reordering the operations.
      do i = 1, imax
         rho_prev = rho
         rho = dproduct(rs, r)
         beta = (rho / rho_prev) * (alpha / omega)
         call xpby(p, -omega, v, tmp)
         call xpby(r, beta, tmp, p)
         call sp_matvec(g, p, v)
         alpha = rho / dproduct(rs, v)
         call xpby(r, -alpha, v, s)
         call sp_matvec(g, s, t)
         theta = dproduct(t, t)
         omega = dproduct(t, s) / theta
         call axpby(alpha, p, omega, s, tmp)
         call xpby(f0(:, g), 1.0_dp, tmp, r)
         call xew(r, f0(:, g))
         call xpby(s, -omega, t, r)
         converge = sqrt(dproduct(r, r))
         if (verbose) then
           write(*,'(a,i0,1x,a,es13.6)') 'it=', i, 'converge=', converge
         end if
         if ((converge < inner_atol) .or. (converge < inner_rtol*xsto)) exit
      end do
      !$acc exit data delete(b)
      !$acc update self(f0(:,g))

      ! NOTE: consider checking if the solution actually converged before
      ! the number of iterations reached imax.

   end subroutine bicg

   !****************************************************************************!

   subroutine sp_matvec(g, x, vx)

      use sdata, only: a, ind, nnod

      !purpose: to perform matrix vector multiplication Axb. A is a square matrix.
      ! Sparse matrix saved in A and indexed in ind

      integer, intent(in) :: g                    ! group number of FDM matrix
      real(dp), dimension(:), intent(in) :: x   ! vector
      real(dp), dimension(:) :: vx   ! resulting vector

      integer :: i, j
      integer :: row_start, row_end
      real(dp) :: tmpsum

      !$acc kernels present(A(g)%elmn, ind, x, vx)
      !$acc loop device_type(nvidia) gang worker(32)
      do i = 1, nnod
         tmpsum = 0.
         row_start = ind%row(i)
         row_end = ind%row(i + 1) - 1
         ! vector(8) because at each row there are 7 non zero elements
         !$acc loop device_type(nvidia) vector(8)
         do j = row_start, row_end
            tmpsum = tmpsum + a(g)%elmn(j) * x(ind%col(j))
         end do
         vx(i) = tmpsum
      end do
      !$acc end kernels

   end subroutine sp_matvec

   !****************************************************************************!

   real(dp) function dproduct(a, b)
      !purpose: to perform dot product

      real(dp), dimension(:), intent(in) :: a, b   ! vector
      real(dp) :: x   ! resulting vector

      integer :: n


      !$acc kernels present(a, b)
      x = 0._dp
      do n = 1, size(a)
         x = x + a(n) * b(n)
      end do
      !$acc end kernels

      dproduct = x

   end function dproduct

   !****************************************************************************!

   subroutine axpby(alpha, x, beta, y, w)
      implicit none
      real(dp) :: alpha, beta, x(:), y(:)
      real(dp) :: w(:)
      integer :: i, length

      length = size(x)
      !$acc kernels present(x,y,w)
      do i = 1, length
         w(i) = alpha * x(i) + beta * y(i)
      end do
      !$acc end kernels

   end subroutine

   !****************************************************************************!

   subroutine xpby(x, beta, y, w)
      implicit none
      real(dp) :: beta, x(:), y(:)
      real(dp) :: w(:)
      integer :: i, length

      length = size(x)
      !$acc kernels present(x,y,w)
      do i = 1, length
         w(i) = x(i) + beta * y(i)
      end do
      !$acc end kernels

   end subroutine

   !****************************************************************************!

   subroutine xmuly(x, y, z)
      implicit none
      real(dp), intent(in) :: x(:), y(:)
      real(dp), intent(out) :: z(:)
      integer :: i, length

      length = size(x)
      !$acc kernels present(x,y,z)
      do i = 1, length
        z(i) = x(i) * y(i)
      end do
      !$acc end kernels

   end subroutine

   !****************************************************************************!

   subroutine xew(x, w)
      implicit none
      real(dp) :: x(:)
      real(dp) :: w(size(x))
      integer :: i, length

      length = size(x)
      !$acc kernels present(x,w)
      do i = 1, length
         w(i) = x(i)
      end do
      !$acc end kernels

   end subroutine

end module
