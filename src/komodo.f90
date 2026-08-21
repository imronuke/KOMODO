program main

   use sdata, only: dp, mode, tranw, fdm_time, nod_time, xs_time, &
   inp_time, th_time, get_time
   use io, only: ounit, scr, inp_read, bther
   use control, only: forward, adjoint, fixedsrc, cbsearch, cbsearcht
   use trans, only: rod_eject, rod_eject_th
   use hdf5_output, only: hdf5_close

   implicit none

   real(dp) :: st, fn, tot_time

   ! Read input
   st = get_time()
   call inp_read()
   fn = get_time()
   inp_time = fn - st

   if (scr) then
      write(*, *)
      write(*, *) ' reading input ... done'
   end if

   select case(mode)
      case('FIXEDSRC')
         call fixedsrc()
      case('ADJOINT')
         call adjoint()
      case('RODEJECT')
         if (bther == 0) then
            call rod_eject()
         else
            call rod_eject_th()
         end if
      case('BCSEARCH')
         if (bther == 0) then
            call cbsearch()
         else
            call cbsearcht()
         end if
      case default
         call forward()
   end select

   if (tranw) then
      write(ounit, *)
      write(ounit, *) '  WARNING: ONE OR MORE OUTER ITERATIONS DID NOT CONVERGE.'&
      // 'YOU MAY NEED TO REDUCE TIME STEP'
      write(*, *)
      write(*, *) '  WARNING: ONE OR MORE OUTER ITERATIONS DID NOT CONVERGE.'&
      // 'YOU MAY NEED TO REDUCE TIME STEP'
   end if

   tot_time = fdm_time + th_time + nod_time + xs_time + inp_time
   write(ounit, *) ;
   write(ounit, *) ;
   write(ounit, 1123)
   write(ounit, 1124) inp_time, inp_time / tot_time * 100.
   write(ounit, 1125) xs_time, xs_time / tot_time * 100.
   write(ounit, 1126) fdm_time, fdm_time / tot_time * 100.
   write(ounit, 1127) nod_time, nod_time / tot_time * 100.
   write(ounit, 1128) th_time, th_time / tot_time * 100.
   write(ounit, 1129)
   write(ounit, 1130) tot_time
   if (scr) then
      write(*, *) ;
      write(*, *) ;
      write(*, 1123)
      write(*, 1124) inp_time, inp_time / tot_time * 100.
      write(*, 1125) xs_time, xs_time / tot_time * 100.
      write(*, 1126) fdm_time, fdm_time / tot_time * 100.
      write(*, 1127) nod_time, nod_time / tot_time * 100.
      write(*, 1128) th_time, th_time / tot_time * 100.
      write(*, 1129)
      write(*, 1130) tot_time
   end if
   1123 format (2x, 'CPU time breakdown in seconds')
   1124 format (4x, 'Input reading time   :', f10.4, '  (', f5.1, '%)')
   1125 format (4x, 'XSEC processing time :', f10.4, '  (', f5.1, '%)')
   1126 format (4x, 'CMFD time            :', f10.4, '  (', f5.1, '%)')
   1127 format (4x, 'Nodal update time    :', f10.4, '  (', f5.1, '%)')
   1128 format (4x, 'T-H time             :', f10.4, '  (', f5.1, '%)')
   1129 format (4x, '------------------------------------------')
   1130 format (4x, 'Total time           :', f10.4)

   write(*, *)
   write(*, *) "  KOMODO EXIT NORMALLY"

   call hdf5_close()

end program
