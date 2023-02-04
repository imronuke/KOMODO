program main

    use data
    use control
    use read, only: inp_read
    use time
    
    implicit none
    
    type(timer) :: total
    
    call total % on
    
    ! Read input
    call inp_time % on
    call inp_read()
    call inp_time % off
    
    write(output_unit,*)
    write(output_unit,*) 'READING INPUT ... DONE'
    
    select case(mode)
        case('FIXEDSRC')
            call fixedsrc()
        case('ADJOINT')
            call adjoint()
        case('RODEJECT')
            call rod_eject()
        case('BCSEARCH')
            call cbc_search()
        case default
            call forward()
    end select
    
    ! if (transient_warning) then
    !     write(ounit,*)
    !     write(ounit,*) '  WARNING: ONE OR MORE OUTER ITERATIONS DID NOT CONVERGE.'&
    !                    // 'YOU MAY NEED TO REDUCE TIME STEP'
    !     write(output_unit,*)
    !     write(output_unit,*) '  WARNING: ONE OR MORE OUTER ITERATIONS DID NOT CONVERGE.'&
    !                   // 'YOU MAY NEED TO REDUCE TIME STEP'
    ! end if

    call total % off
    
    write(ounit,*); write(ounit,*); write(ounit,1123)
    write(ounit,1124)inp_time   % elapsed_time, inp_time   % elapsed_time/total % elapsed_time*100.
    write(ounit,1125)xs_time    % elapsed_time,  xs_time   % elapsed_time/total % elapsed_time*100.
    write(ounit,1126)fdm_time   % elapsed_time, fdm_time   % elapsed_time/total % elapsed_time*100.
    write(ounit,1127)nodal_time % elapsed_time, nodal_time % elapsed_time/total % elapsed_time*100.
    write(ounit,1128)th_time    % elapsed_time,  th_time   % elapsed_time/total % elapsed_time*100.
    write(ounit,1129)
    write(ounit,1130)total % elapsed_time

    write(output_unit,*); write(output_unit,*); write(output_unit,1123)
    write(output_unit,1124)inp_time   % elapsed_time, inp_time   % elapsed_time/total % elapsed_time*100.
    write(output_unit,1125)xs_time    % elapsed_time,  xs_time   % elapsed_time/total % elapsed_time*100.
    write(output_unit,1126)fdm_time   % elapsed_time, fdm_time   % elapsed_time/total % elapsed_time*100.
    write(output_unit,1127)nodal_time % elapsed_time, nodal_time % elapsed_time/total % elapsed_time*100.
    write(output_unit,1128)th_time    % elapsed_time,  th_time   % elapsed_time/total % elapsed_time*100.
    write(output_unit,1129)
    write(output_unit,1130)total % elapsed_time

    1123 format (2X,'CPU time breakdown in seconds')
    1124 format (4X,'Input reading time   :', F10.4, '  (', F5.1,'%)')
    1125 format (4X,'XSEC processing time :', F10.4, '  (', F5.1,'%)')
    1126 format (4X,'CMFD time            :', F10.4, '  (', F5.1,'%)')
    1127 format (4X,'Nodal update time    :', F10.4, '  (', F5.1,'%)')
    1128 format (4X,'T-H time             :', F10.4, '  (', F5.1,'%)')
    1129 format (4X,'------------------------------------------')
    1130 format (4X,'Total time           :', F10.4)
    

    
    write(output_unit,*)
    write(output_unit,*) "  KOMODO EXIT NORMALLY"

end PROGRAM
