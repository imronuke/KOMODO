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
    
    write(*,*)
    write(*,*) 'READING INPUT ... DONE'
    
    select case(mode)
        ! case('FIXEDSRC')
        !     CALL fixedsrc()
        case('ADJOINT')
            CALL adjoint()
        ! case('RODEJECT')
        !     IF (bther == 0) THEN
        !         CALL rod_eject()
        !     ELSE
        !         CALL rod_eject_th()
        !     END IF
        ! case('BCSEARCH')
        !     IF (bther == 0) THEN
        !         CALL cbsearch()
        !     ELSE
        !         CALL cbsearcht()
        !     END IF
        case DEFAULT
            CALL forward()
    END select
    
    ! IF (transient_warning) THEN
    !     write(ounit,*)
    !     write(ounit,*) '  WARNING: ONE OR MORE OUTER ITERATIONS DID NOT CONVERGE.'&
    !                    // 'YOU MAY NEED TO REDUCE TIME STEP'
    !     write(*,*)
    !     write(*,*) '  WARNING: ONE OR MORE OUTER ITERATIONS DID NOT CONVERGE.'&
    !                   // 'YOU MAY NEED TO REDUCE TIME STEP'
    ! END IF

    call total % off
    
    write(ounit,*); write(ounit,*); write(ounit,1123)
    write(ounit,1124)inp_time   % elapsed_time, inp_time   % elapsed_time/total % elapsed_time*100.
    write(ounit,1125)xs_time    % elapsed_time,  xs_time   % elapsed_time/total % elapsed_time*100.
    write(ounit,1126)fdm_time   % elapsed_time, fdm_time   % elapsed_time/total % elapsed_time*100.
    write(ounit,1127)nodal_time % elapsed_time, nodal_time % elapsed_time/total % elapsed_time*100.
    write(ounit,1128)th_time    % elapsed_time,  th_time   % elapsed_time/total % elapsed_time*100.
    write(ounit,1129)
    write(ounit,1130)total % elapsed_time

    write(*,*); write(ounit,*); write(ounit,1123)
    write(*,1124)inp_time   % elapsed_time, inp_time   % elapsed_time/total % elapsed_time*100.
    write(*,1125)xs_time    % elapsed_time,  xs_time   % elapsed_time/total % elapsed_time*100.
    write(*,1126)fdm_time   % elapsed_time, fdm_time   % elapsed_time/total % elapsed_time*100.
    write(*,1127)nodal_time % elapsed_time, nodal_time % elapsed_time/total % elapsed_time*100.
    write(*,1128)th_time    % elapsed_time,  th_time   % elapsed_time/total % elapsed_time*100.
    write(*,1129)
    write(*,1130)total % elapsed_time

    1123 format (2X,'CPU time breakdown in seconds')
    1124 format (4X,'Input reading time   :', F10.4, '  (', F5.1,'%)')
    1125 format (4X,'XSEC processing time :', F10.4, '  (', F5.1,'%)')
    1126 format (4X,'CMFD time            :', F10.4, '  (', F5.1,'%)')
    1127 format (4X,'Nodal update time    :', F10.4, '  (', F5.1,'%)')
    1128 format (4X,'T-H time             :', F10.4, '  (', F5.1,'%)')
    1129 format (4X,'------------------------------------------')
    1130 format (4X,'Total time           :', F10.4)
    

    
    write(*,*)
    write(*,*) "  KOMODO EXIT NORMALLY"

END PROGRAM
