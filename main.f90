!
! Main for the 3D CG translocon simulations; it will source all necessary 
! modules and call them in the correct order.
!

program main

!===================================================================================
! Source all the modules we will use here.
!===================================================================================

    use sys_param   ! All parameters are defined here
    use polymer     ! Contains subroutines related to polymer properties
    use forceField  ! Contains subroutines to calculate the forces
    use integrator  ! Used to update the positions (velocities are not used in BD)
    use analysis    ! Subroutines to analyze the current state of the system
    use channel ! Translocon Channel properties are defined here
    use iolib   ! Functions related to I/O

!===================================================================================
! Define some variables used in main
!===================================================================================

    implicit none
    integer :: clock_start, clock_end, clock_rate   ! Used to track performance
    integer :: exitflag=0 ! Used to evaluate if the simulation is finished
    double precision :: elapsed_time, disp, xt      ! disp is displacement per
                            ! translation step, xt
                            ! the desired position of the translated bead.
    integer :: ib           ! Used to iterate over loops
    integer(8) :: it=1 ! More memory required for the timestep
    integer :: last_bead ! The last bead to be translated; different if we are using force pulling protocol
    ! *** Added for force pulling
    double precision :: ftmp(NDIM)


!===================================================================================
! Initialize objects
!===================================================================================

    ! First those used to track performance
    call system_clock(COUNT_RATE=clock_rate)
    call system_clock(COUNT=clock_start)

    call read_user_input()  ! Find out the settings desired by the user, iolib.f90
    call init_io()      ! Setup the files that need to be accessed, iolib.f90
    if (logFile.eq.1) then
            open(unit=6,file='log.txt',status='replace')
    endif
    call init_polymer() ! Initialize the polymer, polymer.f90
    call read_channel() ! Read the channel info, iolib.f90
    call allocate_tables()  ! Allocate memory for the force tables, forcefield.f90
    call init_integrator()  ! Initialize the integrator, integrator.f90
    call init_analysis()    ! Pre-create output variables for analysis, analysis.f90
    call init_ff()      ! Set some of the ff parameters, forcefield.f90
    write(6,*) 'Done initializing'

        ! Check if this is a restart
        if (restart/=0) call read_checkpoint(pos(:,:),it,bipBinding(:)) ! iolib.f90

    ! Write the initial coordinates to the trajectory file as well. This way
    ! bonds are displayed correctly in VMD
        if (restart.eq.0) then
                call write_traj(pos(1:nBeads,1:NDIM),aname(1:nBeads)) ! iolib.f90
                call write_LGstatus() ! iolib.f90
        endif

!===================================================================================
! Seperate loop for translation
!===================================================================================



    write(6,*) 'Starting translation!'
    disp = 1.d0/trSteps ! displacement per translation step

    if (ForcePulling .eq. 1) then
        last_bead = nBeads-2
    else 
        last_bead = nBeads-1
    endif

    ! Because I loop to nBeads-2 the last bead is simply released at the end
    do nFree = nInit,last_bead ! ib here counts the last free bead
        nForce = nFree + 1 ! Always calculate forces based on the free beads + 2
!                xt = exitX ! Initially the bead that is being translated is at the ribosomal exit
        do while (it.lt.trSteps+1) ! Timesteps of this translation cycle
            ! We calc forces including a nonfree bead
            call calc_forces(pos(1:nForce,1:NDIM),force(1:nForce,1:NDIM)) ! Calculate forces, forcefield.f90
            ! Update positions based on these forces for all beads that are free
            call brownian_integrate(pos(freeStart:nFree,1:NDIM),force(freeStart:nFree,1:NDIM)) ! integrator.f90
            ! The last bead is simply moved, alternative use the
            ! commented out harmonic potential to move it 
            pos(nForce,1) = pos(nForce,1) + disp
            ! Do LG dynamics based on these new positions
            if (FIXEDLG .eq. 0 ) then
                call LG_dynamics(pos(1:nBeads,1:NDIM)) ! Always use all positions, channel.f90
            endif
            if (modulo(it,write_frequency).lt.1e-5) then
                call write_trajF(pos(1:nBeads,1:NDIM),force(1:nBeads,1:NDIM),aname(1:nBeads))   ! Write out the current positions, iolib.f90
                call write_LGstatus() ! Write the current status of the lateral gate, iolib.f90
            endif
            if (modulo(it,checkpoint_frequency).lt.1e-5) then
                call write_checkpoint(pos,it,bipBinding) ! iolib.f90, for restarts
            endif
            it = it + 1 ! Update the timestep
        enddo
        it = 1 ! *** Re-writing it this way is better for restarts
    enddo
    ! If it is not equal to 1 then we are restarting post translation
    if (it.eq.1) then
            call write_trajF(pos(1:nBeads,1:NDIM),force(1:nBeads,1:NDIM),aname(1:nBeads))
            call write_LGstatus() ! iolib.f90
    endif

    nForce = nBeads ! All beads are also included in the force calculation now

!===================================================================================
! Now simply iterate updating forces and positions
!===================================================================================

write(6,*) 'Translation has finished, continuing dynamics of the free polymer chain'
    do ! Exit conditions are user defined see, read_user_input() in iolib.f90
        ! More steps can be added here
        force = 0.0
   !     bondStatus = 0
        force = 0.0
        call calc_forces(pos(1:nForce,1:NDIM),force(1:nForce,1:NDIM))       
        ! output warning if bond length is too long
   !     if (bondStatus .gt. 0 .and. LGstatus .eq. 0 ) then
   !         write(6, *) 'WARNING: Bond between ',bondStatus,' and &
   !                     ',bondStatus+1,'greater than 1.6 with LG closed at &
   !                     timestep ',it
   !     endif
        call brownian_integrate(pos(freeStart:nFree,1:NDIM),force(freeStart:nFree,1:NDIM))  ! Update the positions, see integrator.f90
        if (FIXEDLG .eq. 0 ) then
            call LG_dynamics(pos(1:nBeads,1:NDIM)) ! Evolve LG state; Always use all positions, channel.f90
        endif
        if (ForcePulling .eq. 1) then
            call do_analysis(pos(1:nBeads,1:NDIM),force(1:nBeads,1:NDIM),ftmp) ! analysis.f90
        endif
        if (modulo(it,write_frequency).lt.1e-5) then
            call write_trajF(pos(1:nBeads,1:NDIM),force(1:nBeads,1:NDIM),aname(1:nBeads))   ! Write out the current positions, iolib.f90
            call write_LGstatus() ! Write the current status of the lateral gate, iolib.f90
            if (ForcePulling .eq. 1) then
                call write_force(ftmp)
            endif
        endif
        if (modulo(it,checkpoint_frequency).lt.1e-5) then
            call write_checkpoint(pos,it,bipBinding) ! iolib.f90, for restarts
        endif
! *** COMMENTED OUT TO ONLY EXIT AT THE MAX SIMULATION TIME (for TatC)
!        call peptide_status(pos(1:nBeads,1:NDIM),exitflag) ! Check if we can terminate, analysis.f90
        if (exitflag/=0) exit ! If we arrived at a final configuration exit
        if (it.gt.MAXSTEPS) exit ! Took too long
        it = it + 1
    enddo
! Write analysis output
    if (ForcePulling .eq. 1) then
        call peptide_status(pos(1:nBeads,1:NDIM),exitflag)
    endif
    call write_analysis(exitflag,it) ! analysis.f90
    call write_trajF(pos(1:nBeads,1:NDIM),force(1:nBeads,1:NDIM),aname(1:nBeads))   ! Write the configuration, iolib.f90
    call write_LGstatus() ! iolib.f90

!===================================================================================
! Finally write out something on the performance
!===================================================================================

call system_clock(count=clock_end)
elapsed_time = dble((clock_end-clock_start)/clock_rate)
write(6,*) "time used: ", elapsed_time
close(6) ! Close logfile

end program main
