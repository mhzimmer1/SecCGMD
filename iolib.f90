!==========================================================================================
! Module that deals with I/O. Contains subroutines for reading command line arguments,
! defines and creates output files, writes to files.
!==========================================================================================

module iolib

    use sys_param   ! Contains some default fileIDs
    implicit none
    ! Here I define default filenames for files that are always required
    ! if the user provides alternate filenames these defaults will be
    ! overwritten
    character(250) ::   trajFile='traj.xyz',    &! File used to store the trajectory
                pulling_trajFile='pulling_traj.xyz', &
                confInFile='conf0.xyz', &! File with the input coordinates
                confOutFile='conf1.xyz',&! File for output coordinates
                seqInFile='seq.txt',    &! File with the bead properties (name, hydrophobicity)
                checkInFile='checkpoint.dat', & ! Checkpoint data
                checkOutFile='checkpoint.dat', & ! Checkpoint data
                chanInFile='chan.xyz', &    ! File with the channel
                pairwiseFile='pairwise.txt', &    ! File with pairwise interactions between beads
                localizationFile='localization.txt'    ! File with localization restraints for beads
    ! Next default parameters concerning the number of steps and the write frequency
    integer :: nsteps=1000, write_frequency=100, pulling_write_frequency=100, checkpoint_frequency=1000,    &! How many iterations and how often write coords
                        logFile=0               ! Write to file or commandline(0)
        ! For OS calls
        character(250) :: command

    contains

!==========================================================================================
! Subroutine read_user_input() reads the command line arguments provided when calling
! the program.
!==========================================================================================

    subroutine read_user_input()

        ! Space to initialize variables used only here
        implicit none
        integer :: nargs, iarg=0
        character(250) :: arg

        ! Count how many arguments there are
        nargs = command_argument_count()

        ! For each command argument see what it specifies
        do
            call get_command_argument(iarg,arg)
            if (len_trim(arg) == 0) then
                exit ! We have finished reading all arguments
            endif
            if (trim(arg) == '-c') then ! The next input is the initial configuration file
                call get_command_argument(iarg+1,arg)
                read(arg,'(A250)') confInFile
                iarg = iarg + 1 ! Update this since the next argument is not a flag
            endif
            if (trim(arg) == '-seq') then ! The next input is a file with polymer bead properties
                call get_command_argument(iarg+1,arg)
                read(arg,'(A250)') seqInFile
                iarg = iarg + 1
            endif
            if (trim(arg) == '-localization') then ! Parameter for the attractive interaction between TMDs (kT)
                call get_command_argument(iarg+1,arg)
                read(arg,'(A250)') localizationFile ! sys_param.f90, polymer.f90
                iarg = iarg + 1 ! Update this since the next argument is not a flag
            endif
            if (trim(arg) == '-pairwise') then ! Parameter for the attractive interaction between TMDs (kT)
                call get_command_argument(iarg+1,arg)
                read(arg,'(A250)') pairwiseFile ! sys_param.f90, polymer.f90
                iarg = iarg + 1 ! Update this since the next argument is not a flag
            endif
            if (trim(arg) == '-tmpotU') then ! Parameter for magnitude of the transmembrane potential
                call get_command_argument(iarg+1,arg)
                read(arg,*) tmpU ! sys_param.f90, forceField.f90
                iarg = iarg + 1 ! Update this since the next argument is not a flag
            endif
            if (trim(arg) == '-tmpotL') then ! Parameter for the lengthscale of the transmembrane potential
                call get_command_argument(iarg+1,arg)
                read(arg,*) tmpL ! sys_param.f90, forceField.f90
                iarg = iarg + 1 ! Update this since the next argument is not a flag
            endif
            if (trim(arg) == '-BIP') then ! Change explicit bip binding, 1 is on, 0 is off (default)
                call get_command_argument(iarg+1,arg)
                read(arg,*) BIPON
                iarg = iarg + 1
            endif
            if (trim(arg) == '-F_bip') then ! Change bip on rate 
                call get_command_argument(iarg+1,arg)
                read(arg,*) F_bip
                iarg = iarg + 1
            endif
            if (trim(arg) == '-k_bip_on') then ! Change bip on rate 
                call get_command_argument(iarg+1,arg)
                read(arg,*) k_bip_on
                iarg = iarg + 1
            endif
            if (trim(arg) == '-k_bip_off') then ! Change bip off rate 
                call get_command_argument(iarg+1,arg)
                read(arg,*) k_bip_off
                iarg = iarg + 1
            endif
            if (trim(arg) == '-softcore') then ! Parameter for the LJ softness (0=hardcore)
                call get_command_argument(iarg+1,arg)
                read(arg,*) alphaLJ ! sys_param.f90, forceField.f90
                iarg = iarg + 1 ! Update this since the next argument is not a flag
            endif
            if (trim(arg) == '-chan') then ! The next input is the initial configuration file
                call get_command_argument(iarg+1,arg)
                read(arg,'(A250)') chanInFile
                iarg = iarg + 1 ! Update this since the next argument is not a flag
            endif
            if (trim(arg) == '-log') then ! All output will be written to a logfile instead of the command line
                call get_command_argument(iarg+1,arg)
                read(arg,'(I1)') logFile
                iarg = iarg + 1
            endif
            if (trim(arg) == '-o') then ! The next input is the final configuration file
                call get_command_argument(iarg+1,arg)
                read(arg,'(A250)') confOutFile
                iarg = iarg + 1 ! Next argument is not a flag
            endif
            if (trim(arg) == '-t') then ! The next input is the trajectory file
                call get_command_argument(iarg+1,arg)
                read(arg,'(A250)') trajFile
                iarg = iarg + 1 ! Next argument is not a flag
            endif
            if (trim(arg) == '-dgfillo') then !parameter for filling in open channel
                call get_command_argument(iarg+1,arg)
                read(arg,*) DGFILLO
                iarg = iarg + 1
            endif
            if (trim(arg) == '-dgfillc') then !parameter for filling in closed channel
                call get_command_argument(iarg+1,arg)
                read(arg,*) DGFILLC
                iarg = iarg + 1
            endif
            if (trim(arg) == '-dgempty') then !parameter for solvation energy in open channel
                call get_command_argument(iarg+1,arg)
                read(arg,*) DGEMPTY
                iarg = iarg + 1
            endif
            if (trim(arg) == '-fixedLG') then ! don't allow LG to change state
                FIXEDLG = 1
                !iarg = iarg + 1
            endif
            if (trim(arg) == '-LGstatus') then ! set starting LGstatus
                call get_command_argument(iarg+1,arg)
                read(arg,*) LGstatus
                iarg = iarg + 1
            endif
            if (trim(arg) == '-seed') then ! The next input is number of frames
                call get_command_argument(iarg+1,arg)
                read(arg,'(I15)') extraRanSeed
                iarg = iarg + 1 ! Next argument is not a flag
            endif
            if (trim(arg) == '-nstep') then ! The next input is number of frames
                call get_command_argument(iarg+1,arg)
                read(arg,'(I15)') MAXSTEPS 
                iarg = iarg + 1 ! Next argument is not a flag
            endif
            if (trim(arg) == '-freq') then ! Next input is write frequency
                call get_command_argument(iarg+1,arg)
                read(arg,'(I15)') write_frequency
                iarg = iarg + 1 ! Next argument is not a flag
            endif
            if (trim(arg) == '-pulling_freq') then ! Next input is write frequency
                call get_command_argument(iarg+1,arg)
                read(arg,'(I15)') pulling_write_frequency
                iarg = iarg + 1 ! Next argument is not a flag
            endif
            if (trim(arg) == '-cfreq') then ! Next input is checkpoint frequency
                call get_command_argument(iarg+1,arg)
                read(arg,'(I15)') checkpoint_frequency
                iarg = iarg + 1 ! Next argument is not a flag
            endif
            if (trim(arg) == '-nb') then ! Next input is the number of polymer beads
                call get_command_argument(iarg+1,arg)
                read(arg,'(I15)') nBeads
                iarg = iarg + 1
            endif
            if (trim(arg) == '-ni') then ! Next input is the number of polymer beads already translated
                call get_command_argument(iarg+1,arg)
                read(arg,'(I15)') nInit
                iarg = iarg + 1
            endif
            if (trim(arg) == '-nf') then ! Number of beads on the N-terminal end that are not allowed to move
                call get_command_argument(iarg+1,arg)
                read(arg,'(I15)') freeStart
                freeStart = freeStart+1
                iarg = iarg + 1
            endif
            if (trim(arg) == '-trSteps') then ! The number of steps for the translation of one bead
                call get_command_argument(iarg+1,arg)
                read(arg,'(I15)') trSteps
                iarg = iarg + 1
            endif
            if (trim(arg) == '-ci') then ! We are starting from a checkpoint
                    call get_command_argument(iarg+1,arg)
                    read(arg,'(A250)') checkInFile
                    restart=1 ! Set the flag that this is a restart
                    iarg = iarg+1
            endif
            if (trim(arg) == '-co') then ! Change the default checpoint outfile name
                    call get_command_argument(iarg+1,arg)
                    read(arg,'(A250)') checkOutFile
                    iarg = iarg+1
            endif
            if (trim(arg) == '-ForcePull') then ! 0 for no force pulling protocol (default), 1 to use
                    call get_command_argument(iarg+1,arg)
                    read(arg,*) ForcePulling
                    iarg = iarg + 1
            endif
            if (trim(arg) == '-memX') then ! Parameter for lipid thickness
                call get_command_argument(iarg+1,arg)
                read(arg,*) MBR ! sys_param.f90, forceField.f90
                ! Set dependent variables
                MBL = -MBR
                LEFTCUT = MBL-2.d0
                RIGHTCUT = MBR+2.d0 
                iarg = iarg + 1 ! Update this since the next argument is not a flag
            endif
            if (trim(arg) == '-memR') then ! Parameter for pore through lipid radius
                call get_command_argument(iarg+1,arg)
                read(arg,*) rChannel ! sys_param.f90, forceField.f90
                iarg = iarg + 1 ! Update this since the next argument is not a flag
            endif
            iarg = iarg + 1 ! Always update
        enddo
        write(6,*) 'These are the input files and options defined by the user:'
        write(6,*) 'Input coordinate file: ',trim(confInFile)
        write(6,*) 'Input polymer property file: ',trim(seqInFile)
        write(6,*) 'Input channel coordinate and type file: ',trim(chanInFile)
        write(6,*) 'Output coordinate file: ',trim(confOutFile)
        write(6,*) 'Trajectory file: ',trim(trajFile)
        write(6,*) 'Number of frames: ',nsteps
        write(6,*) 'Trajectory saved every ',write_frequency,' frames'
        write(6,*) 'Translation will take ',trSteps,' frames per bead.'
        write (6,*) 'This corresponds to a rate of ',3.d0/(trSteps*dtRealTime*1E-9),' res/sec'
        write (6,*) 'dgempty is',DGEMPTY,' and dGfillO is ',DGFILLO,' and dGfillC is ', DGFILLC
        write (6,*) 'FIXEDLG is set at', FIXEDLG

    end subroutine read_user_input

!==========================================================================================
! Subroutine init_io() loads the initial coordinates and pre-creates the output files. Please
! note that this overwrites all pre-excisting files with the same name!
!==========================================================================================

    subroutine init_io()

        ! Space to initialize variables used only here
        implicit none
        logical :: exists

        ! If this is a restart check if we are going to append the
        ! trajectory file
        if (restart/=0) then ! If there is a traj file we want to append to it
            inquire(file=trim(trajFile),exist=exists)
            if (exists) then
                open(unit=trajFID,file=trim(trajFile),status='old',position='append')
            else
                open(unit=trajFID,file=trim(trajFile),status='new')
            endif
            inquire(file='lg.txt',exist=exists)
            if (exists) then
                open(unit=77,file='lg.txt',status='old',position='append')
            else
                open(unit=77,file='lg.txt',status='new')
            endif
            inquire(file='pulling_lg.txt',exist=exists)
            if (exists) then
                open(unit=77,file='pulling_lg.txt',status='old',position='append')
            else
                open(unit=77,file='pulling_lg.txt',status='new')
            endif
            inquire(file=trim(pulling_trajFile),exist=exists)
            if (exists) then
                open(unit=pulling_trajFID,file=trim(pulling_trajFile),status='old',position='append')
            else
                open(unit=pulling_trajFID,file=trim(pulling_trajFile),status='new')
            endif
        else
            ! This is not a restart, if a traj file exists we want
            ! to replace it 
            open(unit=trajFID,file=trim(trajFile),status='replace')
            open(unit=pulling_trajFID,file=trim(pulling_trajFile),status='replace')
            open(unit=77,file='lg.txt',status='replace')
        endif
        close(trajFID)
        close(pulling_trajFID)
        close(77)
        ! Polymer files, restart and final configuration
        open(unit=polyFID,file=trim(confOutFile),status='replace')
        close(polyFID)
        ! TEMPORARY FILES FOR BUGTESTING
        open(unit=77,file='time.txt',status='replace') ! center of mass file
        close(77)
        open(unit=77,file='rog.txt',status='replace') ! radius of gyration file
        close(77)

    end subroutine init_io

!==========================================================================================
! Subroutine read_polymer(pos,dGtransfer,bLambda,charge,aname) reads in polymer info from a file.
!==========================================================================================

    subroutine read_polymer(pos,dG,cLamb,oLamb,charge,aname,mBead)

        implicit none
        double precision, intent(inout) :: pos(:,:), dG(:), cLamb(:), oLamb(:)
        integer, intent(inout) :: charge(:), mBead(:)
        character(3), intent(inout) :: aname(:)
        integer :: ib, err ! Used to iterate over beads
        character(1) :: dummy ! Gets rid of info I dont need

        ! Open the coordinate file for reading
        open(unit=77,file=trim(confInFile))
        ! First line contains the number of free polymer beads
        write(6,*) 'Reading polymer information; coordinate file contains ',nInit,' positions'
        read(77,*) dummy
        read(77,*)
        ! Read positions
        do ib=1,nInit
            read(77,*,iostat=err) dummy, pos(ib,x), pos(ib,y), pos(ib,z)
            if (err/=0) then
                write(6,*) 'There was an error reading polymer initial coordinates'
                stop
            endif
        enddo
        write(6,*) 'Successfully read polymer initial positions'
        close(77)

        ! Next get the properties of the polymer beads, there have to be
        ! properties defined for each bead!
        ! Open the file with properties
        open(unit=77,file=trim(seqInFile))
        ! Start reading
        do ib=1,nBeads
            read(77,*,iostat=err) aname(ib), dG(ib), cLamb(ib), oLamb(ib), charge(ib), mBead(ib)
            if (err/=0) then
                write(6,*) 'ERROR reading polymer bead properties for bead ',ib
                stop
            endif
            if (charge(ib).ne.0) then
                nCharge = nCharge +1
                chBeads(nCharge) = ib
            endif
        enddo
        write(6,*) 'Successfully read polymer bead properties from file'
        close(77)

    end subroutine read_polymer 

!==========================================================================================
! Subroutine read_channel(ctype,cpos,ceps,ccharge) reads in channel info from a file.
!==========================================================================================

    subroutine read_channel()

        implicit none
        integer :: cb, err, i ! Used to iterate over beads
        character(1) :: dummy

        ! Open the file for reading
        open(unit=77,file=trim(chanInFile))
        ! First line containts the number of channel beads
        read(77,*) nBeadChan
        write(6,*) 'Reading channel information; the channel contains ',nBeadChan,' particles'
        ALLOCATE(chanType(nBeadChan),chanBeads(nBeadChan,NDIM),chanP(nBeadChan,12), &
        chanCharge(nBeadChan),STAT=err)
        if (err/=0) then
            write(6,*) 'Could not allocate memory for channel bead properties'
            STOP
        endif
        !if (nBeadChan.gt.MAXNUMBEADS) then
        !        write(6,*) 'Channel has too many beads, adjust MAXNUMBEADS in sys_param.f90'
        !        STOP
        !endif
        read(77,*) dummy ! The comment line
        ! Read the positions and properties
        do cb=1,nBeadChan
            read(77,*,iostat=err) chanType(cb),chanBeads(cb,x),chanBeads(cb,y), &
            chanBeads(cb,z),chanCharge(cb), (chanP(cb,i),i=1,12)
            chanP(cb,7:10) = chanP(cb,7:10)**2 ! Cut-off squared for convenience
            chanP(cb,11:12) = chanP(cb,11:12)**2 ! Sigma squared for convenience
            if (err/=0) then
                write(6,*) 'ERROR reading channel info at ',cb
                stop
            endif
        enddo
        ! Set the cutoffs for the grid, only one grid-size is used for both electrostatic as well as NB
        xCutChanL = minval(chanBeads(1:nBeadChan,x))-max(ljRcrLG,dhRcrLG)
        xCutChanR = maxval(chanBeads(1:nBeadChan,x))+max(ljRcrLG,dhRcrLG)
        yCutChanL = minval(chanBeads(1:nBeadChan,y))-max(ljRcrLG,dhRcrLG)
        yCutChanR = maxval(chanBeads(1:nBeadChan,y))+max(ljRcrLG,dhRcrLG)
        zCutChanL = minval(chanBeads(1:nBeadChan,z))-max(ljRcrLG,dhRcrLG) 
        zCutChanR = maxval(chanBeads(1:nBeadChan,z))+max(ljRcrLG,dhRcrLG)
        write(6,*) 'GRID INFO:'
        write(6,*) xCutChanL, xCutChanR, yCutChanL, yCutChanR, zCutChanL, zCutChanR
        write(6,*) 'Successfully read the channel coordinates and properties'
        close(77)

    end subroutine read_channel

!==========================================================================================
! Subroutine read_checkpoint(pos) reads the info from a checkpoint file to restart
! a simulation. Currently I do NOT continue the quasi-random sequence of random
! numbers, I make a new seed and start anew.
!==========================================================================================

    subroutine read_checkpoint(pos,it,bip)
        !use rand_num
        integer :: dummy, i, err  ! use this to skip all the RNG lines
        double precision, intent(inout) :: pos(:,:), bip(:)
        integer(8), intent(inout) :: it ! iteration

        open(unit=77,file=trim(checkInFile))
        ! *** NEXT BLOCKS SKIPS THE RNG RESTART
        read(77,*) dummy
        do i=1,dummy+3 ! These are all the lines related to the RNG
                read(77,*)
        enddo
        ! *** UNCOMMENT THE NEXT, AND COMMENT THE PREVIOUS IF YOU WANT
        ! FULL RNG RESTART. ALSO UNCOMMENT use rand_num at the start
        !call init_random_seed_restart(77)
        ! *** The next lines read info on the polymer for restart
        read(77,*) it
        read(77,*) LGstatus
        read(77,*) nInit
        read(77,*) nForce
        do i=1,nBeads
            read(77,*,iostat=err) pos(i,x), pos(i,y), pos(i,z), bip(i)
            if (err/=0) then
                    write(6,*) 'ERROR: Incomplete coordinates in restart file'
                    STOP
            endif
        enddo
        close(77)
    end subroutine read_checkpoint

!==========================================================================================
! Subroutine write_checkpoint(pos) writes a checkpoint file that can be used to
! restart the simulation. Random number info is written as well as the
! coordinates of the polymer.
!==========================================================================================

    subroutine write_checkpoint(pos,it,bip)
        use rand_num 
        integer :: i, err
        double precision, intent(in) :: pos(:,:), bip(:)
        integer(8), intent(in) :: it ! iteration

        ! Random number stuff
        open(unit=77,file="checkpoint.temp",status='replace')
        write(77,'(I10)') seed_size
        do i=1,seed_size
                write(77,*) seed(i)
        enddo
        write(77,*) random_count
        write(77,*) iset
        write(77,*) gset
        ! Next polymer stuff
        write(77,'(I20)') it
        write(77,'(I1)') LGstatus
        write(77,'(I10)') nFree
        write(77,'(I10)') nForce
        do i=1,nBeads
                write(77,'(3F12.4)') pos(i,x), pos(i,y), pos(i,z), bip(i)
        enddo
!        write(77,*) ' '! Comment this out for normal usage
        close(77)
        ! Finally copy the temp file to overwrite the actual checkpoint
        command = "cp checkpoint.temp " // trim(checkOutFile)
        !command = "cat checkpoint.temp >> " // trim(checkOutFile)
        call system(command)

    end subroutine write_checkpoint

!==========================================================================================
! Subroutine traj_rerun(pos) reads the current positions from the trajectory file.
!==========================================================================================

    subroutine traj_rerun(pos)

        ! Space to initiliaze variables used only here
        implicit none
        double precision, intent(inout) :: pos(:,:)
        integer :: ib
        character(3) :: dummy ! Used to contain output I dont need

        ! Write necessary comment lines for .xyz format
        read(trajFID,*) dummy
        read(trajFID,*) dummy
        ! Write coordinates
        do ib = 1,nBeads
            read(trajFID,'(A3,6F12.4)') dummy,pos(ib,x),pos(ib,y),pos(ib,z) ! *** Dimensionality hardcoded
        enddo

    end subroutine traj_rerun

!==========================================================================================
! Subroutine traj_rerunF(pos,force) reads the current positions/forces from the trajectory file.
!==========================================================================================

    subroutine traj_rerunF(pos,force)

        ! Space to initiliaze variables used only here
        implicit none
        double precision, intent(inout) :: pos(:,:), force(:,:)
        integer :: ib
        character(3) :: dummy ! Used to contain output I dont need

        ! Write necessary comment lines for .xyz format
        read(trajFID,*) dummy
        read(trajFID,*) dummy
        ! Write coordinates
        do ib = 1,nBeads
            read(trajFID,'(A3,6F12.4)') dummy,pos(ib,x),pos(ib,y),pos(ib,z),force(ib,x),force(ib,y),force(ib,z) ! *** Dimensionality hardcoded
        enddo

    end subroutine traj_rerunF

!==========================================================================================
! Subroutine write_trajF(pos,force,aname) writes the current positions/forces to the trajectory file.
!==========================================================================================

    subroutine write_trajF(pos,force,aname)

        ! Space to initiliaze variables used only here
        implicit none
        double precision, intent(in) :: pos(:,:), force(:,:)
        character(3), intent(in) :: aname(:)
        integer :: ib

        ! Open the file for writing
        open(unit=trajFID,file=trim(trajFile),position='append')
        ! Write necessary comment lines for .xyz format
        write(trajFID,*) nBeads
        write(trajFID,'(A11)') 'Commentline'
        ! Write coordinates
        do ib = 1,nBeads
            write(trajFID,'(A3,6F12.4)') aname(ib),pos(ib,x),pos(ib,y),pos(ib,z),force(ib,x),force(ib,y),force(ib,z) ! *** Dimensionality hardcoded
        enddo
        close(trajFID) ! Close the file again

    end subroutine write_trajF


    subroutine write_pulling_trajF(pos,force,aname)

        ! Space to initiliaze variables used only here
        implicit none
        double precision, intent(in) :: pos(:,:), force(:,:)
        character(3), intent(in) :: aname(:)
        integer :: ib

        ! Open the file for writing
        open(unit=pulling_trajFID,file=trim(pulling_trajFile),position='append')
        ! Write necessary comment lines for .xyz format
        write(pulling_trajFID,*) nBeads
        write(pulling_trajFID,'(A11)') 'Commentline'
        ! Write coordinates
        do ib = 1,nBeads
            write(pulling_trajFID,'(A3,6F12.4)') aname(ib),pos(ib,x),pos(ib,y),pos(ib,z),force(ib,x),force(ib,y),force(ib,z) ! *** Dimensionality hardcoded
        enddo
        close(pulling_trajFID) ! Close the file again

    end subroutine write_pulling_trajF

!==========================================================================================
! Subroutine write_force(force) writes the current forces to the debug file.
!==========================================================================================

        subroutine write_force(f)
                implicit none
                double precision, intent(in) :: f(:)
                write(debugFID,*) f
        end subroutine write_force

!==========================================================================================
! Subroutine write_trajD(pos,force,aname) writes the current positions/forces to
! the debug file.
!==========================================================================================

    subroutine write_trajD(pos,force,aname)

        ! Space to initiliaze variables used only here
        implicit none
        double precision, intent(in) :: pos(:,:), force(:,:)
        character(3), intent(in) :: aname(:)
        integer :: ib

        ! Open the file for writing
       ! open(unit=trajFID,file=trim(trajFile),position='append')
        ! Write necessary comment lines for .xyz format
        write(trajDFID,*) nBeads
        write(trajDFID,'(A11)') 'Commentline'
        ! Write coordinates
        do ib = 1,nBeads
            write(trajDFID,'(A3,6F12.4)') aname(ib),pos(ib,x),pos(ib,y),pos(ib,z),force(ib,x),force(ib,y),force(ib,z)
!*** Dimensionality hardcoded
        enddo
!       close(trajFID) ! Close the file again

    end subroutine write_trajD

!==========================================================================================
! Subroutine write_traj(pos,aname) writes the current positions to the trajectory file.
!==========================================================================================

    subroutine write_traj(pos,aname)

        ! Space to initiliaze variables used only here
        implicit none
        double precision, intent(in) :: pos(:,:)
        character(3), intent(in) :: aname(:)
        integer :: ib

        ! Open the file for writing
        open(unit=trajFID,file=trim(trajFile),position='append')
        ! Write necessary comment lines for .xyz format
        write(trajFID,*) nBeads
        write(trajFID,'(A11)') 'Commentline'
        ! Write coordinates
        do ib = 1,nBeads
            write(trajFID,'(A3,3F12.4)') aname(ib),pos(ib,x),pos(ib,y),pos(ib,z) ! *** Dimensionality hardcoded
        enddo
        close(trajFID) ! Close the file again

    end subroutine write_traj

!==========================================================================================
! Subroutine write_LGstatus() writes the current state of the lateral gate and its FE of opening 
! used in testing LG dynamics
!==========================================================================================

    subroutine write_LGstatus()

        implicit none

        ! Write to file
        open(unit=77,file='lg.txt',position='append')
        write(77,'(I1)') LGstatus ! Save the state of the LG
        close(77)

    end subroutine write_LGstatus

    subroutine write_pulling_LGstatus()

        implicit none

        ! Write to file
        open(unit=77,file='pulling_lg.txt',position='append')
        write(77,'(I1)') LGstatus ! Save the state of the LG
        close(77)

    end subroutine write_pulling_LGstatus

!==========================================================================================
! Subroutine write_com(pos) writes the position of the center of mass to com.txt, I used 
! this subroutine mainly to ensure that the nascent chain dynamics are ok.
!==========================================================================================

    subroutine write_com(pos)

        implicit none
        double precision, intent(in) :: pos(:,:) ! position array
        double precision :: com(NDIM) ! coords for COM

        ! Determine COM coordinate
        com = sum(pos(:,:),dim=1)/nBeads
        ! Write to file
        open(unit=77,file='com.txt',position='append')
        write(77,'(3F12.4)') com(1),com(2),com(3) ! *** ASSUMING 3D
        close(77)

    end subroutine write_com

!==========================================================================================
! Subroutine write_frameN(it) writes the current frame number to time.txt, I used this
! subroutine to check when my simulations complete
!==========================================================================================

    subroutine write_frameN(it)

        implicit none
        integer, intent(in) :: it ! Iteration

        ! Write to file
        open(unit=77,file='time.txt',position='append')
        write(77,'(I10)') it
        close(77)

    end subroutine write_frameN

!==========================================================================================
! Subroutine write_rg(pos) writes the radius of gyration squared to rog.txt, I used this
! subroutine to ensure that the nascent chain dynamics are ok.
!==========================================================================================

    subroutine write_rg(pos)

        implicit none
        double precision, intent(in) :: pos(:,:) ! position array
        double precision :: rog2 ! radius of gyration squared

        rog2 = sum((pos(nBeads,:)-pos(1,:))**2)
        ! Write to file
        open(unit=77,file='rog.txt',position='append')
        write(77,'(F12.4)') rog2
        close(77)

    end subroutine write_rg

end module iolib
