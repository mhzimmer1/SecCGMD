!=========================================================================================
! This module contains subroutines and variables related to the state of the polymer.
!=========================================================================================

module polymer
    
    use sys_param   ! From this we get the parameters used

    ! First define the variables used in this module
    implicit none
    character(3) :: aname(MAXNUMBEADS) ! Will contain the name for each atom
    double precision :: pos(MAXNUMBEADS,NDIM), force(MAXNUMBEADS,NDIM) ! Position and force array
    double precision :: dGtransfer(MAXNUMBEADS), cLambda(MAXNUMBEADS), oLambda(MAXNUMBEADS) ! Contains the transfer free energy for each bead
    integer :: bCharge(MAXNUMBEADS), assignTM(MAXNUMBEADS) ! charge and TM region
    double precision :: transSide(MAXNUMBEADS), bipBinding(MAXNUMBEADS) ! BiP Binding
    integer, allocatable :: membBead(:) 
    double precision :: pairwise(MAXNUMBEADS,MAXNUMBEADS) ! Used to define TMDs and specific pairwise interactions
    double precision :: intracellular(MAXNUMBEADS), extracellular(MAXNUMBEADS) ! Used to define localization restraints
    integer :: nMemb

    contains
    ! The sub routines go here
    
!=========================================================================================
! Subroutine to initialize the polymer, the initial configuration is read from a user
! input file.
!=========================================================================================

    subroutine init_polymer()

        use iolib, only : read_polymer ! Subroutine for reading polymer information
        implicit none
        ! Defaults, will be used for beads not contained in the file
        pos(:,x) = exitX; pos(:,y) = exitY; pos(:,z) = exitZ; ! All beads at ribosome exit
        dGtransfer(:) = 0.d0 ! No phobicity
        cLambda(:) = 0.d0 ! No phobicity
        oLambda(:) = 0.d0 ! No phobicity
        transSide(:) = 0.d0 ! Bead is translocated
        bipBinding(:) = 0.d0 ! Bead is bound to BiP
        aname(:) = 'C' ! This name only affects the color of the beads
        assignTM(:)=-1
        write(6,*) 'Reading polymer info'
        call read_polymer(pos,dGtransfer,cLambda,oLambda, bCharge,aname, assignTM) ! Read info from file
        write(6,*) 'Read polymer info, defining TMDs'
        ! *** Added for specific NC-NC interactions
!        write(6,*) 'Going to allocate an array for pairwise interactions'
!        allocate(pairwise(nBeads,nBeads)) ! Used for specific pairwise interactions
!        write(6,*) 'Successfully allocated the pairwise array'
        pairwise(:,:) = 0.d0 ! All repulsive
        write(6,*) 'Zerorized the pairwise array'
        intracellular(:) = 0.d0 ! No restraint forces
        extracellular(:) = 0.d0 
        write(6,*) 'Zerorized the localization restraint arrays'
        call define_memb_beads()
        write(6,*) 'Defined TMDs'
        call define_pairwise()
        write(6,*) 'Defined pairwise interactions'
        call define_localization()
        write(6,*) 'Defined localization restraints'
        ! If the file does not contain information I use standard values.
        if (nBeads.gt.MAXNUMBEADS) then
            write(6,*) 'Polymer has too many beads, adjust MAXNUMBEADS in sys_param.f90'
            STOP
        endif

    end subroutine init_polymer

!==========================================================================================
! Subroutine define_localization, defines localization restraints on
! user-defined beads. Updates entries in the intracellular and extracellular
! arrays with the restraint force on the bead
!==========================================================================================

    subroutine define_localization()

        use iolib, only : localizationFile
        implicit none
        integer :: ib, n, nres ! For loop
        integer :: err, err_open, err_read ! For file-reading errors
        double precision :: resForceI, resForceE ! Restraint forces

        err_read = 1
        open(unit=77,file=trim(localizationFile),iostat=err)
        if (err .eq. 0) then
            read(77,*, iostat=err_read) nres
        endif

        if (err_read .eq. 0) then
            read(77,*) nres ! First line has the number of restraints
            do n=1,nres
                read(77,*,iostat=err) ib, resForceI, resForceE ! bead, intracellular restraint, extracellular restraint
                if (err/=0) then
                    write(6,*) 'ERROR: unable to read localization restraint file'
                    STOP
                endif
                intracellular(ib) = resForceI
                extracellular(ib) = resForceE
                write(6,*) 'Bead ',ib,' location restraint with constant forces: ', resForceI,' and ',resForceE
            enddo
            write(6,*) 'Successfully read localization restraints from file ',localizationFile
        else
            write(6,*) 'ERROR: unable to read localization restraint file'
        endif
        close(77)

    end subroutine define_localization

!==========================================================================================
! Subroutine define_pairwise, defines the pairwise interactions between beads,
! as set by the user in an optional input file. This subroutine updates entries
! in the pairwise(:,:) array to indicate attractive interactions between beads.
!==========================================================================================

    subroutine define_pairwise()

        use iolib, only : pairwiseFile
        implicit none
        integer :: npairs, n, ib, jb ! For loop
        integer :: err, err_open, err_read ! Keep track of file-reading errors
        double precision :: eps_ib_jb

        err_read = 1
        open(unit=77,file=trim(pairwiseFile),iostat=err_open)
        if (err_open .eq. 0) then
            read(77,*, iostat=err_read) npairs ! First line contains the number of pairs
        endif

        if(err_read .eq. 0) then            
            do n=1,npairs
                read(77,*,iostat=err) ib, jb, eps_ib_jb ! bead1, bead2, epsilon
                if (err/=0) then
                    write(6,*) 'ERROR: unable to read pairwise interaction file'
                    STOP
                endif
                pairwise(ib,jb) = eps_ib_jb
                pairwise(jb,ib) = eps_ib_jb ! In case user provides the higher index bead first
                write(6,*) 'Pairwise interaction between ',ib,' and ',jb,' with epsilon ',eps_ib_jb
            enddo
            write(6,*) 'Successfully read pairwise interactions from file ',pairwiseFile
        else
            write(6,*) 'No pairwise interactions specified by user'
        endif
        close(77)

    end subroutine define_pairwise

!==========================================================================================
! Subroutine to identify which set of beads will define the TM region
! Assumes only one, discrete group.
!==========================================================================================

    subroutine define_memb_beads()
        implicit none
        integer :: ib, jb, curBead

        do ib=1, nBeads
            ! identify current loop / TM and increment bead count accordingly
            if (assignTM(ib) .gt. 0) then
                nMemb= nMemb + 1
            endif
        enddo

        print *, "Beads assigned to TM ", nMemb

        ! allocate array
        allocate(membBead(nMemb))

        ! now save relevant indices
        curBead = 1
        do ib = 1, nBeads
            if (assignTM(ib) .gt. 0) then
                membBead(curBead) = ib
                curBead = curBead + 1
            endif
        end do

    end subroutine define_memb_beads


!=========================================================================================
! Subroutine to initialize the polymer, in this case specifically for free polymer
! translocation simulations. These simulations were used to determine the diffusion coefficient.
!=========================================================================================

    subroutine init_polymer_diffusion()

        implicit none
        integer :: ib       ! Used to iterate over the beads

        ! *** Here I make sure the last bead is inside the channel at (0.5,0)
        do ib = 1,nBeads
            pos(ib,x) = ib*1.0  ! purely as a test case we start with a
                        ! stretched chain, initial distances are sigma
            pos(ib,y) = 0.0
            pos(ib,z) = 0.0
            dGtransfer(ib) = -2.0 ! *** hydrophilic polypeptide
        end do
        pos(:,x) = pos(:,x) - nBeads + 0.5 ! *** Place the last bead inside the channel

        ! Set the initial forces to zero
        force = 0.d0

    end subroutine init_polymer_diffusion

end module polymer
