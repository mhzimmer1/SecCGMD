!====================================================================================================
! Contains subroutines used for analysis of my simulations.
!====================================================================================================

module analysis

    use sys_param ! All parameters come from this module
    use iolib, only : nsteps ! Number of iterations

    implicit none
    ! Space for defining variables known throughout the module
    integer :: bhist(201) ! Contain bond data, I use the fact here that the max bondL is determined
                ! feneR2
    double precision, parameter :: integrateD2=100.0 ! TM beads at least 10 sigma from channel center
    double precision :: fixForce(NDIM)
    integer(8) :: nfbin=20001, fhist(20001)
    double precision :: bs=0.01d0, f0=-100.d0

    contains

!====================================================================================================
! Subroutine init_analysis() is used to initialize the analysis module. Uncomment what should be
! analyzed.
!====================================================================================================

    subroutine init_analysis()

        implicit none
        ! Local variables

        ! Initialize for bondL histogram
!       bhist = 0 ! zerorize
!       open(unit=77,file='bhist.txt',status='replace')
!       close(77)
                ! This file will be used in do_analysis(pos,force) to write
                ! numerous properties during runtime
        open(unit=analysisFID,file='analysis.txt',status='replace')
        open(unit=energyFID,file='energies.txt',status='replace')
! *** DEBUG FILES, UNCOMMENT WHEN REQUIRED
        ! open(unit=debugFID,file='debug.txt',status='replace')
        ! open(unit=trajDFID,file='debugTraj.xyz',status='replace')
        fixForce(:) = 0.d0

    end subroutine init_analysis

!====================================================================================================
! Subroutine peptide_status(pos,flag) used to asses if the nascent chain reached a terminal topology.
!====================================================================================================

subroutine peptide_status(r,pstat)

        use polymer, only : membBead, nMemb
        implicit none
        double precision, intent(in) :: r(:,:) ! Nascent chain position
        integer, intent(inout) :: pstat ! status flag, 0 if no exit status
        integer :: ib, nb, channel_empty
        double precision :: rho2

        ! Check if the entire chain has been translated
        pstat = 1
        channel_empty = 1
                ! *** I altered the index here, I start checking for
                ! translocated beads at the putative TMD. This way I avoid
                ! issues with including the anchor, and I can use a greater
                ! distance from the membrane. ***
                ! membBead(1) = 1 if there is no TMD
        do ib = membBead(1),nBeads ! It is vital to not include constrained beads!
            if (r(ib,1).lt.2.5) pstat = 0 ! This bead is not translocated
        enddo
        ! RVL EDITS: Terminates successfully if TMD has 1 of two orientations
        ! and either A) the TMD is laterally > integrateD2 away, or B) the LG is
        ! closed, and no beads are in the channel. 
        if ((LGstatus.eq.0 .or. FIXEDLG.eq.1).and.(pstat.eq.0)) then ! If its not translocated, check if integrated
            if ((r(membBead(1)-1,1).lt.MBL).and.(r(membBead(nMemb)+2,1).gt.MBR)) pstat=2 ! Ncyt/Cer
            if ((r(membBead(1)-1,1).gt.MBR).and.(r(membBead(nMemb)+2,1).lt.MBL)) pstat=3 ! Ner/Ccyt
            ! Correct orientations for TM segment; find reason to disqualify now
            if (pstat/=0) then
                ! check distance of TMD beads radially from center of channel.
                do ib = 1,nBeads
                   rho2=sum(r(ib,2:3)**2)
                   ! check if in channel, defined as < 2.0 radially and >-2, <2
                   ! axially
                   if (abs(r(ib,1)) .lt. 2.0 .and. rho2 .lt. 4.0) then
                      ! channel is not empty, so set flag
                      channel_empty = 0
                   endif
                enddo
                do ib = 1,nMemb
                    nb = membBead(ib)
                    rho2=sum(r(nb,2:3)**2)
                    if (abs(r(nb,1)).gt.2.0) pstat=0 ! Not in membrane
                    ! Not finished if rho2 < integrateD2 and the channel is
                    ! empty; if rho2 > integrateD2 OR the channel is empty we
                    ! can end
                    if (rho2.lt.integrateD2 .and. channel_empty .eq. 0 ) pstat=0 ! Not far enough
                enddo
            endif
        endif

    end subroutine peptide_status

!====================================================================================================
! Subroutine do_analysis(pos,force) does a bunch of analysis as specified here. If you want to 
! change the analysis that is done just change what is called here.
!====================================================================================================

    subroutine do_analysis(pos,force,ftmp)
        use forceField ! For force calculations
        implicit none
        double precision, intent(in) :: pos(:,:), force(:,:) ! Positions and forces
        double precision, intent(out) :: ftmp(NDIM)
        double precision :: dr(NDIM), r2, ri2 ! for force calculation
        integer :: bi

! *** ENERGY CALCULATION MOVED TO MAIN.f90 ***
!                double precision :: energies(6) ! Used to write the energy out
!                integer :: i ! Used to loop

!       call calc_bondL(pos,bhist) ! Add bondL to histogram
!                tgmd = mindist(pos(1:75,:),pos(78:153,:)) ! minimum distance between TatC and GFP
!                comv = comvector(pos(1:75,:),pos(78:153,:)) ! Vector between center of mass, from tatC to GFP
!                tolink = comvector(pos(1:75,:),pos(75:76,:)) ! Vector between TatC COM and the linker
!                tatx = sum(pos(1:75,1))/75.d0 ! The average X coord of TatC                

                ! Write what I want to analyze to the analysis file
        ftmp(:) = 0.d0; 
        dr(:) = pos(nForce,:) - pos(nForce-1,:); 
        r2 = sum(dr(:)**2);
        ri2 = 1.d0/r2;
        ftmp(:) = ftmp(:) - dr(:)*feneK/(1.d0-r2/feneR2) ! Attractive part
        if (r2.lt.ljRcr2) then
            ftmp(:) = ftmp(:) + 48.d0*dr(:)*ri2**4*ljEps*(ri2**3-0.5d0) ! Repulsive
        endif
        fixForce(:) = fixForce(:) + ftmp
        ! Write what I want to analyze to the analysis file
        ! Add the current force to the histogram of forces
        bi = floor((ftmp(1)-f0)/bs)
        if (bi.lt.1) bi=1
        if (bi.gt.nfbin) bi = nfbin
        fhist(bi) = fhist(bi) + 1       

    end subroutine do_analysis

!====================================================================================================
! Subroutine force_check(pos,force,aname) checks the particle that is feeling the
! highest force, print the particle ID, positions -1:1:+1 and the force. This is
! meant to be done every iteration for simulations that are experiencing trouble
!====================================================================================================

        subroutine force_check(pos,force,aname)

                use iolib
                use forceField
!                use scaffold

                implicit none
                double precision, intent(in) :: pos(:,:), force(:,:) ! Position, force
                character(3), intent(in) :: aname(:) ! names
                double precision :: fdum(MAXNUMBEADS,NDIM) ! dummy
                double precision :: maxF, ftmp
                integer :: ib, idmaxF ! For loop, index of highest force bead

                ! Initialize
                maxF = 0.d0 ! Small initially
                ! Find the max force
                do ib = 1,nForce
                        ftmp = sum(force(ib,:)**2)
                        if (ftmp.gt.maxF) then
                                maxF = ftmp
                                idmaxF = ib
                        endif
                enddo
                maxF = sqrt(maxF) ! Take root now
                ! Write to the debug file
                !write(debugFID,*) maxF, idmaxF
                ! In the case where the maximum force is larger than a certain
                ! value, set a flag to force a full position/force dump
                if (maxF.gt.1000.d0) then 
                        CALL write_trajF(pos,force,aname)
                        CALL write_trajD(pos,force,aname)
                        fdum(:,:) = 0.d0
                        CALL bonded_forces(pos,fdum)
                        ftmp = sqrt(sum(fdum(idmaxF,:)**2))
                        ! write(debugFID,*) 'Bonded force: ',ftmp
                        fdum(:,:) = 0.d0
                        CALL nb_forces(pos,fdum)
                        ftmp = sqrt(sum(fdum(idmaxF,:)**2))
                        ! write(debugFID,*) 'Non-bonded force: ',ftmp
                        fdum(:,:) = 0.d0
 !                       CALL scaffoldF(pos,fdum)
 !                       ftmp = sqrt(sum(fdum(idmaxF,:)**2))
 !                       write(debugFID,*) 'Scaffold force: ',ftmp
                endif
                if (maxF.gt.1000.d0) then
                        write(6,*) 'ERROR: forces are too large'
                        STOP 
                endif

        end subroutine force_check

!====================================================================================================
! Subroutine write_FE_axis() can be used to write the FE profile for translocation through the  
! channel. Modify the subroutine if the profile for a different channel is desired (or along another
! axis). 
!====================================================================================================

    subroutine write_FE_axis()

                use forceField, only : table_ljC, table_ljO, table_ljCpho, table_ljOpho, hNGRIDX, hNGRIDY, hNGRIDZ
        implicit none
        integer :: xb, yb, zb ! Bin IDs
                double precision :: x,y,z,r ! coordinates
                double precision :: FE(NGRIDX,4) ! Used to store the FE values
                double precision :: ftmpC, ftmpO, ftmpCP, ftmpOP, u ! Temporary store FE, PE

                FE(:,:) = 0.0 ! Initialize
        do xb = 1,NGRIDX ! For formatting so that Matlab can read it easily
                        x = (xb-hNGRIDX)*DGRID
                        ftmpC = 0.0; ftmpO = 0.0; ftmpCP = 0.0; ftmpOP = 0.0;
                        do yb = 1,NGRIDY
                                y = (yb-hNGRIDY)*DGRID
                                do zb = 1,NGRIDZ
                                        z = (zb-hNGRIDZ)*DGRID
                                        r = sqrt(z**2 + y**2)
                                        if (r.le.2.5) then ! Its in the cylinder
                                                u = table_ljC(xb,yb,zb)
                                                ftmpC = ftmpC + exp(-beta*u)
                                                u = 0.7d0*table_ljCpho(xb,yb,zb) &
                                                + 0.3d0*table_ljC(xb,yb,zb)
                                                ftmpO = ftmpO + exp(-beta*u)
                                                u = table_ljCpho(xb,yb,zb) 
                                                ftmpCP = ftmpCP + exp(-beta*u)
                                                u = 0.3d0*table_ljCpho(xb,yb,zb) &
                                                + 0.7d0*table_ljC(xb,yb,zb)
                                                ftmpOP = ftmpOP + exp(-beta*u)
                                        endif                                                
                                enddo
            enddo
                        FE(xb,1) = ftmpC ! Store the boltzmann weight
                        FE(xb,2) = ftmpO ! Store the boltzmann weight
                        FE(xb,3) = ftmpCP ! Store the boltzmann weight
                        FE(xb,4) = ftmpOP ! Store the boltzmann weight
        enddo
                ! Normalize
                FE(:,1) = FE(:,1)/FE(1,1) ! This bin should be free solvent
                FE(:,2) = FE(:,2)/FE(1,2) ! This bin should be free solvent
                FE(:,3) = FE(:,3)/FE(1,3) ! This bin should be free solvent
                FE(:,4) = FE(:,4)/FE(1,4) ! This bin should be free solvent
        ! Open the files for writing
        open(unit=81, file='translocation_FE_closed.txt', status='replace')
        open(unit=82, file='translocation_FE3phil7pho_closed.txt', status='replace')
        open(unit=83, file='translocation_FEpho_closed.txt', status='replace')
        open(unit=84, file='translocation_FE7phil3pho_closed.txt', status='replace')
                do xb = 1,NGRIDX
                        write(81,*) -kBT*log(FE(xb,1))
                        write(82,*) -kBT*log(FE(xb,2))
                        write(83,*) -kBT*log(FE(xb,3))
                        write(84,*) -kBT*log(FE(xb,4))
                enddo
                close(81); close(82); close(83); close(84);

    end subroutine write_FE_axis

!====================================================================================================
! Subroutine write_channel(sr) can be used to write the current channel PE landscape to a 
! file for open and closed channel; in 3-dimensions (takes alot of space). Since the grid tends to 
! be very fine it can be resampled taking only every sr points.
!====================================================================================================

    subroutine write_channel(sr)

                use forceField, only : table_ljC, table_ljO, table_chC, table_chO
        implicit none
        integer, intent(in) :: sr ! Every sr gridpoints are written to the files
        integer :: xb, yb, zb ! Bin IDs

        ! Open the files for writing
        open(unit=81, file='channel_chC.txt', status='replace')
        open(unit=82, file='channel_ljC.txt', status='replace')
        open(unit=83, file='channel_chO.txt', status='replace')
        open(unit=84, file='channel_ljO.txt', status='replace')
        do xb = 1,NGRIDX,sr ! For formatting so that Matlab can read it easily
            do yb = 1,NGRIDY,sr
                                do zb = 1,NGRIDZ,sr
                        write(82,"(F12.6)") table_ljC(xb,yb,zb)
                        write(84,"(F12.6)") table_ljO(xb,yb,zb)
                        write(81,"(F12.6)") table_chC(xb,yb,zb)
                        write(83,"(F12.6)") table_chO(xb,yb,zb)
                                enddo
            enddo
        enddo
        close(81); close(82); close(83); close(84);

    end subroutine write_channel

!====================================================================================================
! Subroutine write_channelXY(sr) can be used to write the current channel PE landscape to a 
! file for open and closed channel; in the XY plane. Since the grid tends to be very fine it 
! can be resampled taking only every sr points.
!====================================================================================================

    subroutine write_channelXY(sr)

                use forceField, only : table_ljC, table_ljO, table_chC, table_chO
        implicit none
        integer, intent(in) :: sr ! Every sr gridpoints are written to the files
        integer :: xb, yb, zb ! Bin IDs

        ! Open the files for writing
        open(unit=81, file='channel_chC_xy.txt', status='replace')
        open(unit=82, file='channel_ljC_xy.txt', status='replace')
        open(unit=83, file='channel_chO_xy.txt', status='replace')
        open(unit=84, file='channel_ljO_xy.txt', status='replace')
        do xb = 1,NGRIDX,sr ! For formatting so that Matlab can read it easily
            do yb = 1,NGRIDY,sr
                write(82,"(F12.6)") table_ljC(xb,yb,NGRIDZ/2)
                write(84,"(F12.6)") table_ljO(xb,yb,NGRIDZ/2)
                write(81,"(F12.6)") table_chC(xb,yb,NGRIDZ/2)
                write(83,"(F12.6)") table_chO(xb,yb,NGRIDZ/2)
            enddo
        enddo
        close(81); close(82); close(83); close(84);

    end subroutine write_channelXY

!====================================================================================================
! Subroutine write_channelXZ(sr) can be used to write the current channel PE landscape to a 
! file for open and closed channel; in the XZ plane. Since the grid tends to be very fine it 
! can be resampled taking only every sr points.
!====================================================================================================

    subroutine write_channelXZ(sr)

                use forceField, only : table_ljC, table_ljO, table_chC, table_chO
        implicit none
        integer, intent(in) :: sr ! Every sr gridpoints are written to the files
        integer :: xb, yb, zb ! Bin IDs

        ! Open the files for writing
        open(unit=81, file='channel_chC_xz.txt', status='replace')
        open(unit=82, file='channel_ljC_xz.txt', status='replace')
        open(unit=83, file='channel_chO_xz.txt', status='replace')
        open(unit=84, file='channel_ljO_xz.txt', status='replace')
        do zb = 1,NGRIDZ,sr ! For formatting so that Matlab can read it easily
            do xb = 1,NGRIDX,sr
                write(82,"(F12.6)") table_ljC(xb,NGRIDY/2,zb)
                write(84,"(F12.6)") table_ljO(xb,NGRIDY/2,zb)
                write(81,"(F12.6)") table_chC(xb,NGRIDY/2,zb)
                write(83,"(F12.6)") table_chO(xb,NGRIDY/2,zb)
            enddo
        enddo
        close(81); close(82); close(83); close(84);

    end subroutine write_channelXZ

!====================================================================================================
! Subroutine write_channelYZ(sr) can be used to write the current channel PE landscape to a 
! file for open and closed channel; in the YZ plane. Since the grid tends to be very fine it 
! can be resampled taking only every sr points.
!====================================================================================================

    subroutine write_channelYZ(sr)

                use forceField, only : table_ljC, table_ljO, table_chC, table_chO
        implicit none
        integer, intent(in) :: sr ! Every sr gridpoints are written to the files
        integer :: xb, yb, zb ! Bin IDs

        ! Open the files for writing
        open(unit=81, file='channel_chC_yz.txt', status='replace')
        open(unit=82, file='channel_ljC_yz.txt', status='replace')
        open(unit=83, file='channel_chO_yz.txt', status='replace')
        open(unit=84, file='channel_ljO_yz.txt', status='replace')
        do zb = 1,NGRIDZ,sr ! For formatting so that Matlab can read it easily
            do yb = 1,NGRIDY,sr
                write(82,"(F12.6)") table_ljC(NGRIDX/2,yb,zb)
                write(84,"(F12.6)") table_ljO(NGRIDX/2,yb,zb)
                write(81,"(F12.6)") table_chC(NGRIDX/2,yb,zb)
                write(83,"(F12.6)") table_chO(NGRIDX/2,yb,zb)
            enddo
        enddo
        close(81); close(82); close(83); close(84);

    end subroutine write_channelYZ

!====================================================================================================
! Subroutine calc_bondL(pos,bhist) can be used to make a histrogram of the average bond length. The 
! first bond in the sequence is picked as a reference bond.
!====================================================================================================

    subroutine calc_bondL(pos,bhist)

        implicit none
        double precision, intent(in) :: pos(:,:) ! Array with particle positions
        integer, intent(inout) :: bhist(:) ! Histogram with bondLength values
        double precision :: b ! bondlength

        b = 100.0*sqrt(sum((pos(2,:)-pos(1,:))**2)) ! multiply by 100 for resolution ***
        bhist(int(b)+1) = bhist(int(b)+1) + 1

    end subroutine calc_bondL

!====================================================================================================
! Subroutine write_analysis() writes analysis data to output files
!====================================================================================================

    subroutine write_analysis(exitflag,iteration)

        implicit none
        integer, intent(in) :: exitflag 
                integer(8), intent(in) :: iteration
        integer :: i ! Used for loops

        ! Write out bhist
!       open(unit=77,file='bhist.txt',position='append')
!       do i = 1,201
!           write(77,*) bhist(i)
!       enddo
!       close(77)
        open(unit=77,file='status.txt',status='replace')
        write(77,'(2I20)') exitflag, iteration
        close(77)

    ! *** Added for force pulling
        if (ForcePulling .eq. 1) then
            do i=1,nfbin
                write(analysisFID,*) f0+i*bs, fhist(i)
            enddo
        endif


        close(energyFID) ! Close the energy file
        close(analysisFID) ! Close the main analysis file
        ! close(debugFID) ! Close the debug file

    end subroutine write_analysis

!====================================================================================================
! **** FUNCTIONS GO HERE ****
!====================================================================================================

!====================================================================================================
! function mindist(g1,g2) calculates and returns the minimum distance between two groups of particles.
!====================================================================================================

double precision function mindist(g1,g2)

        double precision, intent(in) :: g1(:,:), g2(:,:)
        double precision :: d2, md2 ! Squared (minimum) distance
        integer :: ib, jb ! For do loops

        md2 = 1000.d0 ! Initialize this large
        do ib = 1,size(g1,dim=1)
                do jb = 1,size(g2,dim=1)
                        d2 = sum((g1(ib,:)-g2(jb,:))**2) ! distance squared
                        if (d2.lt.md2) md2 = d2 ! re-assign
                enddo
        enddo
        mindist = sqrt(md2) ! Finally take the root

end function mindist

!====================================================================================================
! function comvector(g1,g2) calculates and returns the vector betweem the COM of two groups of particles.
!====================================================================================================

function comvector(g1,g2)

        double precision :: comvector(NDIM) ! Its an array
        double precision, intent(in) :: g1(:,:), g2(:,:)
        double precision :: com1(NDIM), com2(NDIM) ! Separate COM
        integer :: ib ! For loop

        com1(:) = sum(g1(:,:),dim=1)/size(g1,dim=1)
        com2(:) = sum(g2(:,:),dim=1)/size(g2,dim=1)

        comvector(:) = com2(:) - com1(:)

end function comvector

end module analysis
