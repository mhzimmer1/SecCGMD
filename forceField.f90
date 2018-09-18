!=======================================================================================
! Module containing subroutines related to calculating the forces on the polymer beads
!=======================================================================================

module forceField

    use sys_param ! Load the parameters required
    use polymer, only : dGtransfer, cLambda, oLambda, bCharge, transSide, bipBinding, pairwise, intracellular, extracellular
                         ! Hydrophobicity, channel interaction, charge, translocation state, bip binding, NC-NC, localization

    ! Define variables that are used in this module
    implicit none
    ! Pre-defined for more efficient calculation
    double precision :: epsCR=0.0, epsCL=0.0, epsCLR=0.0, epsCRR=0.0 
    ! Pre-defined for excluded volume TO REMOVE 
    !double precision :: repulEcut =0.0, alphaRepul=0.0,  repulR2 = 0.0
! The LJ energy at the right cut-off, variable to store excess output
    double precision :: dhCL=0.0, dhCR=0.0 ! Debye-Huckel potential values at cut-offs
    ! Related to the tabulated tanh function used in the solvation potential
    double precision :: tabulated_tanh(NGRID), table_stepsize
    ! The following are used in force table lookup
    double precision :: hNGRIDX, hNGRIDY, hNGRIDZ ! Used in table lookup
    ! Related to the tabulated forces for the closed translocon
    double precision, allocatable :: table_ljC(:,:,:), & ! PE and force in the same array for efficient indexing
                                table_ljO(:,:,:), & ! Open channel
                                table_ljCpho(:,:,:), & ! Closed channel for phobic beads
                                table_ljOpho(:,:,:), & ! Open channel for phobic beads
                                table_chC(:,:,:), & ! Closed channel charges
                                table_chO(:,:,:)  ! Open channel charges
    ! Related to Bip binding
    double precision :: bip_on, bip_off ! Effective on and off probabilities (rate*time)

    contains

!=======================================================================================
! Subroutine init_ff() sets some of the parameters used in the forceField to appropriate
! values. I have to do this because I cannot initialize parameters using real exponents
!=======================================================================================

    subroutine init_ff()

        implicit none

        ljRcr2 = (2.0)**(1.0/3.0) ! Already squared, alphaLJ is 0 for a not soft core potential
        ljRcrA2 = ljRcrA**2 ! Attractive NC-NC interaction cut-off
        ljLGcr2 = ljRcrLG**2
        ljLGcl2 = ljRclLG**2
        dhCR = kBT*exp(-dhRcrLG*dhKappai)/dhRcrLG ! PE shift for DH
        dhCR = kBT*exp(-dhRcrLG*dhKappai)/dhRcrLG ! PE shift for DH
        dhCL = kBT*exp(-dhRclLG*dhKappai)/dhRclLG - dhCR ! PE for DH beyond left cut-off
        tmpUL = tmpU*tmpL ! pre-factor used in trans-membrane potential force calculation

        call build_tanh_table() ! Generate the tanh table used in solvation potential
        call build_table() ! Fill the channel tables

        ! Bip binding/unbinding probabilities
        bip_on = k_bip_on*dt
        bip_off = k_bip_off*dt

    end subroutine init_ff

!=======================================================================================
! Subroutine allocate_tables(), allocate memory for the force tables. Also sets
! some useful parameters for table lookup.
!=======================================================================================

        subroutine allocate_tables()

            implicit none
            integer :: err

            ! First define parameters related to the table lookup
            hNGRIDX = 1 + ceiling((0.0-xCutChanL)/DGRID) ! Location of 0
            NGRIDX = hNGRIDX + ceiling((xCutChanR)/DGRID) ! Total points
            hNGRIDY = 1 + ceiling((0.0-yCutChanL)/DGRID) ! Location of 0
            NGRIDY = hNGRIDY + ceiling((yCutChanR)/DGRID) ! Total points
            hNGRIDZ = 1 + ceiling((0.0-zCutChanL)/DGRID) ! Location of 0
            NGRIDZ = hNGRIDZ + ceiling((zCutChanR)/DGRID) ! Total points
            ! Now allocate tables
            ALLOCATE(table_ljC(NGRIDX+1,NGRIDY+1,NGRIDZ+1),table_ljO(NGRIDX+1,NGRIDY+1,NGRIDZ+1), &
            table_ljCpho(NGRIDX+1,NGRIDY+1,NGRIDZ+1),table_ljOpho(NGRIDX+1,NGRIDY+1,NGRIDZ+1),STAT=err)
            if (err/=0) then
                    write(6,*) 'Error, could not allocate memory for force tables'
                    STOP
            endif
            ALLOCATE(table_chC(NGRIDX+1,NGRIDY+1,NGRIDZ+1),table_chO(NGRIDX+1,NGRIDY+1,NGRIDZ+1),STAT=err)
            if (err/=0) then
                    write(6,*) 'Error, could not allocate memory for charge tables'
                    STOP
            endif
            write(6,*) NGRIDX, hNGRIDX
            write(6,*) NGRIDY, hNGRIDY
            write(6,*) NGRIDZ, hNGRIDZ

        end subroutine allocate_tables

!=======================================================================================
! Subroutine calc_forces(pos,force) is used to calculate the net force on the beads, 
! it calls several other subroutines to calculate the individual force components. 
! If certain force calculations have to change this should be the only place where 
! changes have to be made, simply change which subroutines are being called.
!=======================================================================================

    subroutine calc_forces(pos,force) ! *** Consider moving this to main.f90

        implicit none
        double precision, intent(in) :: pos(:,:) ! Position array
        double precision, intent(inout) :: force(:,:) ! Force array, will be updated
        
        ! Reset the force array
        force = 0.0
        ! I split up the force calculations in a bonded part and a non-bonded part,
        ! this because the bonded forces can be done easily using arrays
        call bonded_forces(pos,force) ! update force array with bonded forces
        !write(6,*) 'Bonded force contribution',force
        !if (nBeads.eq.nFree) write(*,*) 'bond',force(nFree,:)
        call nb_forces(pos,force) ! update force array with non-bonded forces
        !write(6,*) 'Including contribution of pairwise nonbonded',force
        if (LGstatus.eq.0) then ! The LG is closed
            call translocon_forceC(pos,force)
            call position_forcesC(pos,force) ! added the position_forces here due
                                            ! to change in solvation
        else ! The LG is open
            call translocon_forceO(pos,force)
            call position_forcesO(pos,force) ! added the position_forces here due
                                                 ! to change in solvation
        endif
        !write(6,*) 'Including contribution of Ribosome and translocon',force
        ! *** ADDED HERE; TRANSMEMBRANE POTENTIAL FORCE ***
        call transmembrane_potential(pos,force)
        ! BiP binding is evaluated first then the force
        if (BIPON.eq.1) then
            call bip_binding(pos)
            call bip_force(pos,force)
        endif

    end subroutine calc_forces

!=======================================================================================
! Subroutine bip_force(pos,force) updates the force array with forces due to BiP
! binding.
!=======================================================================================

    subroutine bip_force(pos,force)

        implicit none
        double precision, intent(in) :: pos(:,:) ! Position array, used to check distance from translocon exit
        double precision, intent(inout) :: force(:,:) ! Force array, updated here
        integer :: ib ! For iterating over beads

        do ib=1,nForce
            if (bipBinding(ib).eq.1) then ! Only a force if bound to BiP
                if ( ((pos(ib,1)-bipX)**2 + pos(ib,2)**2 + pos(ib,3)**2).lt.4.0 ) then ! Only if near channel exit
                    force(ib,1) = force(ib,1) + F_bip ! Add force in the positive x-direction
                endif
            endif
        enddo

    end subroutine bip_force

!=======================================================================================
! Subroutine bip_binding(pos) updates which beads are bound to BiP based on the
! current positions and the bip binding rates
!=======================================================================================

    subroutine bip_binding(pos)

        implicit none
        double precision, intent(in) :: pos(:,:) ! Positions are not updated
        double precision :: ranNum ! Store a random number
        integer :: ib ! Iterate over beads

        ! First check whether beads are past the membrane
        do ib=1,nForce
            if (pos(ib,1).gt.2.0) then
                transSide(ib) = 1.0
            else
                transSide(ib) = 0.0
            endif
        enddo

        ! Next determine bip binding/unbinding
        do ib=1,nForce,4
            if (sum(transSide(ib:ib+3)).eq.4) then ! All beads in this segment are past the membrane
                ! We will assess whether the binding status changes based on a
                ! random number
                call Random_Number(ranNum)
                random_count = random_count + 1 ! Update the number of random numbers generated, sys_param.f90
                if (sum(bipBinding(ib:ib+3)).eq.0) then ! Currently not bound
                    ! Assess if Bip binds
                    if (ranNum<bip_on) then ! bip_on is initialized as k_on*dt
                        bipBinding(ib:ib+3) = 1.0 ! Now bound
                    endif
                elseif (sum(bipBinding(ib:ib+3)).eq.4) then ! Currently bound
                    ! Assess if Bip unbinds
                    if (ranNum<bip_off) then ! bip_off is initialized as k_off*dt
                        bipBinding(ib:ib+3) = 0.0 ! Now no longer bound
                    endif
                endif
            endif
        enddo

    end subroutine bip_binding

!=======================================================================================
! Subroutine transmembrane_potential(pos,force) is used to calculate the force
! on charged beads due to the transmembrane potential. 
!=======================================================================================

    subroutine transmembrane_potential(pos,force)

        implicit none
        double precision, intent(in) :: pos(:,:) ! Position array
        double precision, intent(inout) :: force(:,:) ! Force array
        integer :: ib ! bead index
        double precision :: xb, qb, etmp, ftmp ! x pos, bead charge, exponent, force

        ! Loop only over charged beads
        do ib = 1,nCharge
            xb = pos(chBeads(ib),1) ! Position of the bead along the membrane axis
            qb = bCharge(chBeads(ib)) ! Bead charge
            etmp = exp(tmpL*xb) ! Pre-calculate exponent
            ftmp = qb*tmpUL*etmp/((1+etmp)**2)
            force(chBeads(ib),1) = force(chBeads(ib),1) + ftmp ! Update force
        enddo

    end subroutine transmembrane_potential

!=======================================================================================
! Subroutine bonded_forces(pos,force) updates the force array with bonded forces that
! depend on the current position. If we want to change the kind of interactions for 
! bonded particles then the calls in this subroutine need to be updated.
!=======================================================================================

    subroutine bonded_forces(pos,force)

        implicit none
        double precision, intent(in) :: pos(:,:) ! Position array
        double precision, intent(inout) :: force(:,:) ! force array
        double precision :: dr(nForce-1,NDIM), r2(nForce-1) ! bond vector and magnitude^2 

        ! Create an array with bondvectors
        dr = pos(2:nForce,:) - pos(1:nForce-1,:) ! bondvectors
        r2 = sum(dr**2,dim=2) ! Bond vector magnitude squared
                                      
        ! Now update the force array with the bonded forces
        call bond_fene(dr,r2,force) ! Need the bond vectors to calculate the force
        call bond_lj_repulsive(dr,r2,force) ! Use bond vectors to calculate forces

    end subroutine bonded_forces

!=======================================================================================
! Subroutine nb_forces(pos,force) updates the force array with the non-bonded forces
! that depend on the current position. To change the interactions update the calls here.
!=======================================================================================

    subroutine nb_forces(pos,force)

        implicit none
        double precision, intent(in) :: pos(:,:)
        double precision, intent(inout) :: force(:,:)   ! by using inout Im overwriting
                                ! the force array
        integer :: ib,jb    ! Used for do loops over the number of particles
        double precision :: dr(NDIM), ftmp(NDIM) ! Stores bond vector, force vector
        double precision :: q1, q2 ! Relate to charges

        ! Now loop over all particles and calculate the relevant forces.
        ! Note that forces between neighboring atoms are already calculated
        ! so we only need to loop over non-neighbors
        do ib = 1,nForce-2  ! nForce, only beads that are in the system 
            do jb = ib+2,nForce
                dr(:) = pos(jb,:) - pos(ib,:)
                ! *** IF statement added to allow for weak attraction between TMDs
                if (pairwise(ib,jb).eq.0.d0) then
                    call lj_repulsive(dr,ljEps,ljRcr2,ftmp) ! Parameters eps and cut-off necessary
                else
                    call lj_repulsive(dr,pairwise(ib,jb),ljRcrA2,ftmp) ! Parameters for attractive LJ (cut-off beyond 2^(1/6) sigma)
                endif
                force(ib,:) = force(ib,:) - ftmp(:)
                force(jb,:) = force(jb,:) + ftmp(:)
            enddo
        enddo

        ! Electrostatics between NC beads
        do ib = 1,nCharge-1 ! Only charged beads are included
            do jb = ib+1,nCharge
                ! Only interact if both the beads are translated
                if (chBeads(jb).gt.nForce) cycle ! Skip this pair
                dr(:) = pos(chBeads(jb),:) - pos(chBeads(ib),:)
                q1 = bCharge(chBeads(ib)); q2 = bCharge(chBeads(jb));
                call dh_elec(dr,dhcl2,q1,q2,ftmp)
                force(chBeads(ib),:) = force(chBeads(ib),:) - ftmp(:)
                force(chBeads(jb),:) = force(chBeads(jb),:) + ftmp(:)
            enddo
        enddo

    end subroutine nb_forces

!=======================================================================================
! Subroutine position_forcesO(pos,force) updates the force array with forces that depend
! only on the position of each bead. This includes forces due to solvation, and
! interaction with the ribosome/translocon.
!=======================================================================================

    subroutine position_forcesO(pos,force)

        implicit none
        double precision, intent(in) :: pos(:,:) ! positions are not changed
        double precision, intent(inout) :: force(:,:) ! forces are updated
        integer :: ib ! used for loops
        double precision :: r(NDIM),f(NDIM) ! position vector
        double precision :: rsym(NDIM-1), fsym(NDIM-1) ! Symmetric basis
        
        do ib = 1,nForce ! Position forces only for free beads!***!
            r = pos(ib,:)
            f(:) = 0.0 ! Initialize
! **** USE THIS RIBOSOME IF I WANT TO PUT SOME SPHERICAL EXCLUDED VOLUME SOMEWHERE, 
! prevents NC from moving back into the ribosomal exit tunnel which isnt
! explicitely included in the simulations
            call ribosome_force(r,f) ! Ribosome
            force(ib,:) = force(ib,:) + f ! Add to the force array
            ! The next few contributions are radially symmetric so I switch basis here
            fsym(:) = 0.0 ! Initialize
            rsym(1) = abs(r(1)) ! X-coordinate absolute symmetry 
            rsym(2) = sqrt(r(2)**2 + r(3)**2) ! Spherically symmetric wrt y,z
            call solvation_forceO(rsym,dGtransfer(ib),fsym) ! solvation force
            ! Force contributions have to be projected back on the original basis
            if (rsym(1).gt.machine_eps) then
                force(ib,1) = force(ib,1) + fsym(1)*r(1)/rsym(1)
            endif
            if (rsym(2).gt.machine_eps) then
                force(ib,2) = force(ib,2) + fsym(2)*r(2)/rsym(2)
                force(ib,3) = force(ib,3) + fsym(2)*r(3)/rsym(2)
            endif
        ! *** ADDED HERE; RESTRAINT TO KEEP BEADS IN INTRA OR EXTRA CELLULAR
            if (pos(ib,1).gt.-2.d0) force(ib,1) = force(ib,1) - intracellular(ib)
            if (pos(ib,1).lt.2.d0) force(ib,1) = force(ib,1) + extracellular(ib)
        enddo

    end subroutine position_forcesO
!=======================================================================================
! Subroutine position_forcesC(pos,force) updates the force array with forces that depend
! only on the position of each bead. This includes forces due to solvation, and
! interaction with the ribosome/translocon.
!=======================================================================================

    subroutine position_forcesC(pos,force)

        implicit none
        double precision, intent(in) :: pos(:,:) ! positions are not changed
        double precision, intent(inout) :: force(:,:) ! forces are updated
        integer :: ib ! used for loops
        double precision :: r(NDIM),f(NDIM) ! position vector
        double precision :: rsym(NDIM-1), fsym(NDIM-1) ! Symmetric basis
        
        do ib = 1,nForce ! Position forces only for free beads!***!
            r = pos(ib,:)
            f(:) = 0.0 ! Initialize
! **** USE THIS RIBOSOME IF I WANT TO PUT SOME SPHERICAL EXCLUDED VOLUME SOMEWHERE
            call ribosome_force(r,f) ! Ribosome
            force(ib,:) = force(ib,:) + f ! Add to the force array
            ! The next few contributions are radially symmetric so I switch basis here
! **** I HAVE TO REVISIT IF THIS ACTUALLY IS MORE EFFICIENT, IT MAY NOT BE ****
            fsym(:) = 0.0 ! Initialize
            rsym(1) = abs(r(1)) ! X-coordinate absolute symmetry 
            rsym(2) = sqrt(r(2)**2 + r(3)**2) ! Spherically symmetric wrt y,z
            call solvation_forceC(rsym,dGtransfer(ib),fsym) ! solvation force
            ! Force contributions have to be projected back on the original basis
            if (rsym(1).gt.machine_eps) then
                force(ib,1) = force(ib,1) + fsym(1)*r(1)/rsym(1)
            endif
            if (rsym(2).gt.machine_eps) then
                force(ib,2) = force(ib,2) + fsym(2)*r(2)/rsym(2)
                force(ib,3) = force(ib,3) + fsym(2)*r(3)/rsym(2)
            endif
        ! *** ADDED HERE; RESTRAINT TO KEEP BEADS IN INTRA OR EXTRA CELLULAR
            if (pos(ib,1).gt.-2.d0) force(ib,1) = force(ib,1) - intracellular(ib)
            if (pos(ib,1).lt.2.d0) force(ib,1) = force(ib,1) + extracellular(ib)
        enddo

    end subroutine position_forcesC

!=======================================================================================
! Subroutine on_grid(r,ir) assess wether a CG bead is within range of the translocon.
!=======================================================================================

    subroutine on_grid(r,ir)
    
        implicit none
        double precision, intent(in) :: r(:) ! position
        integer, intent(out) :: ir
        ir = 0
        if ((r(1).gt.xCutChanL).and.(r(1).lt.xCutChanR)) then
            if ((r(2).gt.yCutChanL).and.(r(2).lt.yCutChanR)) then
                if ((r(3).gt.zCutChanL).and.(r(3).lt.zCutChanR)) then
                    ir = 1 ! Inside grid
                endif
            endif
        endif

    end subroutine on_grid

!=======================================================================================
! Subroutine table_force(r,f,table) returns position dependent forces on a
! particle from a table. LINEAR INTERPOLATION OF GRADIENT ALONG TABLE EDGES.
!=======================================================================================

! *** THIS IS AN ALTERNATE VERSION THAT DIRECTLY CALCULATES THE PE GRADIENT,
! INSTEAD OF INTERPOLATING FORCES *** !

        subroutine table_force(r,f,table)

            implicit none
            double precision, intent(in) :: r(:), table(:,:,:)
            double precision, intent(inout) :: f(:) ! Contain force + PE
            double precision :: xbr, ybr, zbr, phix, phiy, phiz ! interpolation
            double precision :: phixi, phiyi, phizi ! 1-lambda
            integer :: xb, yb, zb, ic ! Bin indices

            ! First check if the particle is on the grid
            call on_grid(r,ic) ! 0 if the position is not on the grid
            if (ic.eq.1) then ! This bead feels the channel
                xbr = hNGRIDX + DGRIDi*r(1) ! x in bin coordinates
                ybr = hNGRIDY + DGRIDi*r(2)
                zbr = hNGRIDZ + DGRIDi*r(3)
                xb = floor(xbr); yb = floor(ybr); zb = floor(zbr)  ! Discrete indices
                phix = xbr-xb; phiy = ybr-yb; phiz = zbr-zb;! Measure of how close to each of the table elements we are
                phixi = 1.d0-phix; phiyi = 1.d0-phiy; phizi = 1.d0-phiz;
                ! Determine forces one direction at a time
                f(1) = f(1) - ((table(xb+1,yb,zb)-table(xb,yb,zb))*phiyi*phizi &
                    + (table(xb+1,yb+1,zb)-table(xb,yb+1,zb))*phiy*phizi &
                    + (table(xb+1,yb,zb+1)-table(xb,yb,zb+1))*phiyi*phiz &
                    + (table(xb+1,yb+1,zb+1)-table(xb,yb+1,zb+1))*phiy*phiz)*DGRIDi
                f(2) = f(2) - ((table(xb,yb+1,zb)-table(xb,yb,zb))*phixi*phizi &
                    + (table(xb+1,yb+1,zb)-table(xb+1,yb,zb))*phix*phizi &
                    + (table(xb,yb+1,zb+1)-table(xb,yb,zb+1))*phixi*phiz &
                    + (table(xb+1,yb+1,zb+1)-table(xb+1,yb,zb+1))*phix*phiz)*DGRIDi
                f(3) = f(3) - ((table(xb,yb,zb+1)-table(xb,yb,zb))*phiyi*phixi &
                    + (table(xb,yb+1,zb+1)-table(xb,yb+1,zb))*phiy*phixi &
                    + (table(xb+1,yb,zb+1)-table(xb+1,yb,zb))*phiyi*phix &
                    + (table(xb+1,yb+1,zb+1)-table(xb+1,yb+1,zb))*phiy*phix)*DGRIDi
            endif

    end subroutine table_force

!=======================================================================================
! Subroutine table_energy(r,f,table) returns position dependent potential energy of a
! particle from a table. LINEAR INTERPOLATION OF VALUES IN TABLE.
!=======================================================================================

       subroutine table_energy(r,u,table)

           implicit none
           double precision, intent(in) :: r(:), table(:,:,:)
           double precision, intent(inout) :: u ! PE
           double precision :: xbr, ybr, zbr, phix, phiy, phiz ! interpolation
           double precision :: phixi, phiyi, phizi ! 1-lambda
           integer :: xb, yb, zb, ic ! Bin indices

           ! First check if the particle is on the grid
           call on_grid(r,ic) ! 0 if the position is not on the grid
           if (ic.eq.1) then ! This bead feels the channel
               xbr = hNGRIDX + DGRIDi*r(1) ! x in bin coordinates
               ybr = hNGRIDY + DGRIDi*r(2)
               zbr = hNGRIDZ + DGRIDi*r(3)
               xb = floor(xbr); yb = floor(ybr); zb = floor(zbr)  ! Discrete indices
               phix = xbr-xb; phiy = ybr-yb; phiz = zbr-zb;! Measure of how close to each of the table elements we are
               phixi = 1.d0-phix; phiyi = 1.d0-phiy; phizi = 1.d0-phiz;
               u = u + table(xb,yb,zb)*phixi*phiyi*phizi & ! For PE its just a linear interpolation
                   + table(xb+1,yb,zb)*phix*phiyi*phizi   &
                   + table(xb,yb+1,zb)*phixi*phiy*phizi  &
                   + table(xb+1,yb+1,zb)*phix*phiy*phizi   &
                   + table(xb,yb,zb+1)*phixi*phiyi*phiz & 
                   + table(xb+1,yb,zb+1)*phix*phiyi*phiz   &
                   + table(xb,yb+1,zb+1)*phixi*phiy *phiz  &
                   + table(xb+1,yb+1,zb+1)*phix*phiy*phiz
           endif

   end subroutine table_energy

!=======================================================================================
! Subroutine translocon_forceO(pos,force) calculates forces on a CG particle due to
! interactions with the translocon, for the OPEN translocon.
!=======================================================================================

    subroutine translocon_forceO(pos,force)
        
            implicit none
            double precision, intent(in) :: pos(:,:)
            double precision, intent(inout) :: force(:,:)
            double precision :: r(NDIM), f(NDIM)
            integer :: ib, cb

            do ib=1,nForce
                r(:) = pos(ib,:)
                f(:) = 0.0 ! re-set
                call table_force(r,f,table_ljO(:,:,:)) ! LJ forces
                force(ib,:) = force(ib,:) + f*(1.d0-oLambda(ib)) ! Update
                f(:) = 0.0 ! re-set
                call table_force(r,f,table_ljOpho(:,:,:))
                force(ib,:) = force(ib,:) + f*oLambda(ib) ! Update
            enddo
            do ib=1,nCharge ! Also the contribution due to electrostatics
                cb = chBeads(ib)
                if (cb.le.nForce) then ! Only if it is feeling force
                    r(:) = pos(cb,:); f(:) = 0.0; ! reset
                    call table_force(r,f,table_chO(:,:,:))
                    force(cb,:) = force(cb,:) + bCharge(cb)*f(:)
                endif
            enddo

        end subroutine translocon_forceO

!=======================================================================================
! Subroutine translocon_forceC(pos,force) calculates forces on a CG particle due to
! interactions with the translocon, for the CLOSED translocon.
!=======================================================================================

        subroutine translocon_forceC(pos,force)
        
            implicit none
            double precision, intent(in) :: pos(:,:)
            double precision, intent(inout) :: force(:,:)
            double precision :: r(NDIM), f(NDIM)
            integer :: ib, cb

            do ib=1,nForce
                r(:) = pos(ib,:)
                f(:) = 0.0 ! re-set
                call table_force(r,f,table_ljC(:,:,:)) ! LJ forces
                force(ib,:) = force(ib,:) + f*(1.d0-cLambda(ib)) ! Update
                f(:) = 0.0 ! re-set
                call table_force(r,f,table_ljCpho(:,:,:))
                force(ib,:) = force(ib,:) + f*cLambda(ib) ! Update
            enddo
            do ib=1,nCharge ! Also the contribution due to electrostatics
                cb = chBeads(ib)
                if (cb.le.nForce) then ! Only if it is feeling force
                    r(:) = pos(cb,:); f(:) = 0.0; ! reset
                    call table_force(r,f,table_chC(:,:,:))
                    force(cb,:) = force(cb,:) + bCharge(cb)*f(:)
                endif
            enddo

        end subroutine translocon_forceC

!=======================================================================================
! Subroutine solvation_forceO(pos,phobicity,force) calculates the position dependent force on a 
! particle due to the solvation environment (lipid/water).
!=======================================================================================

    subroutine solvation_forceO(dr,g,f)

        implicit none
        double precision, intent(in) :: dr(:), g ! Position and hydrophobicity
        double precision, intent(inout) :: f(:) ! Force will be updated
        double precision ::  ftmp
        double precision :: Sx, dSdx, Sx2, dS2dx, Sy, dSdy

        ! Check if we are near the membrane
        if (dr(1)<RIGHTCUT) then
            ! switch function calculates Sx, dSdx for the switch function
            ! as described in the Cell Reports (2012) paper
            call switch_function(dr(1), MBL, MBR, Sx, dSdx, 0.25d0)
            call switch_function(dr(2), -rChannel, rChannel, Sy, dSdy, 0.25d0)
            Sy = 1-Sy; dSdy = -dSdy ! More convenient
            f(1) = f(1) - g * dSdx * Sy 
            f(2) = f(2) - g * Sx * dSdy
            ! this is the filling channel part
            f(1) = f(1) - DGFILLO * dSdx  + dSdx * Sy * DGFILLO  
            f(2) = f(2) +  Sx * dSdy * DGFILLO 
!           endif
        endif

    end subroutine solvation_forceO
!=======================================================================================
! Subroutine solvation_forceC(pos,phobicity,force) calculates the position dependent force on a 
! particle due to the solvation environment (lipid/water).
!=======================================================================================

    subroutine solvation_forceC(dr,g,f)

        implicit none
        double precision, intent(in) :: dr(:), g ! Position and hydrophobicity
        double precision, intent(inout) :: f(:) ! Force will be updated
        double precision ::  ftmp
        double precision :: Sx, dSdx, Sx2, dS2dx, Sy, dSdy, Sy2, dS2dy

        ! Check if we are near the membrane
        if (dr(1)<RIGHTCUT) then
            ! switch function calculates Sx, dSdx for the switch function
            ! as described in the Cell Reports (2012) paper
            call switch_function(dr(1), MBL, MBR, Sx, dSdx, 0.25d0)
            call switch_function(dr(2), -rChannel, rChannel, Sy, dSdy, 0.25d0)
            Sy = 1-Sy; dSdy = -dSdy ! More convenient
            f(1) = f(1) -g * dSdx * Sy
            f(2) = f(2) -g * Sx * dSdy 
            ! this is the filling channel part
            f(1) = f(1) - dSdx * DGFILLC + dSdx * Sy * DGFILLC
            f(2) = f(2) +  Sx * dSdy * DGFILLC
!           endif
        endif

    end subroutine solvation_forceC

!=======================================================================================
! Subroutine ribosome_force(pos,force) calculates the forces due to interaction with
! the ribosome. These are modeled to be purely repulsive.
!=======================================================================================

! *** At the moment the ribosome is a sphere with a center and radius defined in sys_param

    subroutine ribosome_force(r,f)

        implicit none
        double precision, intent(in) :: r(:) ! Dont change here
        double precision, intent(inout) :: f(:) ! Update forces
        double precision :: dr(NDIM), rmag, r2, ri2, ri6
        ! *** This subroutine is very simple now, may want to add some
        ! *** parameters instead
        dr = r-riboCenter
        rmag = sqrt(sum(dr**2))
        dr = (rmag-riboR)*dr/rmag
        r2 = sum(dr**2)
        f(:) = 0.0 ! Zero if outside of cutoff
        if (r2.lt.2.d0**(1.d0/3.d0)) then ! Currently the force is purely repulsive
            ri2 = 1.0/r2
            ri6 = ri2*ri2*ri2
            f = 48.0*dr*ri2*ri6*(ri6-0.5)
        endif

    end subroutine ribosome_force

!=======================================================================================
! Subroutine bond_fene(dr,r2,force) calculates the attractive component of the bonded
! forces. The repulsive component is modelled as the repulsive part of a LJ potential
! and it is calculated by subroutine bond_lj_repulsive(dr,r2,force).
!=======================================================================================

    subroutine bond_fene(dr,r2,force)

        implicit none
        double precision, intent(in) :: dr(:,:), r2(:)
        double precision, intent(inout) :: force(:,:)
        double precision :: f(nForce-1) ! Local force vector ***
        integer :: nd ! Local integer for loop

        f = -feneK/(1.0-r2/feneR2)
        do nd = 1,NDIM
            force(1:nForce-1,nd) = force(1:nForce-1,nd) - f*dr(:,nd) ! *** nForce-1 isnt good
            force(2:nForce,nd) = force(2:nForce,nd) + f*dr(:,nd) ! *** Do something instead of
                                    ! *** nBeads if not all beads move
        enddo

    end subroutine bond_fene

!=======================================================================================
! Subroutine bond_lj_repulsive(dr,r2,force) calculates the repulsive component of the
! bonded forces. It is the repulsive part of a 12-6 LJ potential.
!=======================================================================================

    subroutine bond_lj_repulsive(dr,r2,force)

        implicit none
        double precision, intent(in) :: dr(:,:), r2(:)
        double precision, intent(inout) :: force(:,:)
        double precision :: f(NDIM), ri2(nForce-1), ri6(nForce-1) ! Local vectors ***
        integer :: ib ! Local integer for loops

        ri2 = 1.0/r2
        ri6 = ri2*ri2*ri2
        do ib = 1,nForce-1
            if (r2(ib).lt.ljRcr2) then
                f = 48.0*dr(ib,:)*ri2(ib)*ri6(ib)*ljEps*(ri6(ib)-0.5)
                force(ib,:) = force(ib,:) - f! Update force
                force(ib+1,:) = force(ib+1,:) + f
            endif
        enddo

    end subroutine bond_lj_repulsive

!=======================================================================================
! Subroutine lj_repulsive(dr,eps,rcut2,f) calculates the short-ranged non-bonded 
! interactions between pairs of particles.
!=======================================================================================

    subroutine lj_repulsive(dr,eps,rcut2,f)

        implicit none
        double precision, intent(in) :: dr(:), eps, rcut2
        double precision, intent(out) :: f(NDIM)
        double precision :: r2,ri2,ri6 ! r^2, 1/r^2, 1/r^6

        ! Calculate the repulsive part of the bond potential
        r2 = sum(dr**2)
        f = 0.0 ! First zerorize
        if (r2 .lt. rcut2) then ! Now assign a nonzero value if inside cutoff
            ri2 = 1.0/r2
            ri6 = ri2*ri2*ri2
            f = 48.0*ri2*ri6*eps*(ri6-0.5)*dr ! Force
        endif

    end subroutine lj_repulsive

!=======================================================================================
! Subroutine dh_elec(dr,lcut2,q1,q2,f) calculates the electrostatic interactions
! between a pair of CG particles using the Debye-Huckel potential. Returns the
! force only and no right cut-off is applied, meant for NC-NC interactions.
!=======================================================================================

        subroutine dh_elec(dr,lcut2,q1,q2,f)

                implicit none
                double precision, intent(in) :: dr(:), lcut2, q1, q2
                double precision, intent(out) :: f(:)
                double precision :: r2, r

                r2 = sum(dr**2)
                f(:) = 0.0; ! Initialize
                if (r2.gt.lcut2) then ! Beyond the left cut-off, no force otherwise
                        r = sqrt(r2)
                        f(:) = q1*q2*kBT*exp(-dhKappai*r)*(1.0/r + dhKappai)*dr(:)/r2
                endif
                       

        end subroutine dh_elec

!=======================================================================================
! Subroutine dh_pe(dr,lcut2,rcut2,q1,q2,f,u) calculates the electrostatic interactions
! between a pair of CG particles using the Debye-Huckel potential. Returns the
! force as well as the PE (needed for LG dynamics).
!=======================================================================================

    subroutine dh_pe(dr,lcut2,rcut2,q1,q2,f,u)

        implicit none
        double precision, intent(in) :: dr(:), lcut2, rcut2, q1, q2
        double precision, intent(out) :: f(:), u
        double precision :: r2, r

        r2 = sum(dr**2)
        f(:) = 0.0; u = 0.0; ! Initialize
        if (r2.lt.rcut2) then ! We are within the right cut-off
            if (r2.gt.lcut2) then ! Beyond the left cut-off
                r = sqrt(r2)
                u = kBT*exp(-dhKappai*r)/r
                f = q1*q2*u*(1.0/r2 + dhKappai/r)*dr
                u = q1*q2*(u-dhCR) ! Now shift the potential to make it go to zero
            else
                u = q1*q2*dhCL ! Potential is left cut-off max
                ! Force is zero here
            endif
        endif
               

    end subroutine dh_pe

!=======================================================================================
! Subroutine lj_soft_pe(dr,eps,rcut2,f,u) calculates the short-ranged non-bonded 
! interactions between pairs of particles. Returns the force as well as the PE.
! SOFT-CORE LENNARD JONES
!=======================================================================================

    subroutine lj_soft_pe(dr,eps,rcut2,f,u,sigma2)

        implicit none
        double precision, intent(in) :: dr(:), eps, rcut2, sigma2
        double precision, intent(inout) :: f(:), u
        double precision :: r2,ari6 ! r^2, 1/(a+r^6)
        double precision :: alpha, ecut ! Locally determined now to allow changes in epsilon
        double precision :: cut2

        ! Determine alpha and ecut, *** INEFFICIENT BUT FLEXIBLE ***
        alpha = 2.d0*eps*(sqrt(1+maxE/eps)-1)/maxE
        ! Check if this is a purely repulsive bead first
        if (rcut2.lt.machine_eps) then
            cut2 = (2.d0-alpha)**(1.0/3.0)
        !    sigma2=1.44 ! the repulsive beads have a larger radius
        else
            cut2=rcut2
        endif
        r2= sum(dr**2) / sigma2
        ! Next calculate the shift at the cutoff
        ! Calculate the repulsive part of the bond potential
        u = 0.d0; f(:) = 0.d0 ! First zerorize
        if (r2 .lt. cut2) then ! Now assign a nonzero value if inside cutoff
            ari6 = 1/(alpha+r2*r2*r2)
            f = 48.d0*r2*r2*eps*ari6*ari6*(ari6-0.5d0)*dr/sigma2 ! Force
            ecut = (1/(cut2**3+alpha) - 1)/(cut2**3+alpha)
            u = 4.d0*eps*(ari6*(ari6-1.d0) - ecut) ! PE
        endif

    end subroutine lj_soft_pe

!=======================================================================================
! Subroutine build_table(cType,cBeads,epsB,charges) Populates PE and Force 
! tables for a channel. Separate tables are made for the LJ and electrostatic
! forces.
!=======================================================================================

    subroutine build_table()

        implicit none
        integer :: ib, xb, yb, zb ! Used to iterate
        double precision :: xv, yv, zv, rv! Keep track of grid point coords
        double precision :: lambdaM ! Switching parameter between water and lipid and derivatives
        double precision :: Sx, Sy, dSx, dSy ! used to calculate switching parameters
        double precision :: dr(NDIM), fp(NDIM), up, fh(NDIM), uh , fr(NDIM), ur ! Used in force calculation
        double precision :: fw(NDIM), fm(NDIM), uw, um ! Water/Membrane part

        fh(:) = 0.d0; fp(:) = 0.d0; dr(:) = 0.d0; ! Zerorize
        ! We go through every point in the force/energy table
        do zb=1,NGRIDZ
write(6,*) 'Calculating channel force table, progress at ',(zb*100.d0/NGRIDZ),'%'
            do yb=1,NGRIDY
                do xb=1,NGRIDX
                    ! Initialize for this gridpoint
                    xv = (xb-hNGRIDX)*DGRID ! coords of gridpoint
                    yv = (yb-hNGRIDY)*DGRID
                    zv = (zb-hNGRIDZ)*DGRID
                    rv = sqrt(yv**2+zv**2) ! Radial coordinate
                    ! Determine the type of solvent we are in (water/lipid)
                    call switch_function(xv,MBL,MBR,Sx,dSx, 0.25d0)
                    call switch_function(rv,-rChannel,rChannel,Sy,dSy, 0.25d0)
                    lambdaM = Sx*(1.d0-Sy) ! lambdaM = 1 (lipid), 0 (water)
                    ! For every bead in the channel/ribosome add its
                    ! contribution to the force/energy table
                    do ib=1,nBeadChan
                        ! Vector between the gridpoint and bead
                        dr(1) = xv-chanBeads(ib,1)
                        dr(2) = yv-chanBeads(ib,2)
                        dr(3) = zv-chanBeads(ib,3)
                        ! Get force and PE and update the correct table
                        ! Interaction with a phobic particle
                        call lj_soft_pe(dr,chanP(ib,3),chanP(ib,7),fw,uw,chanP(ib,11))
                        call lj_soft_pe(dr,chanP(ib,5),chanP(ib,9),fm,um,1.d0)
                        up = lambdaM*um + (1.d0-lambdaM)*uw
                        call lj_soft_pe(dr,chanP(ib,4),chanP(ib,8),fw,uw,chanP(ib,12))
                        call lj_soft_pe(dr,chanP(ib,6),chanP(ib,10),fm,um,1.d0)
                        uh = lambdaM*um + (1.d0-lambdaM)*uw
                        if (chanP(ib,1).gt.1e-5) then ! Closed channel bead
                            table_ljC(xb,yb,zb) = table_ljC(xb,yb,zb) + uh
                            table_ljCpho(xb,yb,zb) = table_ljCpho(xb,yb,zb) + up
                        endif
                        if (chanP(ib,2).gt.1e-5) then ! Open channel bead
                            table_ljO(xb,yb,zb) = table_ljO(xb,yb,zb) + uh
                            table_ljOpho(xb,yb,zb) = table_ljOpho(xb,yb,zb) + up
                        endif                                      
                        ! Next charge interactions
                        if (chanCharge(ib)/=0) then
                            call dh_pe(dr,dhRclLG2,dhRcrLG2,chanCharge(ib),1.d0,fp,up)
                            if (chanP(ib,1).gt.1e-5) then ! Closed channel bead
                                table_chC(xb,yb,zb)=table_chC(xb,yb,zb)+up
                            endif
                            if (chanP(ib,2).gt.1e-5) then ! open channel bead
                                table_chO(xb,yb,zb)=table_chO(xb,yb,zb)+up
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo
        write(6,*) 'Channel tables generated'
        write(6,*) 'de-allocating channel property data'
        DEALLOCATE(chanBeads,chanCharge,chanP,chanType)
                                
        end subroutine build_table

!=======================================================================================
! The following code block is taken directly from Bins code *** MODIFY IT ***
! It is used for generating tables for the tanh function as well as for the LJ
! interactions with the closed channel (Some parts are modified wrt the channel
! interactions, modifications are marked with ***[description of change]***.
!=======================================================================================

! Function S(x,psi,phi) from the Cell Reports (2012) paper
    subroutine switch_function(x, LB, RB, Sx, dSdx, CUTB)
        implicit none
        double precision, intent(in) :: x, LB, RB, CUTB
        double precision, intent(out) :: Sx, dSdx

        double precision :: tanh_mbl, tanh_mbr

        tanh_mbl = tanh_from_table((x-LB)/CUTB)
        tanh_mbr = tanh_from_table((x-RB)/CUTB)
        Sx = 0.25*(1+tanh_mbl) * (1-tanh_mbr)
        dSdx = 0.25 * (1-tanh_mbl**2) / CUTB * (1-tanh_mbr) + &
               0.25 * (1+tanh_mbl) * (-1) * (1-tanh_mbr**2) / CUTB

    end subroutine

! Simpler version of the same function, doesnt calculate dSdx

! Generate a table for the tanh function, it will be used extensively in potentials
! containing a switch function
    subroutine build_tanh_table()

        implicit none
        integer :: ig
        double precision :: r

        table_stepsize =  (RIGHTCUT-LEFTCUT)/ (NGRID-1)

        r = LEFTCUT
        do ig=1, NGRID
            r = r + table_stepsize
            tabulated_tanh(ig) = tanh(r)
        enddo

    end subroutine

! Extract a tanh value from the generated table
    function tanh_from_table(r )

        implicit none
        double precision, intent(in) :: r
        double precision :: tanh_from_table

        integer :: tablepos, iindex
        double precision :: dindex, phi

        if (r < LEFTCUT) then
            tanh_from_table = -1.0
            return
        endif
        if (r > RIGHTCUT) then
            tanh_from_table = 1.0
            return
        endif

        dindex = (r-LEFTCUT) / table_stepsize
        iindex = floor(dindex)
        phi = dindex - iindex
        tablepos = iindex + 1 ! b/c fortran index starts from 1
        tanh_from_table = tabulated_tanh(tablepos)*(1-phi) + &
                          tabulated_tanh(tablepos+1)*phi

        return
    end function



end module forceField
