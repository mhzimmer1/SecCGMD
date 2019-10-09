!=========================================================================================
! Module containing subroutines related to channel dynamics. This is completely independent
! of the brownian dynamics by which the polymer chain is evolved.
!=========================================================================================

module channel

    use sys_param
    use forceField
    use polymer, only : dGtransfer,cLambda,oLambda ! TFE for each bead 

    implicit none
    double precision :: LGfricti=dtRealTime/(dt*tauLG) ! 1/tau_{LG} in cell reports 2012

    contains

!=========================================================================================
! Subroutine LG_dynamics(pos) time evolution of the lateral gate based on current position
!=========================================================================================

    subroutine LG_dynamics(r)

        implicit none
        double precision, intent(in) :: r(:,:) ! polymer chain positions
        double precision :: dG, ranNum, maxP ! The change in free energy for opening the LG, a random number, and the maximum possible probability of changing gate state

        call Random_Number(ranNum) ! Need a random number for this stochastic opening/closing
        maxP = dt*LGfricti
        if (ranNum < maxP) then
            dG = DGEMPTY ! Initialize ! Empty Channel value
            call LGenergy(r,dG,nBeads) ! Add a contribution due to polymer - translocon interaction energy
            call LGsolv(r,dG,nBeads) !solvation change due to channel conf. change
            call LGopen(dG, ranNum) ! Evolve the LG status
        endif        

    end subroutine LG_dynamics  

!=========================================================================================
! LGenergy(pos,dU) updates the change in energy associated with opening the LG based on 
! interaction between beads and the lateral gate.
!=========================================================================================

    subroutine LGenergy(rvec,dGopen,nb)

        implicit none
        double precision, intent(in) :: rvec(:,:) ! Positions
        double precision, intent(inout) :: dGopen ! dG for LG opening
        integer, intent(in) :: nb ! Number of positions
        double precision :: r(NDIM), Uo, Uc, utmp ! location and Uopen, Uclosed
        integer :: ib, cb ! Used for loop

        ! For each bead we have an open and a closed contribution to dGopen,
        ! both are from table lookup
        Uo = 0.0; Uc = 0.0; ! Initialize
        do ib = 1,nb
            r = rvec(ib,1:NDIM)
            utmp = 0.d0
            call table_energy(r,utmp,table_ljO(:,:,:)) ! LJ for open
            Uo = Uo + (1.d0-oLambda(ib))*utmp; utmp = 0.d0
            call table_energy(r,utmp,table_ljC(:,:,:)) ! LJ for closed
            Uc = Uc + (1.d0-cLambda(ib))*utmp; utmp = 0.d0
            call table_energy(r,utmp,table_ljOpho(:,:,:)) ! LJ for open
            Uo = Uo + oLambda(ib)*utmp; utmp = 0.d0
            call table_energy(r,utmp,table_ljCpho(:,:,:)) ! LJ for closed
            Uc = Uc + cLambda(ib)*utmp; utmp = 0.d0
        enddo
        do ib=1,nCharge ! Also the contribution due to electrostatics
            cb = chBeads(ib)
            if (cb.le.nForce) then ! Only if it is feeling force
                r(:) = rvec(cb,1:NDIM)
                call table_energy(r,utmp,table_chO(:,:,:)) !open
                Uo=Uo+bCharge(cb)*utmp; utmp = 0.d0
                call table_energy(r,utmp,table_chC(:,:,:)) !closed
                Uc=Uc+bCharge(cb)*utmp; utmp = 0.d0
            endif
        enddo
        dGopen = dGopen + (Uo-Uc) ! Update dGopen

    end subroutine LGenergy


!=========================================================================================
! LGsolv(pos,dU, nb) updates the change in energy associated with opening the LG based on 
! the change in solvent between open and closed channel states
!=========================================================================================

    subroutine LGsolv(rvec, dGopen, nb)
        implicit none
        double precision, intent(in)::rvec(:,:)
        double precision, intent(inout) :: dGopen 
        integer, intent(in) :: nb
        double precision :: r(NDIM), Uo(1), Uc(1), utmp(1)
        integer :: ib, cb ! Used for loop
        Uo = 0.0; Uc = 0.0; ! Initialize
        do ib = 1,nb
            r = rvec(ib,1:NDIM)
            call solvation_energyO(r,dGtransfer(ib), utmp(1)) !open 
            Uo = Uo + utmp; utmp=0.d0
            call solvation_energyC(r,dGtransfer(ib), utmp(1)) !closed
            Uc = Uc + utmp; utmp=0.d0
        enddo
        dGopen = dGopen + (Uo(1)-Uc(1)) ! Update

    end subroutine LGsolv
!=======================================================================================
! Subroutine solvation_energyO(r, g, energy) calculates the solvation energy 
!  at position r  for bead of hydrophobicity g in the open channel
!=======================================================================================

    subroutine solvation_energyO(r, g, etmp)
        implicit none
        double precision , intent (in) :: r(:), g
        double precision, intent (inout) :: etmp
        double precision :: Sx, dSdx, Sy, dSdy, rsym, alpha, fill

        etmp = 0.0
        ! check if at the membrane
        if (r(1) < RIGHTCUT .and. r(1) > LEFTCUT) then
            ! calculate the radial distance from y and z
            rsym = sqrt(r(2)*r(2)+r(3)*r(3))
            call switch_function(r(1), MBL, MBR, Sx, dSdx, 0.25d0)
            call switch_function(rsym, -rChannel, rChannel, Sy, dSdy, 0.25d0)
            Sy = 1-Sy;
            ! g-dependent term for membrane solvation
            ! and an in-channel dgfill_open term
            etmp = g*(Sx*Sy) + Sx * (DGFILLO) - (Sx*Sy)*(DGFILLO) 
        endif
     end subroutine solvation_energyO
!=======================================================================================
! Subroutine solvation_energyC(r, g, energy) calculates the solvation energy 
!  at position r for bead of hydrophobicity g in the closed channel
!=======================================================================================

    subroutine solvation_energyC(r, g, etmp)
        implicit none
        double precision , intent (in) :: r(:), g
        double precision, intent (inout) :: etmp
        double precision :: Sx, dSdx, Sy, dSdy, rsym, fill

        etmp = 0.0
        ! check if at the membrane
        if (r(1) < RIGHTCUT .and. r(1) > LEFTCUT) then
            ! calculate the radial distance from y and z
            rsym = sqrt(r(2)**2+r(3)**2)
            call switch_function(r(1), MBL, MBR, Sx, dSdx, 0.25d0)
            call switch_function(rsym, -rChannel, rChannel, Sy, dSdy, 0.25d0)
            Sy = 1-Sy;
            ! g-dependent term for membrane solvation
            ! and an in-channel dgfill_closed term
            etmp = g*(Sx*Sy) + Sx * (DGFILLC) - (Sx*Sy)*(DGFILLC)
        endif
     end subroutine solvation_energyC

!=========================================================================================
! LGopen(dGopen) controls the opening and closing of the lateral gate. Updates the 
! LGstatus if an opening or closing event takes place.
!=========================================================================================

    subroutine LGopen(dG, ranNum)

    implicit none
    double precision, intent(in) :: dG, ranNum! dG for opening, random number
    double precision :: prob, changeP ! weight of open, chance of change

    prob = exp(-dG/kBT)
    random_count = random_count + 1
    if (LGstatus.eq.0) then ! LG is currently closed, evaluate if it should open
        changeP = dt*LGfricti*prob/(1+prob)
        if (ranNum < changeP) then ! check if we change the status
            LGstatus = 1 ! Now open
        endif
    else
        changeP = dt*LGfricti/(1+prob)
        if (ranNum < changeP) then ! check if we change status
            LGstatus = 0 ! Now closed
        endif
    endif

    end subroutine LGopen

end module channel
