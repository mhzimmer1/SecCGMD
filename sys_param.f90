!=========================================================================================
! This module contains all the parameters used throughout the program, only 
! parameters are included. Variables are defined where they are required.
!=========================================================================================

module sys_param

    integer(8) :: MAXSTEPS=20e12 ! Terminate the simulation if it does not complete in this number of steps
    ! Global parameters, from the Cell Reports paper
    integer, parameter :: NDIM = 3, x=1, y=2, z=3   ! Define the dimensionality
    double precision, parameter :: kBT = 1.0    ! Energy in reduced units
    double precision :: beta = 1.0/kBT  ! 1/kBT 
    integer :: ForcePulling = 0 ! Are we running the force pulling protocol?
    double precision, parameter :: pi = dacos(-1.d0)
    double precision, parameter :: T = 300.0    ! Temperature
! *** timestep 3 times larger than previously, this is acceptable because the Diffusion coefficient is 3x smaller
    double precision, parameter :: dt = 3e4, dtRealTime = 300.0 ! I did this so I can change the timestep during equilibration etc
                            ! 100ns *** in the 2D code, have to determine what it 
                            ! really corresponds to here by checking diffusion through SecY
    double precision :: machine_eps=EPSILON(kBT)    ! Machine precision
    ! Related to the polymer, these should be overwritten when the polymer is initialized
    integer,parameter :: MAXNUMBEADS=500    ! I can probably remove this when the init_polymer subroutine is ready
    integer :: nBeads = 55      ! Number of beads in the polymer
    integer :: nInit=1, nFree=1 ! Initial beads, Number of free beads
    integer :: freeStart=1      ! Beads before this bead are not allowed to move (f.e. Stop-transfer Cell Reports)
    integer :: nForce=2         ! Number of beads on which we calculate forces
    integer :: trSteps=1250000  ! Number of steps to translate 1 bead, corresponds to 24res/sec *** in 2D code
    
    ! Related to the polymer and the different bead types
    integer :: nCharge=0 ! Number of hydorphilic, hydrophobic and charged beads
    integer :: chBeads(MAXNUMBEADS) ! Used to store index of charged beads 
    ! Related to FENE bond potential
    double precision, parameter :: feneK = 5.833 , feneR2 = 4.0 ! Ro is only used squared
    ! Related to BiP binding
    double precision :: k_bip_on = 1e-11, k_bip_off = 0.167e-12 ! On and off rates, Bin
    double precision :: F_bip = 0.39 ! Force along channel axis due to Bip Binding, Bin converted to eps/sigma
    double precision, parameter :: bipX = 2.0 ! Force acts for beads within 2sigma of (bipX,0,0), Bin
    integer :: BIPON=0 ! Determines whether there is explicit Bip binding in the simulation, 1 is yes, 0 is no
    ! Related to the transmembrane potential
    double precision :: tmpU = 0, & ! -3.74, & ! (E coli value) The magnitude of the transmembrane potential (eps/e) 
                        tmpL = 1.6, & ! Lengthscale of the TM potential (1/sigma)
                        tmpUL = 0 ! -5.98 ! (E coli value) Product of the previous, prefactor
    ! *** Added, related to NC-NC LJ potential
    double precision :: ljEpsA = 1.d0, ljRcrA = 2.5d0, ljRcrA2 = 6.25d0 ! Epsilon and cut-off for the attractive NC-NC interaction
    ! Related to LJ potential
    double precision, parameter :: ljEps = 1.0, ljEpsLG = 1.5   ! For these parameters
    double precision, parameter :: ljRcl = 0.5, ljRclLG = 0.1   ! seperate values are used
    double precision, parameter :: ljRcrLG = 2.5            ! normal, for LG
    double precision :: ljRcr2 = 1.0, ljLGcr2 = 1.0, ljLGcl2 = 1.0 
                                    ! I use squared values for cutoffs
                                    ! because the magnitude of the vector
                                    ! between particles is always calculated
                                    ! as r^2
    double precision :: alphaLJ=0.0 ! Soft-ness parameter in the LJ potential; U=4eps[1/(a+r6)2-1/(a+r6)]
    double precision :: maxE=100.0 ! Maximum repulsive energy in the LJ potential
    ! Related to Electrostatic potential ! **** I CHANGED KAPPA WRT BIN 1
    ! instead of 1.4 ***
    double precision, parameter :: dhKappai = 1.0, dhRclLG = 0.5, dhRcrLG = 2.5, &  ! cutoff for open LG
                                    dhcl2 = 0.25, dhRclLG2 = 0.25, dhRcrLG2 = 6.25
    ! Related to the ribosome
    double precision, parameter :: exitX = -6.4d0, exitY = 0.d0, exitZ = 1.4d0 ! Ribosome exittunnel
    double precision, parameter :: riboX = -10.d0, riboY = -0.5, riboZ = 1.0 ! Ribosome centroid
    double precision :: riboCenter(NDIM) = (/riboX, riboY, riboZ/)
    double precision, parameter :: riboR = 2.0 ! Radius of the ribosome
    ! Related to the solvation
    double precision :: LEFTCUT=-4.0, RIGHTCUT=4.0,  & ! *** IN THIS VERSION NOT PARAMETERS BUT INPUTS
                       MBL=-2.0, MBR=2.0, CUTB=0.25,    & ! Parameters in switch function
                       rChannel=1.5 ! Radius of the channel, used for determining if beads are in the channel too
  !  double precision, parameter :: LEFTCUT=-4.0, RIGHTCUT=4.0,  & ! Cutoffs for solvation potential
  !                     MBL=-2.0, MBR=2.0, CUTB=0.25,    & ! Parameters in switch function
  !                     rChannel=1.5 ! Radius of the channel, used for determining if beads are in the channel too
    integer, parameter :: NGRID = 8000 ! Number of gridpoints used for the tabulated tanh function
    ! Related to the geometry of the SecY channel, for more information also check the .ppt files I have on this
    integer :: nBeadChan=48 ! Number of beads in the channel

    ! Will contain channel data read from file
    double precision, allocatable :: chanBeads(:,:), chanCharge(:), chanP(:,:)
    character(1), allocatable :: chanType(:)
    double precision, parameter :: rOpening=2.d0, hChannel=2.d0 ! Relate to closed channel geometry
    double precision, parameter :: tauLG=500.0 ! Time for LG opening in nanoseconds
    ! WHEN CHANGING THE CHANNEL GEOMETRY MAKE SURE TO AVOID TOO STRONG CURVATURE (CHECK YOUR PARAMETERS IN closedC.m)
    integer :: NGRIDX = 91, NGRIDY =  91, NGRIDZ = 91 ! Number of gridpoints for the translocon, recommended at least
    double precision, parameter :: DGRID = 0.1 ! Spacing for the force grid
    double precision :: DGRIDi = 1.d0/DGRID ! Inverse of the spacing, used for bin indexing
                                ! one point per 0.1 sigma??
    double precision :: xCutChan=hChannel+ljRcrLG, yzCutChan=rOpening+ljRcrLG ! Define the range in which particles feel SecY
    ! Cut-offs for channel from file
    double precision :: xCutChanL, xCutChanR
    double precision :: yCutChanL, yCutChanR
    double precision :: zCutChanL, zCutChanR
    ! Flag to track if the translocon is currently open or closed
    integer :: LGstatus = 0 ! 0 is closed, 1 is open
    ! Related to the integrator
! *** THIS VALUE OF 1.2e-8 IS TOO HIGH, I USE IT FOR FASTER TESTING
    double precision, parameter :: diffusion_constant=0.4e-8    ! *** 3D VALUE~1/3 of 2D *** Units are here sigma^2/0.01ns
    !double precision, parameter :: diffusion_constant=1.1855e-8    ! Not rounded value, makes no difference
    ! Related to random numer generator (I copied the one Bin uses, it works fine)
    integer(8) :: random_count=0
    integer, dimension(:), allocatable :: seed
    integer :: seed_size
    integer :: extraRanSeed=1
    ! Flag for restart
    integer :: restart=0
    
    ! File IDs
    integer, parameter ::   polyFID     = 11,   &
                trajFID     = 12,   &
                restartFID  = 13,   &
                debugFID        = 14,   &
                trajDFID        = 16,   &  
                energyFID       = 17,   &  
                analysisFID     = 15    ! These fileIDs are reserved. For
                ! variable files called only once
                ! I will only use the ID 77
    !****Connie LG opening parameters ****!
    double precision :: DGEMPTY=3.0 , DGFILLC=0.0, DGFILLO=0.0  ! modifications to PE surfaces
    integer :: FIXEDLG = 0! for testing fixed LG simulations
end module sys_param
