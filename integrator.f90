!============================================================================================
! Module containing the integrators that can be used to update position (and velocities).
!============================================================================================

module integrator

    use sys_param   ! Contains parameters for the functions used here
    use rand_num    ! random number generator

    implicit none
    double precision :: brownian_factor1, brownian_factor2
    ! Put whatever variables we still need here

    contains

!============================================================================================
! Subroutine to initialize some values that are used in the integrator, this saves time
! because we wont have to recalculate these every iteration.
!============================================================================================

    subroutine init_integrator()

        implicit none

        call init_random_seed() ! rand_num.f90, initialize the seed
        brownian_factor1 = dt*diffusion_constant/kBT
        brownian_factor2 = sqrt(2.0*diffusion_constant*dt)

    end subroutine init_integrator

!============================================================================================
! Subroutine for a Brownian integrator. It uses the parameters dt and diffusivity
! that are defined in sys_param. It also uses a random number generator that samplese a 
! Gaussian distribution with zero mean and unit variance to get variable eta.
!============================================================================================

    subroutine brownian_integrate(pos,force)

        ! Define input and output
        implicit none
        double precision, intent(in) :: force(freeStart:nFree,1:NDIM)
        double precision, intent(inout) :: pos(freeStart:nFree,1:NDIM)
        double precision :: eta
        integer :: ib,coord ! used for loops

        do ib = freeStart,nFree ! *** Here I move all beads, modify if some beads are not translocated
            do coord = 1,NDIM
                eta = gasdev()
                pos(ib,coord) = pos(ib,coord) + brownian_factor1*force(ib,coord) + brownian_factor2*eta
            enddo
        enddo

    end subroutine brownian_integrate

end module integrator
