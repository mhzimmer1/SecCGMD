module rand_num 

  ! random number generator
    use sys_param, only : random_count, seed, seed_size, extraRanSeed
    implicit none
  ! gaussian random number
    integer :: iset = 0
    double precision :: gset

    contains
! ====================================== !
! initialize the random number generator
    subroutine init_random_seed()
        
        implicit none
        integer :: i, clock
        double precision :: ran
          
        call RANDOM_SEED(size = seed_size)
        allocate(seed(seed_size))
          
        call SYSTEM_CLOCK(COUNT=clock)

        clock = clock + extraRanSeed*100
        seed = clock + 37 * (/ (i - 1, i = 1, seed_size) /) ! Change 12 to clock if you want to go back to random!!!
        ! if I submit many jobs really fast, they might all have the
        ! same time. That's why we need the extraRanSeed from the input

        call RANDOM_SEED(PUT = seed)

        random_count = 0
        call RANDOM_NUMBER(ran)
        random_count = random_count + 1
        print *, "random number check: ", ran
          
    end subroutine
      
! ======================================== !
! initialize random number from restart
    subroutine init_random_seed_restart(fid)
        implicit none
        integer :: i, error, fid
        double precision :: ran
        integer(8) :: ic

        read(fid,*,iostat=error) seed_size
        allocate(seed(seed_size))
        i = 1
        do while(.true.)
            read(fid, *, iostat=error) seed(i)
            if (error /= 0) exit
            i = i + 1
        enddo ! done read file
        if (i .ne. seed_size+1) then
            print *, "error reading seeds, !!",i, seed_size
            stop
        endif
        print *, seed
        call RANDOM_SEED(PUT = seed)
        ic = 1.0 ! sometimes random_count can be really big, need the double
                 ! precision counter
        read(fid,*,iostat=error) random_count
        print *, "starting random sequence at: ", random_count
        do while (.true.)
            if (ic > random_count) then
                exit
            endif
            call RANDOM_NUMBER(ran)
            ic = ic+1
            if ( modulo(ic,100000000) .lt. 1e-5) then
                print *, ic
            endif
        enddo
        print *, "finished initializing random number!"
        deallocate(seed)

        read(fid,*) iset
        read(fid,*) gset

    end subroutine
! ======================================== !
! gaussian random number generator
    function gasdev()
        implicit none
        double precision :: gasdev, ran(2)
        double precision :: fac,rsq,v1,v2,ran1

        if (iset == 0) then
1               call Random_Number(ran)
            random_count = random_count + 2
            v1=2.d0*ran(1)-1.d0
            v2=2.d0*ran(2)-1.d0
            rsq=v1**2+v2**2
            if(rsq >= 1. .or. rsq == 0.)goto 1
            fac=sqrt(-2.d0*dlog(rsq)/rsq)
            gset=v1*fac
            gasdev=v2*fac
            iset=1
        else
            gasdev=gset
            iset=0
        endif

        return
    end function
      
end module

