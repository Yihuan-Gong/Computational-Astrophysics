!---------------------------------------------------
!
! National Tsing Hua University
!
! ASTR 660 Computational Astrophysics
!
! Created:  Kuo-Chuan Pan 2020
! Modified: Karen Yang 2022.10.20
!
! Problem:
!
!        Simulating angry bird trajectories
!
program angry_bird
    use constants, only : show_constants, g ,pi, c
    use physics, only : update
    implicit none

    real  :: err=1e99            ! error
    real  :: xs                  ! solution
    real  :: fa, fx              ! f(a) and f(x)
    real  :: velocity = 30.
    real, save :: a = 10.        ! bracking interval [a,b]
    real, save :: b = 90.        ! bracking interval [a,b]
    
    do while (err > 0.1)

        ! The solution is the midpoint
        xs = (a+b)/2.

        ! f(a) and f(x)
        call trajectory(a,  velocity, fa)
        call trajectory(xs, velocity, fx)

        ! Error of the solution = f(x)
        err = abs(fx)

        ! Update the inteval [a,b]
        if (fx*fa > 0) then ! fa and fx have same sign
            a = xs 
        else
            b = xs
        end if

       print*, "Launch angle =", xs, "deg. The angry bird lands at x =", fx+50., "m. Error =", err, "m."

    end do

end program angry_bird




subroutine trajectory(angle,  velocity, f)
    use constants, only : show_constants, g ,pi, c
    use physics, only : update
    implicit none
    real, intent(in)  :: angle, velocity  ! Initial condtion
    real, intent(out) :: f                ! Landing position - 50m (target)  
    real :: dt, time, velx, vely, posx, posy
    real :: angle_rad
    character*40  :: fname     ! output filename

    angle_rad    = angle * pi / 180.0  ! change to rad
    velx = velocity * cos(angle_rad)
    vely = velocity * sin(angle_rad)

    dt   = 0.01  ! sec
    time = 0.0   ! initial time = 0.0

    posx = 0.0
    posy = 0.0


    do while (posy .ge. 0.0)
        call update(time, dt, posx, posy, velx, vely, posx, posy, velx, vely)
        time = time + dt
    end do

    f = posx - 50.

end subroutine trajectory



