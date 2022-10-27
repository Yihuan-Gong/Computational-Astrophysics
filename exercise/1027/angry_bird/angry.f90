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

    real :: angle, velocity
    real :: dt, time, velx, vely, posx, posy
    real :: anal_y, vy0
    character*40  :: fname     ! output filename

    !call show_constants()

    angle    = 60.0                ! degree
    angle    = angle * pi / 180.0  ! change to rad

    velocity = 30 ! m/s 

    if (velocity .gt. c) then
        print *, "Error: the velocity cannot be faster than the speed of the light"
        stop
    endif

    velx = velocity * cos(angle)
    vely = velocity * sin(angle)

    dt   = 1.0  ! sec
    time = 0.0  ! initial time = 0.0

    posx = 0.0
    posy = 0.0

    anal_y = 0.0
    vy0    = vely

    ! Open a file to write the results
    fname = "angry_euler.txt"
    open(unit=1,file=trim(fname))
    ! write the header
    write(1,11) "#", "time", "posx","posy", "velx", "vely", "anal_y", "err_y"
    

    print *, "# time,  posx,  posy,  velx,  vely, anal_y, err_y"
    print *, time, posx, posy, velx, vely, anal_y, 0.0
    write(1,12), time, posx, posy, velx, vely, anal_y, 0.0



    do while (posy .ge. 0.0)

        call update(time, dt,posx,posy,velx,vely,posx,posy,velx,vely)
        time = time + dt
        anal_y = vy0*time - 0.5*g*time**2
        print *, time, posx, posy, velx, vely,anal_y,abs((anal_y-posy)/anal_y)
        write(1,12), time, posx, posy, velx, vely,anal_y,abs((anal_y-posy)/anal_y)
    end do

    close(1)

11  format(a2, 7a24)
12  format(2x,7e24.14)

end program angry_bird

