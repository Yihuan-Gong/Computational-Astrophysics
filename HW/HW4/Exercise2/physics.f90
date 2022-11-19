!---------------------------------------------------
! The physics module
!
module physics
    use Simulation_data
    implicit none
    contains
        
        subroutine initial()
            !
            ! setup initial conditions of each stars
            ! in this example we only have two stars
            !
            use constants, only : au, msun, pi, G
            implicit none
            integer :: i
            real :: m1, m2, force

            m1 = 1.0 * msun
            m2 = 2.0 * msun

            !
            ! Use Kepler's law to evaluate the orbital period
            ! and use Newton's law to evaluate the force between
            ! two stars
            !

            separation = 3. * au !TODO
            period     = 2.*pi * sqrt(separation**3./(G*(m1+m2))) !TODO
            force      = G*m1*m2/(separation**2.) !TODO

            !
            ! setup initial conditions of star M1
            stars(1)%mass = m1
            stars(1)%x    =-2.*au           
            stars(1)%y    = 0.              
            stars(1)%vx   = 0.
            stars(1)%vy   = 2.*au * 2.*pi/period 
            stars(1)%ax   = force/m1
            stars(1)%ay   = 0.

            !
            ! setup initial conditions of star M2
            stars(2)%mass = m2
            stars(2)%x    = 1.*au
            stars(2)%y    = 0.
            stars(2)%vx   = 0.
            stars(2)%vy   =-1.*au * 2.*pi/period
            stars(2)%ax   = force/m2
            stars(2)%ay   = 0.

        end subroutine initial


        subroutine update(dt)
            use constants
            implicit none
            real, intent(in)  :: dt
            integer :: i,j
            real    :: x, y, rx, ry
            real    :: radius, mass

            !
            ! In this example we use a first order scheme (Euler method)
            ! we approximate dx/dt = v  --> x^(n+1) - x^n = v * dt
            ! therefore, x at step n+1 = x^(n+1) = x^n + v * dt
            !
            ! the same approximation can be applied to dv/dt = a
            !

            ! update position to t = t + dt
            stars(1)%x = stars(1)%x + stars(1)%vx * dt
            stars(1)%y = stars(1)%y + stars(1)%vy * dt
            stars(2)%x = stars(2)%x + stars(2)%vx * dt
            stars(2)%y = stars(2)%y + stars(2)%vy * dt

            ! update velocity to t = t + dt
            stars(1)%vx = stars(1)%vx + stars(1)%ax * dt
            stars(1)%vy = stars(1)%vy + stars(1)%ay * dt
            stars(2)%vx = stars(2)%vx + stars(2)%ax * dt
            stars(2)%vy = stars(2)%vy + stars(2)%ay * dt

            ! update accelerations to t = t + dt
            ! Star 1
            rx            = stars(1)%x - stars(2)%x  ! rx = star1 - star2 (source)
            ry            = stars(1)%y - stars(2)%y
            radius        = sqrt(rx**2. + ry**2.)
            mass          = stars(2)%mass   ! source mass (star 2 is the source of 
                                            ! gravity force for star 1)
            stars(1)%ax   = -G*mass*(rx)/(radius**3.)
            stars(1)%ay   = -G*mass*(ry)/(radius**3.)

            ! Star 2
            rx            = stars(2)%x - stars(1)%x 
            ry            = stars(2)%y - stars(1)%y
            radius        = sqrt(rx**2. + ry**2.)
            mass          = stars(1)%mass   ! source mass (star 1 is the source of 
                                            ! gravity force for star 2)
            stars(2)%ax   = -G*mass*(rx)/(radius**3.)
            stars(2)%ay   = -G*mass*(ry)/(radius**3.)


            return
        end subroutine update



        subroutine update_euler(dt)
            implicit none
            real, intent(in)  :: dt
            integer, parameter  :: n=4   ! number of ODEs
            real,dimension(n) :: y_s1, y_s2     ! y for star1 and star2
            real,dimension(n) :: k_s1, k_s2     ! dydt

            ! 
            ! Use Euler method here
            ! 

            ! Pack y(t) = [x, y, vx, vy](t)
            y_s1 =  (/stars(1)%x, stars(1)%y, stars(1)%vx, stars(1)%vy/)
            y_s2 =  (/stars(2)%x, stars(2)%y, stars(2)%vx, stars(2)%vy/)

            ! Update to y(t+dt)
            call my_func(n, y_s1, y_s2, k_s1, k_s2)
            y_s1 = y_s1 + dt*k_s1
            y_s2 = y_s2 + dt*k_s2

            ! Unpack y(t+dt) to star1 and star2
            stars(1)%x  = y_s1(1)
            stars(1)%y  = y_s1(2)
            stars(1)%vx = y_s1(3)
            stars(1)%vy = y_s1(4)
            stars(2)%x  = y_s2(1)
            stars(2)%y  = y_s2(2)
            stars(2)%vx = y_s2(3)
            stars(2)%vy = y_s2(4)

            return            
        end subroutine update_euler



        subroutine update_RK4(dt)
            implicit none
            real, intent(in)    :: dt
            integer, parameter  :: n=4   ! number of ODEs
            real,dimension(n) :: y_s1, y_s2     ! y for star1 and star2
            real,dimension(n) :: y2_s1, y2_s2
            real,dimension(n) :: y3_s1, y3_s2
            real,dimension(n) :: y4_s1, y4_s2    
            real,dimension(n) :: k1_s1, k1_s2
            real,dimension(n) :: k2_s1, k2_s2
            real,dimension(n) :: k3_s1, k3_s2
            real,dimension(n) :: k4_s1, k4_s2

            ! 
            ! Use RK4 method here
            ! 

            ! Pack y(t) = [x, y, vx, vy](t)
            y_s1 =  (/stars(1)%x, stars(1)%y, stars(1)%vx, stars(1)%vy/)
            y_s2 =  (/stars(2)%x, stars(2)%y, stars(2)%vx, stars(2)%vy/)

            ! Compute k1 = f(t, y)
            call my_func(n, y_s1, y_s2, k1_s1, k1_s2)

            ! compute y2 and k2
            y2_s1 = y_s1 + 0.5*dt*k1_s1
            y2_s2 = y_s2 + 0.5*dt*k1_s2
            call my_func(n, y2_s1, y2_s2, k2_s1, k2_s2)

            ! compute y3 and k3
            y3_s1 = y_s1 + 0.5*dt*k2_s1
            y3_s2 = y_s2 + 0.5*dt*k2_s2
            call my_func(n, y3_s1, y3_s2, k3_s1, k3_s2)

            ! compute y4 and k4
            y4_s1 = y_s1 + dt*k3_s1
            y4_s2 = y_s2 + dt*k3_s2
            call my_func(n, y4_s1, y4_s2, k4_s1, k4_s2)

            ! Update to y(t+dt)
            y_s1 = y_s1 + dt/6.*(k1_s1 + 2.*k2_s1 + 2.*k3_s1 + k4_s1)
            y_s2 = y_s2 + dt/6.*(k1_s2 + 2.*k2_s2 + 2.*k3_s2 + k4_s2)

            ! Unpack y(t+dt) to star1 and star2
            stars(1)%x  = y_s1(1)
            stars(1)%y  = y_s1(2)
            stars(1)%vx = y_s1(3)
            stars(1)%vy = y_s1(4)
            stars(2)%x  = y_s2(1)
            stars(2)%y  = y_s2(2)
            stars(2)%vx = y_s2(3)
            stars(2)%vy = y_s2(4)

            return            
        end subroutine update_RK4


        subroutine my_func(n, y_s1, y_s2, k_s1, k_s2)
            use constants, only : au, msun, pi, G
            implicit none
            real :: mass, rx, ry, radius
            integer, intent(in)  :: n   ! number of ODEs
            real,dimension(n),intent(in)  :: y_s1, y_s2     ! y    = [x,  y,  vx, vy]
            real,dimension(n),intent(out) :: k_s1, k_s2     ! dydt = [vx, vy, ax, ay]

            ! Velocity
            k_s1(1)  =  y_s1(3) ! Star1 vx
            k_s1(2)  =  y_s1(4) ! Star1 vy
            k_s2(1)  =  y_s2(3) ! Star2 vx
            k_s2(2)  =  y_s2(4) ! Star2 vy

            ! Accerleration
            ! Star 1
            rx            = y_s1(1) - y_s2(1)  ! x coor: star1 - star2(source)
            ry            = y_s1(2) - y_s2(2)  ! y coor: star1 - star2(source)
            radius        = sqrt(rx**2. + ry**2.)
            mass          = stars(2)%mass   ! source mass (star 2 is the source of 
                                            ! gravity force for star 1)
            k_s1(3)       = -G*mass*(rx)/(radius**3.)   ! ax
            k_s1(4)       = -G*mass*(ry)/(radius**3.)   ! ay

            ! Star 2
            rx            = y_s2(1) - y_s1(1)  
            ry            = y_s2(2) - y_s1(2)  
            radius        = sqrt(rx**2. + ry**2.)
            mass          = stars(1)%mass   ! source mass (star 1 is the source of 
                                            ! gravity force for star 2)
            k_s2(3)       = -G*mass*(rx)/(radius**3.)   ! ax
            k_s2(4)       = -G*mass*(ry)/(radius**3.)   ! ay

            return
        end subroutine my_func

end module physics

