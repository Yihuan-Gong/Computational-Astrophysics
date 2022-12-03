subroutine evolution()
    use Simulation_data
    use IO, only : output
    implicit none

    integer :: n
    real    :: dt, time

    n        = 0
    time     = 0.0

    dt = abs(dx/c)*cfl

    do while(time .le. tend)

        ! reset boundary condition
        call boundary(u)

        ! dump output times with frequency set by io_interval
        if (mod(n,io_interval) .eq. 0) then
            print *, "n =", n ," Time =", time
            call output(n,time)
        endif

        ! update the solution
        call update(time, dt)
        
        n    = n + 1
        time = time + dt
    enddo

end subroutine evolution


subroutine update(time, dt)
    use Simulation_data
    implicit none
    real, intent(in) :: time ,dt
    integer :: i
    real    :: FL, FR

    uold = u

    do i = istart, iend
        call flux(i, dt, FL, FR)
        u(i) = uold(i) - dt/dx*(FR-FL)
        
    enddo

end subroutine update


! 
! Routine to evaluate flux at cell edges
! 
subroutine flux(i, dt, FL, FR)
    use Simulation_data
    implicit none
    integer, intent(in) :: i
    real, intent(in)  :: dt
    real, intent(out) :: FL, FR

    ! ! Arithmetic average method
    ! FL = 0.5*(c*uold(i-1) + c*uold(i))
    ! FR = 0.5*(c*uold(i+1) + c*uold(i))

    ! ! ! The Lax-Friedrich method
    FL = 0.5*(c*uold(i-1) + c*uold(i)) - dx/(2.*dt)*(uold(i)-uold(i-1))
    FR = 0.5*(c*uold(i) + c*uold(i+1)) - dx/(2.*dt)*(uold(i+1)-uold(i))

    ! ! ! Upwind method
    ! FL = c*uold(i-1)
    ! FR = c*uold(i)

    ! ! Lax-Wendroff method
    ! FL = c*(0.5*(uold(i) + uold(i-1)) - dt/(2*dx)*(c*uold(i)-c*uold(i-1)))
    ! FR = c*(0.5*(uold(i+1) + uold(i)) - dt/(2*dx)*(c*uold(i+1)-c*uold(i)))

    
end subroutine flux