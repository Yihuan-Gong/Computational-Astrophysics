subroutine evolution()
    use Simulation_data
    use IO, only : output
    implicit none

    integer :: n
    integer :: interval
    real    :: dt, time


    n        = 0
    time     = 0.0

    dt = abs(dx/c)*cfl

    do while(time .le. tend)

        ! reset boundary condition

        call boundary(u)
        if (mod(n,io_interval) .eq. 0) then
            print *, "n =", n ," Time =", time
            call output(n,time)
        endif
        call update(time, dt)
        
        n    = n + 1
        time = time + dt
    enddo

end subroutine evolution


!!
!! 
!!
subroutine update(time, dt)
    use Simulation_data
    implicit none
    real, intent(in) :: time ,dt
    integer :: i
    real    :: FL, FR

    ! 1st order in time
    uold  = u
    do i = istart, iend
        call flux(i,dt,FL,FR)
        u(i) = uold(i) - dt/dx*(FR-FL)
    enddo

end subroutine update

!
! Routine to evalue flux the cell edge
!
subroutine flux(i,dt,FL, FR)
    use Simulation_data
    implicit none
    integer, intent(in) :: i
    real, intent(in)    :: dt
    real, intent(out)   :: FL, FR

    real :: sig, a, b, qL, qR


    ! Arithmetic average method
    !FL = 0.5*c*(uold(i-1)+uold(i))
    !FR = 0.5*c*(uold(i+1)+uold(i))

    ! The Lax-Friedrichs Method
    !FL = 0.5*c*(uold(i-1)+uold(i)) -0.5*dx/dt*(uold(i)-uold(i-1))
    !FR = 0.5*c*(uold(i+1)+uold(i)) -0.5*dx/dt*(uold(i+1)-uold(i))

    ! The upwind method
    !FL = ! TODO
    !FR = ! TODO

    ! The Lax-Wendroff Method
    !FL = ! TODO
    !FR = ! TODO

    !! Use piecewise linear and slope limiter

    !! left state
    call get_slope(dx,uold(i-2),uold(i-1),uold(i),sig) ! compute sig(i-1)
    FL = c*uold(i-1) + 0.5*c*(dx-c*dt)*sig

    !! right state
    call get_slope(dx,uold(i-1),uold(i),uold(i+1),sig) ! compute sig(i)
    FR = c*uold(i) + 0.5*c*(dx-c*dt)*sig

    return

end subroutine flux

subroutine get_slope(dx,l,m,r,sig)
    implicit none
    real, intent(in)  :: dx
    real, intent(in)  :: l   ! left 
    real, intent(in)  :: m   ! middle
    real, intent(in)  :: r   ! right
    real, intent(out) :: sig ! the slope
    real :: a, b

    ! compute a and b as the left/right slopes 
    a = (m-l)/dx
    b = (r-m)/dx

    ! TODO: implement the minmod limiter
    if ( a*b>0. ) then

        if ( (abs(a) < abs(b))  ) then
            sig = a
        else
            sig = b
        end if

    else
        sig = 0.
    end if

    return
end subroutine get_slope
