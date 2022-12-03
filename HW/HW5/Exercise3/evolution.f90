subroutine evolution()
    use Simulation_data
    use IO, only : output
    implicit none

    integer :: n
    integer :: interval
    real    :: dt, time


    n        = 0
    time     = 0.0

    dt = abs(dx/cx)*cfl

    do while(time .le. tend)

        ! reset boundary condition

        ! call boundary(u)
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
    integer :: i,j
    real    :: FL, FR ! FL: left  FR: right
    real    :: FT, FB ! FT: top   FB: bottom


    ! Do x direction first
    uold  = u             ! copy the solution to uold array
    call boundary_xdir(u) ! update BC in x

    do j = istart, iend
        do i = istart, iend
            call flux_xdir(i,j,dt,FL,FR)
            u(i,j) = uold(i,j) - dt/dx*(FR-FL)
        enddo
    end do

    ! Then repeat for y direction
    uold  = u             ! copy the solution to uold array
    call boundary_ydir(u) ! update BC in x

    do j = istart, iend
        do i = istart, iend
            call flux_ydir(i,j,dt,FB,FT)
            u(i,j) = uold(i,j) - dt/dy*(FT-FB)
        enddo
    end do
    

end subroutine update

!
! Routine to evalue flux the cell edge
!
subroutine flux_xdir(i, j, dt, FL, FR)
    use Simulation_data
    implicit none
    integer, intent(in) :: i, j
    real, intent(in)    :: dt
    real, intent(out)   :: FL, FR

    real :: sig, a, b, qL, qR

    !! Use piecewise linear and slope limiter

    !! left state
    call get_slope(dx,uold(i-2,j),uold(i-1,j),uold(i,j),sig) ! compute sig(i-1)
    FL = cx*uold(i-1,j) + 0.5*cx*(dx-cx*dt)*sig

    !! right state
    call get_slope(dx,uold(i-1,j),uold(i,j),uold(i+1,j),sig) ! compute sig(i)
    FR = cx*uold(i,j) + 0.5*cx*(dx-cx*dt)*sig

    return

end subroutine flux_xdir

subroutine flux_ydir(i, j, dt, FB, FT)
    use Simulation_data
    implicit none
    integer, intent(in) :: i, j
    real, intent(in)    :: dt
    real, intent(out)   :: FB, FT

    real :: sig, a, b, qL, qR


    !! bottom state
    call get_slope(dx,uold(i,j-2),uold(i,j-1),uold(i,j),sig) ! compute sig(i-1)
    FB = cy*uold(i,j-1) + 0.5*cy*(dy-cy*dt)*sig

    !! top state
    call get_slope(dx,uold(i,j-1),uold(i,j),uold(i,j+1),sig) ! compute sig(i)
    FT = cy*uold(i,j) + 0.5*cy*(dy-cy*dt)*sig

    return

end subroutine flux_ydir

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
