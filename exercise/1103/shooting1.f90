!---------------------------------------------------
!
! National Tsing Hua University
!
! ASTR 660 Computational Astrophysics
!
! Created:  Kuo-Chuan Pan 2020
! Modified: Karen Yang 2022.10.25
!
! Problem:
!
!        Solving boundary value problems

program shooting1

    use Solver , only : rk2
    implicit none

    external :: my_func

    real :: h,t,tend
    integer :: i,n
    real :: try_y2
    real, dimension(2) :: y, ynext

    h    = 0.01  ! step size
    t    = 0.0   ! initial t
    tend = 1.0   ! final   t

    ! a trial value
    try_y2 = -1.0

    n    =  2   ! number of ODEs

    ! initial conditions
    y(1) =  1.0      ! y(1) = u
    y(2) =  try_y2   ! y(2) = u' 

    do while(t .lt. tend)
        if ((t+h) .ge. tend) then
            h = tend - t
        endif
        call rk2(n, y, ynext, t, h, my_func)
        y = ynext
        t = t + h
        print *, t, ynext(1)
    enddo
    print *, "y1(1) = ", ynext(1)
    print *, "The desired value is 1."

end program shooting1

!-----------------------------------------------
!
!  Solving Boundary Value problem
!
!  u'' = 6t      0 < t < 1
!
!  with BC
!
!         u(t=0) = 1 and
!         u(t=1) = 1
!
!-----------------------------------------------

subroutine my_func(n, t, yin, k)
    implicit none
    integer, intent(in) :: n  ! number of ODEs
    real, intent(in)    :: t
    real, dimension(n), intent(in)  :: yin
    real, dimension(n), intent(out) :: k    ! dydt

    ! in this example  n = 2
    k(1) =  yin(2)
    k(2) =  6.*t
    return
end subroutine my_func







