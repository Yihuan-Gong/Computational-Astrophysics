module Solver 

    implicit none
    contains

    subroutine euler(n, yin, ynext, t, h, func)
        implicit none
        integer, intent(in) :: n      ! number of ODEs
        real, intent(in)    :: t, h
        real, dimension(n), intent(in)  :: yin
        real, dimension(n), intent(out)  :: ynext
        external      :: func
        integer            :: i
        real, dimension(n) :: k1

        ! call func to obtain the values of dydt
        call func(n, t, yin, k1)

        ! compute ynext using the Euler's method
        ynext = yin + h*k1

    end subroutine euler

    subroutine rk2(n, yin, ynext, t, h, func)
        implicit none
        integer, intent(in) :: n      ! number of ODEs
        real, intent(in)    :: t, h
        real, dimension(n), intent(in)  :: yin
        real, dimension(n), intent(out)  :: ynext
        external      :: func
        integer            :: i
        real, dimension(n) :: k1, k2
        real,dimension(n)  :: y2

        ! compute k1 = func(t, yin)
        call func(n, t, yin, k1)

        ! compute y2 = yin + h*k1
        y2 = yin + h*k1

        ! compute k2 = func(t+h, y2)
        call func(n, t+h, y2, k2)

        ! compute ynext 
        ynext = yin + h/2. * (k1+k2)

    end subroutine rk2

    subroutine rk4(n, yin, ynext, t, h, func)
        implicit none
        integer, intent(in) :: n      ! number of ODEs
        real, intent(in)    :: t, h
        real, dimension(n), intent(in)  :: yin
        real, dimension(n), intent(out)  :: ynext
        external      :: func

        integer :: i
        real              :: h2
        real,dimension(n) :: k1,k2,k3,k4
        real,dimension(n) :: y2,y3,y4

        ! compute k1 = func(t, yin)
        call func(n, t, yin, k1)

        ! compute y2 and k2
        y2 = yin + 0.5*h*k1
        call func(n, t+0.5*h, y2, k2)

        ! compute y3 and k3
        y3 = yin + 0.5*h*k2
        call func(n, t+0.5*h, y3, k3)

        ! compute y4 and k4
        y4= yin + h*k3
        call func(n, t+h, y4, k4)

        ! compute ynext
        ynext = yin + h/6.*(k1 + 2.*k2 + 2.*k3 + k4)

    end subroutine rk4
end module Solver
