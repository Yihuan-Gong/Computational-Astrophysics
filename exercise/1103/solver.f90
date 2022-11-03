module Solver 

    implicit none
    contains

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

end module Solver
