program pi
    implicit none
    
    integer :: N
    real    :: area, error
    real, parameter :: pi= 4.0*atan(1.0)
    integer, parameter :: NMAX = 8
    integer, dimension(NMAX)
    n_iteration = (/10, 100, 10000, &
                    100000, 1000000, 10000000/)


    ! Open a file 'pi_error.dat'
    open(10, 'pi_error.dat')

    do i = 1,NMAX
        ! Call the subroutine
        call compute_integral(N, area)
        ! Compute the relative error
        error = (2*area-pi)/pi
        ! Write results to file
        write(10, *) error
    enddo

    ! Close the file
    close(unit=10)
    
end program pi


subroutine compute_integral(N,A)

    implicit none
    integer, intent(in) :: N 
    real, intent(inout) :: A
    integer :: i 
    real    :: x, dx, h, dA
    real    :: my_func

    dx = 2./real(N)
    A = 0

    do i=1,N 
        ! Midpoint of rectangle i
        x = -1. + (real(i) - 0.5)*dx
        ! Height of rectangle i
        h = my_func(x)
        dA = dx*h 
        A = A + dA
    end do

    return    
end subroutine compute_integral


real function my_func(x)
    real :: x
    my_func = sqrt(1. - x**2.)
    return
end function