program pi4
    implicit none
    
    integer :: i
    real    :: area=0, error
    real, parameter :: pi = 4.0*atan(1.0)
    integer, parameter :: NMAX = 7
    integer, dimension(NMAX) :: n_iteration
    n_iteration = (/10, 100, 1000, 10000, &
                    100000, 1000000, 10000000/)


    ! Open a file 'pi_error.dat'
    open(unit=10, file='pi_error_simpson_8.dat', &
         status='old', position='rewind')

    do i = 1, NMAX
        ! Call the subroutine
        call simpson(n_iteration(i), area)
        ! Compute the relative error
        error = (2*area-pi)/pi
        ! Write results to file
        write(10, *) error
    enddo

    ! Close the file
    close(unit=10)
    
end program pi4

! Midpoint method
subroutine midpoint(N,A)

    implicit none
    integer, intent(in) :: N 
    real, intent(inout) :: A
    integer :: i 
    real    :: x, dx, h, dA
    real    :: my_func

    A = 0
    dx = 2./real(N)

    do i=1,N 
        ! Midpoint of rectangle i
        x = -1. + (real(i) - 0.5)*dx
        ! Height of rectangle i
        h = my_func(x)
        dA = dx*h 
        A = A + dA
    end do

    return    
end subroutine midpoint

! Trapezoid method
subroutine trapezoid(N,A)

    implicit none
    integer, intent(in) :: N 
    real, intent(inout) :: A
    integer :: i 
    real    :: x, dx, dA, x0, x1
    real    :: my_func

    A = 0
    dx = 2./real(N)

    do i=1,N 
        ! Midpoint of rectangle i
        x = -1. + (real(i) - 0.5)*dx
        ! subinterval [x0, x1]
        x0 = x - 0.5*dx
        x1 = x + 0.5*dx

        ! The area within [x0, x1]
        if (i /= N) then
            dA = dx*(my_func(x0)+my_func(x1))/2.
        else
            dA = dx*my_func(x0)/2
        end if

        A = A + dA

    end do

    return    
end subroutine trapezoid

! Simpson method
subroutine simpson(N,A)

    implicit none
    integer, intent(in) :: N 
    real, intent(inout) :: A
    integer :: i 
    real    :: x, dx, h, dA, x0, x1
    real    :: my_func

    A = 0
    dx = 2./real(N)

    do i=1,N 
        ! Midpoint of rectangle i
        x = -1. + (real(i) - 0.5)*dx
        ! subinterval [x0, x1]
        x0 = x - 0.5*dx
        x1 = x + 0.5*dx
        ! The area within [x0, x1]
        if (i /= N) then
            dA = (dx/6.)* (my_func(x0) + 4.*my_func(x) + my_func(x1))
        else
            dA = (dx/6.)* (my_func(x0) + 4.*my_func(x))
        end if

        A = A + dA

    end do

    return    
end subroutine simpson


real function my_func(x)
    real :: x
    my_func = sqrt(1. - x**2.)
    return
end function