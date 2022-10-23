program pi
    implicit none
    
    integer :: N
    real    :: area

    ! Initialize N, area
    N = 100000
    area = 0.

    call compute_integral(N, area)

    print *, 'PI = ', 2.*area

end program pi


subroutine compute_integral(N,A)

    implicit none
    integer, intent(in) :: N 
    real, intent(inout) :: A
    integer :: i 
    real    :: x, dx, h, dA
    real    :: my_func

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
end subroutine compute_integral


real function my_func(x)
    real :: x
    my_func = sqrt(1. - x**2.)
    return
end function