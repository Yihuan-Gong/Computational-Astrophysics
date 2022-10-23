program pi
    implicit none
    
    integer :: i, N
    real    :: x, dx, h, dA, area
    real    :: my_func

    ! Initialize N, dx
    N = 100000
    dx = 2./real(N)

    area = 0.
    do i=1,N 
        ! Midpoint of rectangle i
        x = -1. + (real(i) - 0.5)*dx
        ! Height of rectangle i
        h = my_func(x)
        dA = dx*h 
        area = area + dA
    end do

    print *, 'PI = ', 2.*area

end program pi

real function my_func(x)
    real :: x
    my_func = sqrt(1. - x**2.)
    return
end function