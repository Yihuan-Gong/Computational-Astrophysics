! National Tsing Hua University
!
! ASTR 660 Computational Astrophysics
!
! Created:  Kuo-Chuan Pan 2020.04.14
! Modified: Karen Yang 2022.10.05
!
! Problem:
!
!        Solving non-linear equations
!
!
program nonlinear

    use solver
    implicit none

    integer, parameter  :: NMAX = 50    ! Max iteration number
    real, parameter     :: eps  = 1.e-6 ! tolerance

    integer  :: N              ! number of iteration
    real     :: error          ! error
    real     :: xs             ! solution
    character*40  :: fname     ! output filename
    real, external :: my_func  ! the function to solve
    real, external :: my_dfunc ! the d function to solve

    ! output file name
    ! fname = "bisection.txt"
    ! fname = "newton.txt"
    fname = "secant.txt"

    open(unit=1,file=trim(fname))
    ! write the header
    write(1,11) "#", "N", "solution","Error"

    ! The main iteration loop
    N     = 1
    error = 1e99

    do while (error > eps)
        ! call bisection(my_func, xs, error)
        ! call newton(my_func, my_dfunc, xs, error)
        call secant(my_func, xs, error)

        print *, "N = ",N, " solution is ", xs, " error = ", error
        write(1,12) N, xs, error

        N = N+1
        ! check if we have reached the maximum iteration number
        if (N .gt. NMAX) then
            print *, "The problem is not converged within ", N, " iterations."
            stop
        endif
    enddo

    close(1)

11  format(a2,a4,2a24)
12  format(2x,i4,2e24.14)


    print *, "Done!"
end program nonlinear 

!--------------------------------------------------
!
! my_function : The function to solve
!
!        f(x) = x^2 - 4 sin(x) = 0
!
!--------------------------------------------------
real function my_func(x)
    implicit none
    real, intent(in)  :: x

    my_func = x**2. - 4.*sin(x)

    return
end function my_func


! ---------------------------------------
! return f'(x) for Newton's method
! ---------------------------------------
real function my_dfunc(x)
    implicit none
    real :: x

    my_dfunc = 2.*x - 4.*cos(x)

    return
end function my_dfunc



