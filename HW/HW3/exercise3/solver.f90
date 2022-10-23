module solver

    ! ---------------------------------------------------------
    ! ---------------------------------------------------------

    implicit none
    contains 

        !--------------------------------------------------
        !
        ! Implment the bisection method to solve the func
        !
        !
        ! Inputs:   func
        !
        ! Outputs:
        !           xs   : solution of func
        !           error: relative error
        !
        !--------------------------------------------------


        subroutine bisection(func, xs, err)
        implicit none
        real, external    :: func    ! the function to solve
        real, intent(out) :: xs      ! solution
        real, intent(out) :: err     ! error
        real, save :: a = -5.0        ! bracking interval [a,b]
        real, save :: b = -1.0        ! bracking interval [a,b]
        real  :: fa, fx              ! f(a) and f(x)
        
        ! The solution is the midpoint
        xs = (a+b)/2.
        
        ! f(a) and f(x)
        fa = func(a)
        fx = func(xs)
        
        ! Error of the solution = f(x)
        err = abs(fx)
        
        ! Update the inteval [a,b]
        if (fx*fa > 0) then ! fa and fx have same sign
            a = xs 
        else
            b = xs
        end if

        end subroutine bisection


        !--------------------------------------------------
        !
        ! Implment the Newton's method to solve my_func
        !
        !
        ! Inputs: func, dfunc
        !
        ! Outputs:
        !           xs   : solution of my_func
        !           error: relative error
        !
        !--------------------------------------------------
        subroutine newton(func, dfunc, xs, err)
        implicit none
        real, external    :: func  ! the function to solve
        real, external    :: dfunc ! the first derivative of the function to solve
        real, intent(out) :: xs    ! solution
        real, intent(out) :: err   ! error
        real, save :: x = -3.      ! trial value
        real :: fx
        real :: dfdx
        
        ! f(x) and f'(x)
        fx = func(x)
        dfdx = dfunc(x)
        
        ! Calculate x_{k+1} from x_k and the solution is x_{k+1}
        xs = x - fx/dfdx

        ! Calculate error of the solution
        err = abs(func(xs))

        ! Update x
        x = xs

        end subroutine newton


        !--------------------------------------------------
        !
        ! Implment the Secant method to solve my_func
        !
        !
        ! Inputs: None
        !
        ! Outputs:
        !           xs   : solution of my_func
        !           error: relative error
        !
        !--------------------------------------------------
        subroutine secant(func, xs, err)
        implicit none
        real, external    :: func  ! the function to solve
        real, intent(out) :: xs    ! solution
        real, intent(out) :: err   ! error

        real,save :: x0 = -5.0  ! initial guess
        real,save :: x1 = -3.0  ! initial guess

        real :: fx0, fx1       ! f(x0) and f(x1)
        real :: dfdx           ! slope of the scant line 
                               ! between (x0, f(x0)) and (x1, f(x1))

        fx0 = func(x0)
        fx1 = func(x1)
        dfdx = (fx0 - fx1)/(x0 - x1)
        
        ! Calculate x_{k+1} by x_k and x_{k-1}
        xs = x1 - fx1/dfdx

        ! Calculate error
        err = abs(func(xs))

        ! Update
        x0 = x1
        x1 = xs
        

end subroutine secant


end module solver
