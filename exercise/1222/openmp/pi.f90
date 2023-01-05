program pi

    use OMP_LIB

    implicit none
    integer :: i
    real    :: x, dx , sum, t1, t2
    integer,parameter  :: n = 10000000000

    dx = 1./float(n)

    t1 = omp_get_wtime()
    do i=0,n-1
        x = (i+0.5)*dx
        sum = sum + 4.0/(1.0+x*x)
    enddo
    t2 = omp_get_wtime()
    print *, "Pi   =", dx*sum
    print *, "Time =", (t2-t1)

end program pi


