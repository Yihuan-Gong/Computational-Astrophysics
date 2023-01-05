program hello

    implicit none

    !$omp parallel
    print *, " hello world"
    !$omp end parallel

end program

