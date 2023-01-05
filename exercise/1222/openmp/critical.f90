program critical

    use OMP_LIB

    implicit none
    integer, parameter :: n = 100
    integer :: i, tmp

    tmp = 0

!$omp parallel do 
    do i = 1, n
        
        ! In order to fix the data race problem
        ! Only one thread can excute this line
        !$omp critical
        tmp = tmp + 1
        !$omp end critical

    enddo
!$omp end parallel do

    print*, 'tmp = ', tmp

end program critical


