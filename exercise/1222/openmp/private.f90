program private

    use OMP_LIB

    implicit none
    integer, parameter :: n = 10
    integer :: i, tmp

    tmp = 0

!$omp parallel do private(tmp) 
    do i = 1, n
        tmp = tmp + 1
    enddo
!$omp end parallel do

    print*, 'tmp = ', tmp

end program private



