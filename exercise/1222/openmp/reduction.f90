program reduction

    use OMP_LIB

    implicit none
    integer :: i
    integer, parameter :: n = 100
    real :: a(n), sum

    ! initialization
    do i = 1, n
        a(i) = i * 1.0
    enddo
    sum = 0.0

!$omp parallel do reduction(+:sum)
    do i = 1, n
        sum = sum + a(i)
    enddo

    print*, 'sum = ', sum

end program reduction

