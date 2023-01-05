program workshare

    use OMP_LIB

    implicit none
    integer :: nthreads, tid, chunk, i
    integer, parameter :: n = 100
    integer, parameter :: chunksize = 10
    real :: a(n), b(n), c(n)

    ! initialization
    do i = 1, n
        a(i) = i * 1.0
        b(i) = a(i)
    enddo

    chunk = chunksize

!$omp parallel shared(a,b,c,nthreads,chunk) private(i,tid)

    tid = omp_get_thread_num()
    if (tid .eq. 0) then
        nthreads = omp_get_num_threads()
        print*, 'Number of threads = ', nthreads
    endif
    print*, 'Thread', tid, ' starting...'


! This parallelize the do loop
!$omp do schedule(dynamic,chunk)  
! This does not parallelize the do loop
!!$omp do ordered                   
    do i = 1, n
        c(i) = a(i) + b(i)
        print*, tid, i, c(i)
    enddo
!$omp end do nowait

    print*, 'Thread', tid, ' done.'

!$omp end parallel 

end program workshare

