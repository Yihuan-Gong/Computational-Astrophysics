program hello

    use OMP_LIB

    implicit none
    integer :: nthreads, tid

!! $omp parallel 
!! $omp parallel default(none)
!$omp parallel private(nthreads, tid)
    tid = omp_get_thread_num()
    print*, 'Hello World from thread = ', tid

    if (tid .eq. 0) then
      nthreads = omp_get_num_threads()
      print*, 'Number of threads = ', nthreads
    endif
!$omp end parallel 

end program hello
