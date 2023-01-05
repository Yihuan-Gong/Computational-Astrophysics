program hello
    use mpi
    implicit none
    integer :: rank, size ,ierror

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
    print *, "hello world from rank", rank, '/', size
    call MPI_FINALIZE(ierror)

end program hello


