program hello
    implicit none
    include 'mpif.h'
    integer :: rank, size ,ierror
    integer :: status(MPI_STATUS_SIZE)
    integer :: tag, i
    integer :: msg, msg_sum

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    msg_sum  = 0
    tag      = 1 
    msg      = rank

    if (rank .eq. 0) then
        msg_sum = msg
        do i = 1, (size-1)
            call MPI_Recv(msg, 1, MPI_INTEGER, i, tag, MPI_COMM_WORLD, status, ierror)
            msg_sum = msg_sum + msg
        enddo
    else
        call MPI_Send(msg, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, ierror)
    endif

    if (rank .eq. 0) print *, "sum = ", msg_sum
    call MPI_FINALIZE(ierror)

end program hello


