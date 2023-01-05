program hello
    implicit none
    include 'mpif.h'
    integer :: rank, size ,ierror
    integer :: status(MPI_STATUS_SIZE)
    integer :: tag, i
    integer :: sendbuf, recvbuf, targetRank

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    tag      = 1 
    sendbuf    = (rank+1)*10
    targetRank = mod(rank+1,2)

    ! rank 0: send first and then receive using blocking transfer
    if (rank .eq. 0) then
        call MPI_Send(sendbuf, 1, MPI_INTEGER, targetRank, tag, &
                      MPI_COMM_WORLD, ierror)
        call MPI_Recv(recvbuf, 1, MPI_INTEGER, targetRank, tag, &
                      MPI_COMM_WORLD, status, ierror) 
    ! rank 1: receive first and then send using blocking transfer
    else
        call MPI_Recv(recvbuf, 1, MPI_INTEGER, targetRank, tag, &
                      MPI_COMM_WORLD, status, ierror) 
        call MPI_Send(sendbuf, 1, MPI_INTEGER, targetRank, tag, & 
                      MPI_COMM_WORLD, ierror)
    endif

    print*, "Rank ", rank, ": Send ", sendbuf, ", Recv ", recvbuf 
    call MPI_FINALIZE(ierror)

end program hello


