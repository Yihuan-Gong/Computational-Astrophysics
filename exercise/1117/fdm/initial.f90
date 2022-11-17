subroutine initial()

    use Simulation_data
    implicit none

    integer :: i
    real, parameter :: small = 1.e-99

    dx = (xmax - xmin)/imax

    ! initialize the x array
    do i = istart-ibuf, iend+ibuf
        x(i) = xmin + (i-0.5)*dx
    enddo

    ! initialize u
    do i = istart, iend

        ! (a)
        if ((x(i) .ge. 0.1) .and. (x(i) .le. 0.2)) then
            u(i) = 1.0
        else
            u(i) = 0.01
        endif

    enddo

end subroutine initial
