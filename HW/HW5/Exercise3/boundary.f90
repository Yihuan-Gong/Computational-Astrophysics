subroutine boundary_xdir(v)
    !
    ! apply BC on array v
    !
    use Simulation_data
    implicit none
    real, dimension(istart-ibuf:iend+ibuf, istart-ibuf:iend+ibuf), intent(inout) :: v
    integer :: i,j


    ! apply boundary condition

    do j = istart, iend
        ! left boundary (period)
        do i = 1, ibuf
            v(istart-i, j) = v(iend-i+1, j)
        enddo
        ! right boundary (period)
        do i = 1, ibuf
            v(iend+i, j) = v(istart+i-1, j)
        enddo
    end do
    

end subroutine boundary_xdir


subroutine boundary_ydir(v)
    !
    ! apply BC on array v
    !
    use Simulation_data
    implicit none
    real, dimension(istart-ibuf:iend+ibuf, istart-ibuf:iend+ibuf), intent(inout) :: v
    integer :: i,j


    ! apply boundary condition

    do i = istart, iend
        ! bottom boundary (period)
        do j = 1, ibuf
            v(i, istart-j) = v(i, iend-j+1)
        enddo
        ! top boundary (period)
        do j = 1, ibuf
            v(i, iend+j) = v(i, istart+j-1)
        enddo
    end do
    

end subroutine boundary_ydir