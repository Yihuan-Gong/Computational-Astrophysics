program test
    implicit none
    
    integer, parameter :: N=3
    integer, dimension(N,N) :: A

    A = transpose(reshape((/ 1, 2, 3, 4, 5, 6, 7, 8, 9 /), shape(array)))

end program test