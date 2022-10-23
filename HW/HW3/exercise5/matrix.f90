!
! National Tsing Hua University
!
! ASTR 660 Computational Astrophysics
!
! Created:  Kuo-Chuan Pan 2020
! Modified: Karen Yang 2022.10.15
!
! Problem:
!
!        Solving non-linear equations
!
program linear
    use linalg
    implicit none
    integer, parameter  :: N = 3
    real,dimension(N,N) :: lower, upper, A1, A2, P, Ainv
    real,dimension(N) :: b
    real,dimension(N) :: x
    real,dimension(4,4) :: aa,ll,uu,pp
    integer :: i,j

    ! Example in the lecture: PPT slide p.49 
    A1 = transpose(reshape((/  2.,  4., -2., &
                               4.,  9., -3., &
                              -2., -3.,  7. /), shape(A1)))   

    ! Exercise 4
    A2 = transpose(reshape((/  1.,  2., 2., &
                               4.,  6., 8., &
                               4.,  8., 10. /), shape(A2)))


    ! The vector b
    b(1) =  4.
    b(2) =  6.
    b(3) =  10.

    ! Call LU decomsition from matrix.f90
    call LU_decomposition(N, A2, lower, upper)

    ! Print the L and U matrix
    call mat_print("A", A2)
    call mat_print("L", lower)
    call mat_print("U", upper)

    ! Solve eq 4.1
    call solve_lu(N, A2, b, x)

    call mat_print("A",A2)
    print *, "vector   b = ",b
    print *, "solution x = ",x
    
end program linear 


