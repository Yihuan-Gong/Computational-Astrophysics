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
    real,dimension(N,N) :: lower, upper, A, P, Ainv
    real,dimension(N) :: b
    real,dimension(N) :: x
    real,dimension(4,4) :: aa,ll,uu,pp
    integer :: i,j

    ! lower triangle
    lower(1,1) = -1.0
    lower(1,2) =  0.0
    lower(1,3) =  0.0

    lower(2,1) = -6.0
    lower(2,2) = -4.0
    lower(2,3) =  0.0

    lower(3,1) =  1.0
    lower(3,2) =  2.0
    lower(3,3) =  2.0

    ! the vectore b
    b(1) =  1.0
    b(2) = -6.0
    b(3) =  3.0

    call solve_lower_triangular_matrix(N,lower,b,x)

    call mat_print("L",lower)
    print *, "vector   b = ",b
    print *, "solution x = ",x
    
end program linear 


