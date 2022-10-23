module linalg

    ! ---------------------------------------------------------
    ! ref: https://rosettacode.org/wiki/LU_decomposition#Fortran
    ! ---------------------------------------------------------

    implicit none
    contains 


        subroutine mat_print(amsg,a)
            character(*), intent(in) :: amsg
            class    (*), intent(in) :: a(:,:)
            integer                  :: i
            print*,' '
            print*,amsg
            do i=1,size(a,1)
                select type (a)
                    type is (real(8)) ; print'(100f8.3)',a(i,:)
                    type is (integer) ; print'(100i8  )',a(i,:)
                end select
            end do
            print*,' '
        end subroutine

        subroutine solve_lower_triangular_matrix(N,L,b,x)
            ! Solve Lx = b
            implicit none

            integer, intent(in)  :: N
            real, dimension(N,N), intent(in)  :: L  ! lower triangle
            real, dimension(N)  , intent(in)  :: b  ! vector
            real, dimension(N)  , intent(out) :: x  ! solution
            real, dimension(N)  :: bs              

            integer :: i,j

            bs = b ! Replace b by bs, otherwise the compiler would complain
            do j = 1, N
                if ( L(j,j) == 0 ) then
                    stop
                end if
                x(j) = bs(j)/L(j,j)
                do i = j+1, N
                    bs(i) = bs(i) - L(i,j)*x(j)
                end do
            end do

            return
        end subroutine solve_lower_triangular_matrix

        subroutine solve_upper_triangular_matrix(N,U,b,x)
            implicit none

            integer, intent(in)  :: N
            real, dimension(N,N), intent(in)  :: U  ! upper triangle
            real, dimension(N)  , intent(in)  :: b  ! vector
            real, dimension(N)  , intent(out) :: x  ! solution
            real, dimension(N)  :: bs

            integer :: i,j

            bs = b ! Replace b by bs, otherwise the compiler would complain
            do j = N, 1, -1 
                if ( U(j,j) == 0 ) then
                    stop
                end if
                x(j) = bs(j)/U(j,j)
                do i = 1, j-1
                    bs(i) = bs(i) - U(i,j)*x(j)
                end do
            end do
            
            return
        end subroutine solve_upper_triangular_matrix

        subroutine LU_decomposition(N,A,L,U)
            implicit none
            
            integer, intent(in)  :: N
            real, dimension(N,N), intent(in)  :: A    ! matrix
            real, dimension(N,N), intent(out)    :: L,U  ! matrix
            real, dimension(N,N) :: M, As
            
            integer :: i,j,k

            ! Replace A as As, otherwise the compiler would complain
            As = A

            ! Initialize M
            M = transpose(reshape((/  0.,  0., 0., &
                                      0.,  0., 0., &
                                      0.,  0., 0. /), shape(M)))

            ! Obtain matrix M and A
            do k = 1, N-1
                if ( As(k,k)==0 ) then
                    stop
                end if
                do i = k+1, N
                    M(i,k) = As(i,k)/As(k,k)
                end do
                do j = k+1, N
                    do i = k+1, N
                        As(i,j) = As(i,j) - M(i,k)*As(k,j)
                    end do
                end do
            end do

            ! L = resulting matrix M + identity matrix I
            do j = 1, N
                do i = 1, N
                    if ( i==j ) then
                        L(i,j) = M(i,j) + 1.
                    else
                        L(i,j) = M(i,j)
                    end if                  
                end do
            end do

            ! U = upper triangular part of resulting matrix A
            do j = 1, N
                do i = 1, N
                    ! upper triangular matrix: column index >= row index
                    if ( j>=i ) then 
                        U(i,j) = As(i,j)
                    else
                        U(i,j) = 0.
                    end if                    
                end do
            end do

        end subroutine LU_decomposition

        subroutine solve_lu(N,A,b,x)
            ! solve: A x = b
            integer, intent(in)  :: N
            real, dimension(N,N), intent(in)  :: A  ! upper triangle
            real, dimension(N)  , intent(in)  :: b  ! vector
            real, dimension(N)  , intent(out) :: x  ! solution

            real, dimension(N,N) :: L, U, P
            real, dimension(N)   :: y, pb

            ! Calculate L and U matrix
            call LU_decomposition(N, A, L, U)

            ! Solve Ly = b
            call solve_lower_triangular_matrix(N, L, b, y)

            ! Solve Ux = y
            call solve_upper_triangular_matrix(N, U, y, x)

        end subroutine solve_lu

end module linalg
