program ex1
    implicit none
    
    integer :: sum = 0
    integer :: amin = 1
    integer :: amax = 100
    integer :: n

    do n=amin, amax
        sum = sum + n
    end do

    print *, 'Summation from ', amin, 'to', amax, 'is', sum
    
end program ex1