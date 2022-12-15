function interpolate(a, n, x, xx) 

  implicit none

  real, dimension(n), intent(IN) :: a
  integer, intent(IN) :: n
  real, dimension(n), intent(IN) :: x
  real, intent(IN) :: xx

  integer :: i
  real :: r, interpolate

  r = (log10(xx) - log10(x(1)))*n/(log10(x(n))-log10(x(1)))
  i = aint(r)
  if (a(i) > 0.0) then
     interpolate = a(i)*(a(i+1)/a(i))**(r-i)
  else
     interpolate = 0.0
  endif

  return

end function interpolate
