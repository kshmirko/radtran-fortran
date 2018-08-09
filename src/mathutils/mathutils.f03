module mathutils

  implicit none
  integer, parameter :: dp = kind(1.d0)
  real(kind=dp), parameter  ::  pi=acos(-1.0_dp)
contains
  
  pure function trapz(x, y) result(r)
    !! Calculates the integral of an array y with respect to x using the trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be uniform.
    real(kind=dp), intent(in)  :: x(:)         !! Variable x
    real(kind=dp), intent(in)  :: y(size(x))   !! Function y(x)
    real(kind=dp)              :: r            !! Integral ∫y(x)·dx

    ! Integrate using the trapezoidal rule
    associate(n => size(x))
      r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
    end associate
  end function trapz
  
  subroutine linspace(x0, x1, y)
    real(kind=dp), intent(in) ::  x0, x1
    integer       ::  N, I
    real(kind=dp) ::  y(:)
    real(kind=dp) ::  dx
    
    N = size(y, 1)
    
    dx = (x1-x0)/dble(N-1)
    
    DO I=1, N
      y(I) = x0+I*dx
    END DO
    return 
  end subroutine linspace
  
  
  
end module mathutils
