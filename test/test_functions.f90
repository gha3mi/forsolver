module functions_module
use kinds

implicit none

contains

function F1(x) result(F_val)
   real(rk), intent(in) :: x
   real(rk) :: F_val
   F_val = 5_rk * x**3 + 8_rk * x - 5_rk
end function F1

function dF1dx(x) result(dFdx_val)
   real(rk), intent(in) :: x
   real(rk) :: dFdx_val
   dFdx_val = 15_rk * x**2 + 8_rk
end function dF1dx

function F2(x) result(F_val)
   complex(rk), intent(in) :: x
   complex(rk) :: F_val
   F_val = 5_rk * x**3 + 8_rk * x - 5_rk
end function F2

function F3(x) result(F_val)
   real(rk), dimension(:), intent(in) :: x
   real(rk), dimension(:), allocatable :: F_val

   allocate(F_val(2))

   F_val(1) = 2.0_rk*x(1) - 400.0_rk*x(1) * (x(2) - x(1)**2) - 2
   F_val(2) = 200.0_rk*x(2) - 200.0_rk*x(1)**2

end function F3

function dF3dx(x) result(dFdx_val)
   real(rk), dimension(:), intent(in) :: x
   real(rk), dimension(:,:), allocatable :: dFdx_val

   allocate(dFdx_val(2,2))

   dFdx_val(1,1) = 1200.0_rk*x(1)**2 - 400.0_rk*x(2) + 2
   dFdx_val(1,2) = - 400.0_rk*x(1)

   dFdx_val(2,1) = - 400.0_rk*x(1)
   dFdx_val(2,2) = 200.0_rk

end function dF3dx

function F4(x) result(F_val)
   complex(rk), dimension(:), intent(in)  :: x
   complex(rk), dimension(:), allocatable :: F_val
   integer :: i

   allocate(F_val(2))
   
   F_val(1) = 2.0_rk*x(1) - 400.0_rk*x(1) * (x(2) - x(1)**2) - 2
   F_val(2) = 200.0_rk*x(2) - 200.0_rk*x(1)**2

end function F4

end module functions_module
