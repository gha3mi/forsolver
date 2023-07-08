module functions_module
use kinds
implicit none

contains

function F(x) result(F_val)
   real(rk), intent(in) :: x
   real(rk) :: F_val
   F_val = 5_rk * x**3 + 8_rk * x - 5_rk
end function F

function dFdx(x) result(dFdx_val)
   real(rk), intent(in) :: x
   real(rk) :: dFdx_val
   dFdx_val = 15_rk * x**2 + 8_rk
end function dFdx

end module functions_module



program main
   use kinds
   use forsolver
   use functions_module
   implicit none
   real(rk) :: x0, tol, x_sol
   integer :: maxit

   ! Variable declaration
   x0    = 10.0_rk
   tol   = 1e-8_rk
   maxit = 100

   x_sol = solve(F, dFdx, x0, tol, maxit)

end program main

