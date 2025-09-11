module my_function4
   use kinds, only: rk
   implicit none
contains
   function F1(x) result(F_val)
      real(rk), intent(in) :: x
      real(rk) :: F_val
      F_val = 5.0_rk * x**3 + 8.0_rk * x - 5.0_rk
   end function F1
   function dF1dx(x) result(dFdx_val)
      real(rk), intent(in) :: x
      real(rk) :: dFdx_val
      dFdx_val = 15.0_rk * x**2 + 8.0_rk
   end function dF1dx
end module my_function4

program test_solver4

   use kinds, only: rk
   use forsolver, only: nlsolver
   use my_function4, only: F1, dF1dx
   use forunittest, only: unit_test

   implicit none

   type(nlsolver) :: nls
   real(rk)       :: x, expected_x
   type(unit_test) :: ut

   call nls%set_options(&
      nl_method   = 'newton-modified',&
      nmp         = 2,&
      maxit       = 100,&
      TolFun      = 1e-4_rk,&
      verbosity   = 0)

   call nls%solve(F=F1, dFdx=dF1dx, x0=10.0_rk, x_sol=x)

   ! check if solution is close to ~=0.53128
   expected_x = 0.53128_rk
   call ut%check(x, expected_x, 1.0e-4_rk, 'test_solver4' )

end program test_solver4