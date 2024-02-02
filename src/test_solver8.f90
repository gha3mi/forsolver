module my_function8
   use kinds
   implicit none
contains
   function F2(x) result(F_val)
      complex(rk), intent(in) :: x
      complex(rk) :: F_val
      F_val = 5.0_rk * x**3 + 8.0_rk * x - 5.0_rk
   end function F2
end module my_function8

program test_solver8

   use forsolver
   use my_function8
   use forunittest

   implicit none

   type(nlsolver) :: nls
   complex(rk)    :: x
   real(rk)       :: expected_x
   type(unit_test) :: ut

   call nls%set_options(&
      nl_method   = 'newton-quasi-cs-modified',&
      cs_tol      = tiny(0.0_rk),&
      nmp         = 2,&
      maxit       = 100,&
      TolFun      = 1e-6_rk,&
      verbosity   = 0)

   call nls%solve(F=F2, x0=(0.95_rk, 0.0_rk), x_sol=x)

   ! check if solution is close to ~=0.53128
   expected_x = 0.53128_rk
   call ut%check(x%re, 0.53128_rk, 1.0e-4_rk, 'test_solver8' )

end program test_solver8
