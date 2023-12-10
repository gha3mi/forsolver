program test7

   use kinds
   use functions_module
   use forsolver, only: nlsolver

   implicit none

   type(nlsolver) :: nls
   complex(rk)    :: x_sol

   call nls%set_options(&
      nl_method   = 'newton-quasi-cs',&
      cs_tol      = 1e-200_rk,&
      maxit       = 100,&
      TolFun      = 1e-15_rk,&
      verbosity   = 1)

   call nls%solve(F=F2, x0=(1.0_rk,0.0_rk), x_sol=x_sol)

end program test7
