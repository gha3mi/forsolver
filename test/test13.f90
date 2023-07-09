program test13

   use kinds
   use functions_module
   use forsolver, only : nlsolver

   implicit none

   type(nlsolver)            :: nls
   complex(rk), dimension(2) :: x_sol

   call nls%set_options(&
      nl_method   = 'newton-quasi-cs',&
      cs_tol      = 1e-200_rk,&
      maxit       = 10800,&
      TolFun      = 1e-13_rk,&
      verbosity   = 1)

   call nls%solve(F=F4, x0=[(-1.0_rk,0.0_rk) ,(-1.0_rk,0.0_rk)], x_sol=x_sol)

end program test13
