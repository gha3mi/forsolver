program test3

   use kinds
   use functions_module
   use forsolver, only : nlsolver

   implicit none

   type(nlsolver) :: nls
   real(rk)       :: x_sol

   call nls%set_options(&
      nl_method   = 'newton',&
      maxit       = 100,&
      TolFun      = 1e-15_rk,&
      verbosity   = 1)

   call nls%solve(F=F1, dFdx=dF1dx, x0=10.0_rk, x_sol=x_sol)

end program test3