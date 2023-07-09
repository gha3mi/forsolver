program test6

   use kinds
   use functions_module
   use forsolver, only : nlsolver

   implicit none
   
   type(nlsolver) :: nls
   real(rk)       :: x_sol

   call nls%set_options(&
      nl_method   = 'newton-quasi-fd-modified',&
      fdm_method  = 'central',&
      fdm_tol     = 1e-8_rk,&
      nmp         = 2,&
      maxit       = 100,&
      TolFun      = 1e-15_rk,&
      verbosity   = 1)

   call nls%solve(F=F1, x0=10.0_rk, x_sol=x_sol)

end program test6