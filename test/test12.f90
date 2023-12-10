program test12

   use kinds
   use functions_module
   use forsolver, only: nlsolver

   implicit none
   
   type(nlsolver)         :: nls
   real(rk), dimension(2) :: x_sol

   call nls%set_options(&
      nl_method   = 'newton-quasi-fd-modified',&
      fdm_method  = 'central',&
      fdm_tol     = 1e-6_rk,&
      nmp         = 1,&
      maxit       = 1000000,&
      TolFun      = 1e-13_rk,&
      verbosity   = 1)

      call nls%solve(F=F3, x0=[-1.0_rk,-1.0_rk], x_sol=x_sol)

end program test12