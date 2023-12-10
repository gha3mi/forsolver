program test10

   use kinds
   use functions_module
   use forsolver, only: nlsolver

   implicit none
   
   type(nlsolver)         :: nls
   real(rk), dimension(2) :: x_sol

   call nls%set_options(&
      nl_method   = 'newton-modified',&
      nmp         = 1,&
      maxit       = 100,&
      TolFun      = 1e-15_rk,&
      verbosity   = 1)

      call nls%solve(F=F3, dFdx=dF3dx, x0=[-1.0_rk,-1.0_rk], x_sol=x_sol)

end program test10