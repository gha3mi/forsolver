module my_function12
   use forsolver_kinds, only: rk
   implicit none
contains
   function F3(x) result(F_val)
      real(rk), dimension(:), intent(in) :: x
      real(rk), dimension(:), allocatable :: F_val
      allocate(F_val(2))
      F_val(1) = 2.0_rk*x(1) - 400.0_rk*x(1) * (x(2) - x(1)**2) - 2.0_rk
      F_val(2) = 200.0_rk*x(2) - 200.0_rk*x(1)**2
   end function F3
end module my_function12

program test_solver12

   use forsolver_kinds, only: rk
   use forsolver, only: nlsolver
   use my_function12, only: F3
   use forunittest, only: unit_test

   implicit none

   type(nlsolver)         :: nls
   real(rk), dimension(2) :: x, expected_x
   type(unit_test) :: ut

   call nls%set_options(&
      nl_method   = 'newton-quasi-fd-modified',&
      fdm_method  = 'central',&
      fdm_tol     = 1e-7_rk,&
      nmp         = 1,&
      maxit       = 1000,&
      TolFun      = 1e-2_rk,&
      verbosity   = 0)

   call nls%solve(F=F3, x0=[0.95_rk,0.95_rk], x_sol=x)

   ! check if solution is close to [1,1]
   expected_x = [1.0_rk,1.0_rk]
   call ut%check(x, expected_x, 1.0e-1_rk, 'test_solver12' )

end program test_solver12
