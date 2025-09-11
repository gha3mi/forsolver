module my_function9
   use kinds, only: rk
   implicit none
contains
   function F3(x) result(F_val)
      real(rk), dimension(:), intent(in) :: x
      real(rk), dimension(:), allocatable :: F_val
      allocate(F_val(2))
      F_val(1) = 2.0_rk*x(1) - 400.0_rk*x(1) * (x(2) - x(1)**2) - 2.0_rk
      F_val(2) = 200.0_rk*x(2) - 200.0_rk*x(1)**2
   end function F3
   function dF3dx(x) result(dFdx_val)
      real(rk), dimension(:), intent(in) :: x
      real(rk), dimension(:,:), allocatable :: dFdx_val
      allocate(dFdx_val(2,2))
      dFdx_val(1,1) = 1200.0_rk*x(1)**2 - 400.0_rk*x(2) + 2.0_rk
      dFdx_val(1,2) = - 400.0_rk*x(1)
      dFdx_val(2,1) = - 400.0_rk*x(1)
      dFdx_val(2,2) = 200.0_rk
   end function dF3dx
end module my_function9

program test_solver9

   use kinds, only: rk
   use forsolver, only: nlsolver
   use my_function9, only: F3, dF3dx
   use forunittest, only: unit_test

   implicit none

   type(nlsolver)         :: nls
   real(rk), dimension(2) :: x, expected_x
   type(unit_test) :: ut

   call nls%set_options(&
      nl_method   = 'newton',&
      maxit       = 100,&
      TolFun      = 1e-15_rk,&
      verbosity   = 0)

   call nls%solve(F=F3, dFdx=dF3dx, x0=[0.95_rk,0.95_rk], x_sol=x)

   ! check if solution is close to [1,1]
   expected_x = [1.0_rk,1.0_rk]
   call ut%check(x, expected_x, 1.0e-5_rk, 'test_solver9' )

end program test_solver9
