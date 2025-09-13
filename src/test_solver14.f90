module my_function14
   use forsolver_kinds, only: rk
   implicit none
contains
   function F4(x) result(F_val)
      complex(rk), dimension(:), intent(in)  :: x
      complex(rk), dimension(:), allocatable :: F_val
      allocate(F_val(2))
      F_val(1) = 2.0_rk*x(1) - 400.0_rk*x(1) * (x(2) - x(1)**2) - 2.0_rk
      F_val(2) = 200.0_rk*x(2) - 200.0_rk*x(1)**2
   end function F4
end module my_function14

program test_solver14

   use forsolver_kinds, only: rk
   use forsolver, only: nlsolver
   use my_function14, only: F4
   use forunittest, only: unit_test

   implicit none

   type(nlsolver)            :: nls
   complex(rk), dimension(2) :: x
   real(rk),    dimension(2) :: expected_x
   type(unit_test) :: ut

   call nls%set_options(&
      nl_method   = 'newton-quasi-cs-modified',&
      cs_tol      = tiny(0.0_rk),&
      nmp         = 1,&
      maxit       = 1000,&
      TolFun      = 1e-2_rk,&
      verbosity   = 0)

   call nls%solve(F=F4, x0=[(0.95_rk,0.0_rk) ,(0.95_rk,0.0_rk)], x_sol=x)

   ! check if solution is close to [1,1]
   expected_x = [1.0_rk,1.0_rk]
   call ut%check(real(x,rk), expected_x, 1.0e-1_rk, 'test_solver14' )
   ! call ut%check(x%re, expected_x, 1.0e-1_rk, 'test_solver14' )

end program test_solver14
