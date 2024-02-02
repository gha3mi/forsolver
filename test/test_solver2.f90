program test_solver2

   use kinds
   use forsolver
   use forunittest

   implicit none

   real(rk), dimension(:,:), allocatable :: A
   real(rk), dimension(:)  , allocatable :: x1, x2, x3, expected_x, b
   integer                               :: m,n, i, j
   type(unit_test) :: ut

   m = 4
   n = 4

   allocate(A(m,n),b(m),x1(n),x2(n),x3(n),expected_x(n))

   A(1,:) = [ 1.0_rk,  1.0_rk,  -3.0_rk,  1.0_rk]
   A(2,:) = [-5.0_rk,  3.0_rk,  -4.0_rk,  1.0_rk]
   A(3,:) = [ 1.0_rk,  0.0_rk,   2.0_rk, -1.0_rk]
   A(4,:) = [ 1.0_rk,  2.0_rk,   0.0_rk,  0.0_rk]

   b = [2.0_rk, 0.0_rk, 1.0_rk, 12.0_rk]

   expected_x = [22.0_rk/17.0_rk, 91.0_rk/17.0_rk, 84.0_rk/17.0_rk, 173.0_rk/17.0_rk]

   x1 = solve(A, b)

   ! check if solution is close to expected_x
   call ut%check(x1, expected_x, 1.0e-6_rk, 'test_solver2.1' )

   x2 = solve(A, b, method='gesv')

   ! check if solution is close to expected_x
   call ut%check(x2, expected_x, 1.0e-6_rk, 'test_solver2.2' )

   x3 = solve(A, b, method='gels')

   ! check if solution is close to expected_x
   call ut%check(x3, expected_x, 1.0e-6_rk, 'test_solver2.3' )

   deallocate(A,b,x1,x2,x3,expected_x)

end program test_solver2
