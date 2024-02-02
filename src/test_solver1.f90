program test_solver1

   use kinds
   use forsolver
   use forunittest

   implicit none

   real(rk), dimension(:,:), allocatable :: A
   real(rk), dimension(:)  , allocatable :: x1, x2, expected_x, b
   integer                               :: m,n, i, j
   type(unit_test) :: ut

   m = 3
   n = 2

   allocate(A(m,n),b(m),x1(n),x2(n),expected_x(n))

   A(1,:) = [ 1.0_rk, 5.0_rk]
   A(2,:) = [ 3.0_rk, 1.0_rk]
   A(3,:) = [-2.0_rk, 4.0_rk]

   b = [4.0_rk, -2.0_rk, 3.0_rk]

   expected_x = [-4.0_rk/7.0_rk, 5.0_rk/7.0_rk]

   x1 = solve(A, b)

   ! check if solution is close to expected_x
   call ut%check(x1, expected_x, 1.0e-6_rk, 'test_solver1.1' )

   x2 = solve(A, b, method='gels')

   ! check if solution is close to expected_x
   call ut%check(x2, expected_x, 1.0e-6_rk, 'test_solver1.2' )

   deallocate(A,b,x1,x2,expected_x)

end program test_solver1