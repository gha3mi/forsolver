program test_solver16

   use kinds, only: rk
   use forsolver, only: solve
   use forunittest, only: unit_test

   implicit none

   real(rk), dimension(4,4) :: A
   real(rk), dimension(4)   :: b
   real(rk), dimension(4)   :: x, expected_x
   type(unit_test) :: ut

   ! set the matrix A
   A(1,:) = [ 1.44_rk, -7.84_rk, -4.39_rk,  4.53_rk]
   A(2,:) = [-9.96_rk, -0.28_rk, -3.24_rk,  3.83_rk]
   A(3,:) = [-7.55_rk,  3.24_rk,  6.27_rk, -6.64_rk]
   A(4,:) = [ 8.34_rk,  8.09_rk,  5.28_rk,  2.06_rk]

   ! set the right-hand side
   b = [8.58_rk, 8.26_rk, 8.48_rk,-5.28_rk]

   ! solve the system
   x = solve(A, b) ! X = solve(A, b, method='gesvs')

   ! expected result
   expected_x(1) = -1.0544691129297037_rk
   expected_x(2) = -1.9149827187319857_rk
   expected_x(3) =  2.9192679369935912_rk
   expected_x(4) =  1.7440523733249165_rk

   ! check the result
   call ut%check(x, expected_x, 1.0e-5_rk, 'test_solver16' )

end program test_solver16