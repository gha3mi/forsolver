program test16

   use kinds
   use forsolver, only: solve

   implicit none

   real(rk), dimension(4,4) :: A
   real(rk), dimension(4)   :: b
   real(rk), dimension(4)   :: x, x_expected

   ! set the matrix A
   A(1,:) = [ 1.44_rk, -7.84_rk, -4.39_rk,  4.53_rk]
   A(2,:) = [-9.96_rk, -0.28_rk, -3.24_rk,  3.83_rk]
   A(3,:) = [-7.55_rk,  3.24_rk,  6.27_rk, -6.64_rk]
   A(4,:) = [ 8.34_rk,  8.09_rk,  5.28_rk,  2.06_rk]

   ! set the right-hand side
   b = [8.58_rk, 8.26_rk, 8.48_rk,-5.28_rk]

   ! solve the system
   X = solve(A, b) ! X = solve(A, b, method='gesvs')

   ! expected result
   x_expected(1) = -1.0544691129297037_rk
   x_expected(2) = -1.9149827187319857_rk
   x_expected(3) =  2.9192679369935912_rk
   x_expected(4) =  1.7440523733249165_rk

   ! check the result
   if (norm2(X - x_expected) < 1e-6_rk) then
      print'(a)', 'test16: passed'
   else
      print'(a)', 'test16: failed'
   end if

end program test16