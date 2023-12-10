program test15

   use kinds
   use forsolver, only: solve

   implicit none

   real(rk), dimension(6,4) :: A
   real(rk), dimension(6)   :: b
   real(rk), dimension(4)   :: x, x_expected

   ! set the matrix A
   A(1,:) = [ 1.44_rk, -7.84_rk, -4.39_rk,  4.53_rk]
   A(2,:) = [-9.96_rk, -0.28_rk, -3.24_rk,  3.83_rk]
   A(3,:) = [-7.55_rk,  3.24_rk,  6.27_rk, -6.64_rk]
   A(4,:) = [ 8.34_rk,  8.09_rk,  5.28_rk,  2.06_rk]
   A(5,:) = [ 7.08_rk,  2.52_rk,  0.74_rk, -2.47_rk]
   A(6,:) = [-5.45_rk, -5.70_rk, -1.19_rk,  4.70_rk]

   ! set the right-hand side
   b = [8.58_rk, 8.26_rk, 8.48_rk,-5.28_rk, 5.72_rk, 8.93_rk]

   ! solve the system
   X = solve(A, b) ! X = solve(A, b, method='gels')

   ! expected result
   x_expected(1) = -0.45063713541953410_rk
   x_expected(2) = -0.84915021471399577_rk
   x_expected(3) =  0.70661216240939595_rk
   x_expected(4) =  0.12888575215577794_rk

   ! check the result
   if (norm2(X - x_expected) < 1e-6_rk) then
      print'(a)', 'test15: passed'
   else
      print'(a)', 'test15: failed'
   end if

end program test15