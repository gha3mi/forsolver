program test2

   use kinds
   use forsolver, only: solve

   implicit none

   real(rk), dimension(:,:), allocatable :: A
   real(rk), dimension(:)  , allocatable :: x, b
   integer                               :: m,n, i, j

   m = 4
   n = 4

   allocate(A(m,n),b(m),x(n))

   call random_number(A)
   call random_number(b)
   A = A*10.0_rk
   b = b*10.0_rk

   X = solve(A, b)
   X = solve(A, b, method='gesv')
   X = solve(A, b, method='gels')

   deallocate(A,b,x)

end program test2