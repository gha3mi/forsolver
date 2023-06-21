program test1

   use :: kinds
   use :: forsolver

   implicit none

   real(rk), dimension(:,:), allocatable :: A
   real(rk), dimension(:)  , allocatable :: x, b
   integer                               :: m,n, i, j

   m = 4
   n = 3

   allocate(A(m,n),b(m),x(n))

   call random_number(A)
   call random_number(b)
   A = A*10.0_rk
   b = b*10.0_rk

   X = solver(A, b)

   ! Print A
   print *, "A:"
   print "(4F10.6)", (A(:,j), j = 1, m)

   ! Print b
   print *, "b:"
   print "(4F10.6)", b


   ! Print x
   print *, "x:"
   print "(4F10.6)", x

   deallocate(A,b,x)

end program test1
