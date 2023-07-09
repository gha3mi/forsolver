module functions_module
use kinds

implicit none

contains

function F1(x) result(F_val)
   real(rk), intent(in) :: x
   real(rk) :: F_val
   F_val = 5_rk * x**3 + 8_rk * x - 5_rk
end function F1

function dF1dx(x) result(dFdx_val)
   real(rk), intent(in) :: x
   real(rk) :: dFdx_val
   dFdx_val = 15_rk * x**2 + 8_rk
end function dF1dx

function F2(x) result(F_val)
   complex(rk), intent(in) :: x
   complex(rk) :: F_val
   F_val = 5_rk * x**3 + 8_rk * x - 5_rk
end function F2

function F3(x) result(F_val)
   real(rk), dimension(:), intent(in) :: x
   real(rk), dimension(:), allocatable :: F_val
   integer :: i

   allocate(F_val(3))
   F_val = 0.0_rk
   do i = 1,3
      F_val(i) = x(i)**3
   end do
end function F3

function dF3dx(x) result(dFdx_val)
   real(rk), dimension(:), intent(in) :: x
   real(rk), dimension(:,:), allocatable :: dFdx_val
   integer :: i, j
   real(rk), dimension(3,3) :: Idt

   Idt = 0.0_rk
   Idt(1,1) = 1.0_rk
   Idt(2,2) = 1.0_rk
   Idt(3,3) = 1.0_rk

   allocate(dFdx_val(3,3))
   dFdx_val = 0.0_rk
   do i = 1,3
      do j = 1,3
         dFdx_val(i,j) =  dFdx_val(i,j) + 3.0_rk*Idt(i,j)*x(i)**2
      end do
   end do

end function dF3dx

function F4(x) result(F_val)
   complex(rk), dimension(:), intent(in)  :: x
   complex(rk), dimension(:), allocatable :: F_val
   integer :: i

   allocate(F_val(3))
   F_val = 0.0_rk
   do i = 1,3
      F_val(i) = x(i)**3
   end do
end function F4

end module functions_module
