module forsolver

   !This module provides functions and subroutines for pseudoinverse calculations.

   use :: kinds

   implicit none

   private

   public ::  solver

   !===============================================================================
   interface solver
      procedure :: dgels_rel ! Interface for the pinverse_rel function
   end interface
   !===============================================================================

contains

   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> solves an overdetermined or underdetermined linear system using dgels.
   pure function dgels_rel(A, b) result(x)
      ! inputs:
      real(rk), dimension(:, :), contiguous, intent(in)  :: A    ! input matrix A
      real(rk), dimension(:),    contiguous, intent(in)  :: b    ! right-hand side matrix b

      ! outputs:
      real(rk), dimension(size(b))                        :: x    ! solution matrix x

      ! local variables
      integer                          :: info ! result info
      integer                          :: m, n, lda, ldb, lwork
      real(rk), allocatable            :: work(:)
      real(rk)                         :: work1(1)
      real(rk), dimension(size(b))     :: b_copy

      ! interface for dgels subroutine
      interface
         pure subroutine dgels(ftrans, fm, fn, fnrhs, fa, flda, fb, fldb, fwork, flwork, finfo)
            use kinds
            character(len=1), intent(in)    :: ftrans
            integer,          intent(in)    :: fm, fn, fnrhs, flda, fldb, flwork
            real(rk),         intent(inout) :: fa(flda,*), fb(fldb,*), fwork(*)
            integer,          intent(out)   :: finfo
         end subroutine dgels
      end interface

      b_copy = b

      ! get dimensions
      m    = size(A, 1)
      n    = size(A, 2)
      lda  = max(1, m)
      ldb  = max(1, max(m, n))

      ! calculate the optimal size of the work array
      call dgels('n', m, n, 1, A, lda, b_copy, ldb, work1, -1, info)

      ! allocate work array
      lwork = nint(work1(1))
      allocate(work(lwork))

      ! call dgels subroutine
      call dgels('n', m, n, 1, A, lda, b_copy, ldb, work, lwork, info)

      ! copy the solution matrix
      x = b_copy

      ! deallocate workspace
      deallocate(work)
   end function dgels_rel
   !===============================================================================

end module forsolver
