module forsolver

   ! This module provides functions and subroutines for
   ! solving linear systems and performing Newton's method.

   use :: kinds

   implicit none

   private

   public ::  solver

   !===============================================================================
   interface solver
      procedure :: solver_lin
      procedure :: newton_rel
   end interface
   !===============================================================================

contains

   !===============================================================================
   !> author: Seyed Ali Ghasemi
   pure function solver_lin(A, b, method) result(x)

      ! inputs
      real(rk),     dimension(:, :), contiguous, intent(in) :: A      ! input matrix A
      real(rk),     dimension(:),    contiguous, intent(in) :: b      ! right-hand side matrix b
      character(*), optional,                    intent(in) :: method

      ! outputs:
      real(rk), dimension(size(b))                       :: x    ! solution matrix x

      if (present(method)) then
         select case (method)
          case ('gesv')
            x = dgesv_rel(A, b)
          case ('gels')
            x = dgels_rel(A, b)
         end select
      else
         if (size(A,1)==size(A,2)) then
            x = dgesv_rel(A, b)
         else 
            x = dgels_rel(A, b)
         end if
      end if

   end function solver_lin
   !===============================================================================

   !===============================================================================
   !> author: Seyed Ali Ghasemi
   pure function dgesv_rel(A, b) result(x)
      ! inputs:
      real(rk), dimension(:, :), contiguous, intent(in)  :: A    ! input matrix A
      real(rk), dimension(:),    contiguous, intent(in)  :: b    ! right-hand side matrix b

      ! outputs:
      real(rk), dimension(size(b))                       :: x    ! solution matrix x

      ! local variables
      integer                                  :: info ! result info
      integer                                  :: n, lda, ldb
      integer,  dimension(size(A, 2))          :: ipiv
      real(rk), dimension(size(A,1),size(A,2)) :: a_copy
      real(rk), dimension(size(b))             :: b_copy

      ! interface for dgels subroutine
      interface
         pure subroutine dgesv(fn, fnrhs, fa, flda, fipiv, fb, fldb, finfo)
            use kinds
            integer,  intent(in)    :: fn, fnrhs, flda, fldb
            real(rk), intent(inout) :: fa(flda,*), fb(fldb,*)
            integer,  intent(out)   :: finfo
            integer,  intent(out)   :: fipiv(fn)
         end subroutine dgesv
      end interface

      b_copy = b
      a_copy = a

      ! get dimensions
      n    = size(A, 2)
      lda  = max(1, n)
      ldb  = max(1, n)

      ! call dgels subroutine
      call dgesv(n, 1, a_copy, n, ipiv, b_copy, n, info)

      ! copy the solution matrix
      x = b_copy

   end function dgesv_rel
   !===============================================================================

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
      integer                                  :: info ! result info
      integer                                  :: m, n, lda, ldb, lwork
      real(rk), allocatable                    :: work(:)
      real(rk)                                 :: work1(1)
      real(rk), dimension(size(A,1),size(A,2)) :: a_copy
      real(rk), dimension(size(b))             :: b_copy

      ! interface for dgels subroutine
      interface
         pure subroutine dgels(ftrans, fm, fn, fnrhs, fa, flda, fb, fldb, fwork, flwork, finfo)
            use kinds
            character(len=1), intent(in)    :: ftrans
            integer,          intent(in)    :: fm, fn, fnrhs, flda, fldb, flwork
            real(rk),         intent(inout) :: fa(flda,*), fb(fldb,*)
            real(rk),         intent(in)    :: fwork(*)
            integer,          intent(out)   :: finfo
         end subroutine dgels
      end interface

      a_copy = a
      b_copy = b

      ! get dimensions
      m    = size(A, 1)
      n    = size(A, 2)
      lda  = max(1, m)
      ldb  = max(1, max(m, n))

      ! calculate the optimal size of the work array
      call dgels('n', m, n, 1, a_copy, lda, b_copy, ldb, work1, -1, info)

      ! allocate work array
      lwork = nint(work1(1))
      allocate(work(lwork))

      ! call dgels subroutine
      call dgels('n', m, n, 1, a_copy, lda, b_copy, ldb, work, lwork, info)

      ! copy the solution matrix
      x = b_copy

      ! deallocate workspace
      deallocate(work)
   end function dgels_rel
   !===============================================================================

   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure function newton_rel(F, dFdx, x0, tol, maxit) result(x_sol)
      use kinds
      implicit none

      interface
         impure function Fun(x)
            use kinds
            real(rk), intent(in) :: x
            real(rk)             :: Fun
         end function Fun

         impure function dFun(x)
            use kinds
            real(rk), intent(in) :: x
            real(rk)             :: dFun
         end function dFun
      end interface

      procedure(Fun)  :: F
      procedure(dFun) :: dFdx

      real(rk), intent(in)  :: x0, tol
      integer,  intent(in)  :: maxit
      real(rk)              :: x_sol
      real(rk)              :: x, xnp
      real(rk)              :: F_val, dFdx_val, Krit
      integer               :: it
      logical               :: convergenz

      ! Variable declaration
      x          = x0
      xnp        = x0
      it         = 0
      convergenz = .false.

      write(*, '(a)') '-----------------------------------------------'
      write(*, '(a)') 'maxit             x0                   tol'
      write(*, '(i3, 10x, f12.8, 10x, e12.4)') maxit, x0, tol
      write(*, '(a)') '-----------------------------------------------'
      write(*, '(a)') 'start newton'
      write(*, '(a)') '-----------------------------------------------'
      write(*, '(a)') 'it        xn           F(xn)        dF(xn)/dxn'

      ! Main loop
      do while (.not. convergenz .and. it < maxit)
         F_val    = F(x)
         dFdx_val = dFdx(x)
         xnp      = x - F_val / dFdx_val

         Krit = abs(F_val)

         if (Krit <= tol) then
            convergenz = .true.
            x_sol      = x
         end if

         write(*, '(i3, f12.4, 4x, e12.4, 4x, e12.4)') it, x, F_val, dFdx_val
         it = it + 1

         x = xnp
      end do

      write(*, '(a)') '-----------------------------------------------'
      write(*, '(a)') 'end newton'
      write(*, '(a)') '-----------------------------------------------'
      write(*, '(a, g0)') 'x_sol = ', x_sol

   end function newton_rel
   !===============================================================================

end module forsolver
