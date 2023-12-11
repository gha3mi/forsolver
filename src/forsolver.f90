module forsolver

   ! This module provides functions and subroutines for
   ! solving linear systems and performing Newton's method.

   use kinds
   use fordiff

   implicit none

   private

   public ::  solve, nlsolver

   type :: nlsolver
      character(:), allocatable :: lin_method
      character(:), allocatable :: nl_method
      character(:), allocatable :: fdm_method
      real(rk)                  :: TolFun
      real(rk)                  :: fdm_tol
      real(rk)                  :: cs_tol
      integer                   :: maxit
      integer                   :: nmp
      integer                   :: verbosity
      ! character(:), allocatable :: stepsize
      ! real(rk)                  :: alpha0
      ! real(rk)                  :: c1
      ! real(rk)                  :: c2
   contains
      procedure :: set_options
      procedure :: newton_rel_T0
      procedure :: newton_rel_T1
      procedure :: newton_complex_step_rel_T0
      procedure :: newton_complex_step_rel_T1
      generic   :: solve => newton_rel_T0,&
         newton_rel_T1,&
         newton_complex_step_rel_T0,&
         newton_complex_step_rel_T1
      final     :: deallocate_solver
   end type nlsolver

   !===============================================================================
   interface solve
      procedure :: solver_lin
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
      real(rk), dimension(max(1, size(A, 2))) :: x      ! solution matrix x
      
      ! local variables
      integer :: info ! result info

      ! call solver
      if (present(method)) then
         select case (method)
          case ('gesv')
            call dgesv_rel(A, b, x, info)
          case ('gels')
            call dgels_rel(A, b, x, info)
         end select
      else
         if (size(A,1)==size(A,2)) then
            call dgesv_rel(A, b, x, info)
         else
            call dgels_rel(A, b, x, info)
         end if
      end if

   end function solver_lin
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   pure subroutine dgesv_rel(A, b, x, info)
      ! inputs:
      real(rk), dimension(:, :), contiguous, intent(in) :: A    ! input matrix A
      real(rk), dimension(:),    contiguous, intent(in) :: b    ! right-hand side matrix b

      ! outputs:
      real(rk), dimension(max(1, size(A, 2))), intent(out) :: x    ! solution matrix x
      integer,                                 intent(out) :: info ! result info

      ! local variables
      integer                               :: n, lda, ldb, nrhs
      integer,  dimension(size(A, 2))       :: ipiv
      real(rk), dimension(:,:), allocatable :: a_copy
      real(rk), dimension(:,:), allocatable :: b_copy

      ! interface for dgels subroutine
      interface
         pure subroutine dgesv(fn, fnrhs, fa, flda, fipiv, fb, fldb, finfo)
            import rk
            integer,  intent(in)    :: fn, fnrhs, flda, fldb
            real(rk), intent(inout) :: fa(flda,fn), fb(fldb,fnrhs)
            integer,  intent(out)   :: finfo
            integer,  intent(out)   :: fipiv(fn)
         end subroutine dgesv
      end interface

      ! get dimensions
      nrhs = 1 ! size(b, 2)
      n    = size(A, 2)
      lda  = max(1, n)
      ldb  = max(1, n)

      ! copy the input matrices
      a_copy = a
      allocate(b_copy(ldb, nrhs))
      b_copy(:, 1) = b

      ! call dgels subroutine
      call dgesv(n, nrhs, a_copy, lda, ipiv, b_copy, ldb, info)

      ! copy the solution matrix
      if (info == 0) then
         x = b_copy(1:ldb, 1) ! nrhs = 1
      else
         error stop 'dgesv failed'
      end if

   end subroutine dgesv_rel
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> solves an overdetermined or underdetermined linear system using dgels.
   pure subroutine dgels_rel(A, b, x, info)
      ! inputs:
      real(rk), dimension(:, :), contiguous, intent(in) :: A    ! input matrix A
      real(rk), dimension(:),    contiguous, intent(in) :: b    ! right-hand side matrix b

      ! outputs:
      real(rk), dimension(max(1, size(A, 2))), intent(out) :: x    ! solution matrix x
      integer,                                 intent(out) :: info ! result info

      ! local variables
      character(1)                          :: trans
      integer                               :: m, n, lda, ldb, lwork, nrhs
      real(rk), allocatable                 :: work(:)
      real(rk)                              :: work1(1)
      real(rk), dimension(:,:), allocatable :: a_copy
      real(rk), dimension(:,:), allocatable :: b_copy

      ! interface for dgels subroutine
      interface
         pure subroutine dgels(ftrans, fm, fn, fnrhs, fa, flda, fb, fldb, fwork, flwork, finfo)
            import rk
            character(len=1), intent(in)    :: ftrans
            integer,          intent(in)    :: fm, fn, fnrhs, flda, fldb, flwork
            real(rk),         intent(inout) :: fa(flda,*), fb(fldb,*)
            real(rk),         intent(in)    :: fwork(*)
            integer,          intent(out)   :: finfo
         end subroutine dgels
      end interface

      ! 
      trans = 'n'

      ! get dimensions
      nrhs = 1 ! size(b, 2)
      m    = size(A, 1)
      n    = size(A, 2)
      lda  = max(1, m)
      ldb  = max(1, max(m, n))

      ! copy the input matrices
      a_copy = a
      allocate(b_copy(ldb, nrhs))
      b_copy(:, 1) = b

      ! calculate the optimal size of the work array
      call dgels(trans, m, n, nrhs, a_copy, lda, b_copy, ldb, work1, -1, info)

      ! allocate work array
      lwork = nint(work1(1))
      allocate(work(lwork))

      ! call dgels subroutine
      call dgels(trans, m, n, nrhs, a_copy, lda, b_copy, ldb, work, lwork, info)

      ! copy the solution matrix
      if (info == 0) then
         if (trans == 'n') x = b_copy(1:n, 1) ! nrhs = 1
         if (trans == 't') x = b_copy(1:m, 1) ! nrhs = 1
      else
         error stop 'dgels failed'
      end if

      ! deallocate workspace
      deallocate(work)
   end subroutine dgels_rel
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine newton_rel_T0(this, F, dFdx, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            real(rk), intent(in) :: x
            real(rk)             :: Fun
         end function Fun

         impure function dFun(x)
            import rk
            real(rk), intent(in) :: x
            real(rk)             :: dFun
         end function dFun
      end interface

      procedure(Fun)            :: F
      procedure(dFun), optional :: dFdx

      class(nlsolver), intent(inout) :: this
      real(rk),        intent(in)    :: x0
      real(rk),        intent(out)   :: x_sol

      if (this%verbosity == 1) then
         print '(a)', '-----------------------------------------------'
         print '(a)', 'maxit             x0                   tol'
         print '(g0, 10x, f12.8, 10x, e12.4)', this%maxit, x0, this%TolFun
         print '(a)', '-----------------------------------------------'
         print '(a)', 'start newton'
         print '(a)', '-----------------------------------------------'
         print '(a)', 'it        xn           F(xn)         dF(xn)/dxn'
      end if

      select case (this%nl_method)
       case ('newton')
         call newton_method_T0(this, F, dFdx, x0,  x_sol)
       case ('newton-modified')
         call modified_newton_method_T0(this, F, dFdx, x0,  x_sol)
       case ('newton-quasi-fd')
         call quasi_fd_newton_method_T0(this, F, x0,  x_sol)
       case ('newton-quasi-fd-modified')
         call modified_quasi_fd_newton_method_T0(this, F, x0,  x_sol)
         !  case ('newton-quasi-bfgs')
         !    call quasi_bfgs_newton_method_T0(this, F, x0,  x_sol)
         !  case ('newton-quasi-bfgs-modified')
         !    call modified_quasi_bfgs_newton_method_T0(this, F, x0,  x_sol)
      end select

      if (this%verbosity == 1) then
         print '(a)', '-----------------------------------------------'
         print '(a)', 'end newton'
         print '(a)', '-----------------------------------------------'
         print '(a, g0)', 'x_sol = ', x_sol
      end if


   end subroutine newton_rel_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine newton_rel_T1(this, F, dFdx, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            real(rk), dimension(:), intent(in)  :: x
            real(rk), dimension(:), allocatable :: Fun
         end function Fun

         impure function dFun(x)
            import rk
            real(rk), dimension(:), intent(in)    :: x
            real(rk), dimension(:,:), allocatable :: dFun
         end function dFun
      end interface

      procedure(Fun)            :: F
      procedure(dFun), optional :: dFdx

      class(nlsolver),               intent(inout) :: this
      real(rk), dimension(:),        intent(in)    :: x0
      real(rk), dimension(size(x0)), intent(out)   :: x_sol
      integer                                      :: i

      if (this%verbosity == 1) then
         print '(a)', '-----------------------------------------------'
         print '(a)', 'maxit             tol'
         print '(g0, 10x, e12.4)', this%maxit, this%TolFun
         print '(a)', '-----------------------------------------------'
         print '(a)', 'start newton'
         print '(a)', '-----------------------------------------------'
         print '(a)', 'it     ||F||'
      end if

      select case (this%nl_method)
       case ('newton')
         call newton_method_T1(this, F, dFdx, x0,  x_sol)
       case ('newton-modified')
         call modified_newton_method_T1(this, F, dFdx, x0,  x_sol)
       case ('newton-quasi-fd')
         call quasi_fd_newton_method_T1(this, F, x0,  x_sol)
       case ('newton-quasi-fd-modified')
         call modified_quasi_fd_newton_method_T1(this, F, x0,  x_sol)
         ! case ('newton-quasi-bfgs')
         !    call quasi_bfgs_newton_method_T1(this, F, x0,  x_sol)
         ! case ('newton-quasi-bfgs-modified')
         !    call modified_quasi_bfgs_newton_method_T1(this, F, x0,  x_sol)
      end select

      if (this%verbosity == 1) then
         print '(a)', '-----------------------------------------------'
         print '(a)', 'end newton'
         print '(a)', '-----------------------------------------------'
         do i = 1,size(x_sol)
            print '(a, g0)', 'x_sol = ', x_sol(i)
         end do
      end if

   end subroutine newton_rel_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine newton_complex_step_rel_T0(this, F, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            complex(rk), intent(in) :: x
            complex(rk)             :: Fun
         end function Fun
      end interface

      procedure(Fun) :: F

      class(nlsolver), intent(inout) :: this
      complex(rk),     intent(in)    :: x0
      complex(rk),     intent(out)   :: x_sol

      if (this%verbosity == 1) then
         print '(a)', '-----------------------------------------------'
         print '(a)', 'maxit             x0                   tol'
         print '(g0, 10x, f12.8, 10x, e12.4)', this%maxit, real(x0, kind=rk), this%TolFun
         print '(a)', '-----------------------------------------------'
         print '(a)', 'start newton'
         print '(a)', '-----------------------------------------------'
         print '(a)', 'it        xn           F(xn)        dF(xn)/dxn'
      end if

      select case (this%nl_method)
       case ('newton-quasi-cs')
         call quasi_cs_newton_method_T0(this, F, x0,  x_sol)
       case ('newton-quasi-cs-modified')
         call modified_quasi_cs_newton_method_T0(this, F, x0,  x_sol)
      end select

      if (this%verbosity == 1) then
         print '(a)', '-----------------------------------------------'
         print '(a)', 'end newton'
         print '(a)', '-----------------------------------------------'
         print '(a, g0)', 'x_sol = ', x_sol
      end if

   end subroutine newton_complex_step_rel_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine newton_complex_step_rel_T1(this, F, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            complex(rk), dimension(:), intent(in)  :: x
            complex(rk), dimension(:), allocatable :: Fun
         end function Fun
      end interface

      procedure(Fun)  :: F

      class(nlsolver),                  intent(inout) :: this
      complex(rk), dimension(:),        intent(in)    :: x0
      complex(rk), dimension(size(x0)), intent(out)   :: x_sol
      integer                                         :: i

      if (this%verbosity == 1) then
         print '(a)', '-----------------------------------------------'
         print '(a)', 'maxit             tol'
         print '(g0, 10x, f12.8, e12.4)', this%maxit, this%TolFun
         print '(a)', '-----------------------------------------------'
         print '(a)', 'start newton'
         print '(a)', '-----------------------------------------------'
         print '(a)', 'it     ||F||'
      end if

      select case (this%nl_method)
       case ('newton-quasi-cs')
         call quasi_cs_newton_method_T1(this, F, x0,  x_sol)
       case ('newton-quasi-cs-modified')
         call modified_quasi_cs_newton_method_T1(this, F, x0,  x_sol)
      end select

      if (this%verbosity == 1) then
         print '(a)', '-----------------------------------------------'
         print '(a)', 'end newton'
         print '(a)', '-----------------------------------------------'
         do i = 1,size(x_sol)
            print '(a, g0)', 'x_sol = ', real(x_sol(i), kind=rk)
         end do
         ! print '(a, g0)', 'x_sol = ', x_sol
      end if

   end subroutine newton_complex_step_rel_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine set_options(this,&
      nl_method, lin_method, maxit, TolFun, alpha0, c1, c2, nmp, fdm_method, fdm_tol, cs_tol, stepsize, verbosity)
      class(nlsolver), intent(inout)        :: this
      character(*),    intent(in), optional :: nl_method
      character(*),    intent(in), optional :: lin_method
      character(*),    intent(in), optional :: stepsize
      character(*),    intent(in), optional :: fdm_method
      real(rk),        intent(in), optional :: TolFun
      real(rk),        intent(in), optional :: fdm_tol
      real(rk),        intent(in), optional :: cs_tol
      integer,         intent(in), optional :: maxit
      real(rk),        intent(in), optional :: alpha0
      real(rk),        intent(in), optional :: c1
      real(rk),        intent(in), optional :: c2
      integer,         intent(in), optional :: nmp
      integer,         intent(in), optional :: verbosity

      if (present(nl_method)) then
         this%nl_method = nl_method
      else
         this%nl_method = 'newton'
      end if

      if (present(lin_method)) then
         this%lin_method  = lin_method
      else
         ! this%lin_method  =
      end if

      if (present(fdm_method)) then
         this%fdm_method  = fdm_method
      else
         this%fdm_method  = 'forward'
      end if

      if (present(maxit)) then
         this%maxit  = maxit
      else
         this%maxit  = 100
      end if

      if (present(TolFun)) then
         this%TolFun  = TolFun
      else
         this%TolFun  = 1e-4_rk
      end if

      if (present(fdm_tol)) then
         this%fdm_tol  = fdm_tol
      else
         this%fdm_tol  = 1e-4_rk
      end if

      if (present(cs_tol)) then
         this%cs_tol  = cs_tol
      else
         this%cs_tol  = 1e-100_rk
      end if

      if (present(nmp)) then
         this%nmp  = nmp
      else
         this%nmp  = 2
      end if

      if (present(verbosity)) then
         this%verbosity  = verbosity
      else
         this%verbosity  = 1
      end if

      ! if (present(stepsize))   this%stepsize   = stepsize
      ! if (present(alpha0))     this%alpha0     = alpha0
      ! if (present(c1))         this%c1         = c1
      ! if (present(c2))         this%c2         = c2

   end subroutine set_options
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   elemental pure subroutine deallocate_solver(this)
      type(nlsolver), intent(inout)     :: this

      if (allocated(this%nl_method)) deallocate(this%nl_method)
      if (allocated(this%fdm_method)) deallocate(this%fdm_method)
   end subroutine deallocate_solver
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine newton_method_T0(this, F, dFdx, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            real(rk), intent(in)  :: x
            real(rk) :: Fun
         end function Fun

         impure function dFun(x)
            import rk
            real(rk), intent(in)    :: x
            real(rk) :: dFun
         end function dFun
      end interface

      procedure(Fun)  :: F
      procedure(dFun) :: dFdx

      class(nlsolver), intent(inout) :: this
      real(rk),        intent(in)    :: x0
      real(rk),        intent(out)   :: x_sol
      real(rk)                       :: xk
      real(rk)                       :: F_val
      real(rk)                       :: dFdx_val
      real(rk)                       :: criteriaFun
      integer                        :: k
      logical                        :: convergenz
      real(rk)                       :: dk
      real(rk)                       :: alphak

      k          = 0
      xk         = x0
      convergenz = .false.

      do while (.not. convergenz .and. k < this%maxit)
         F_val = F(xk)
         dFdx_val = dFdx(xk)

         criteriaFun = abs(F_val)

         if (this%verbosity == 1) then
            print '(g0, f12.4, 4x, e12.4, 4x, e12.4)', k, xk, F_val, dFdx_val
         end if

         if (criteriaFun <= this%TolFun) then
            convergenz = .true.
            x_sol      = xk
            return
         else
            dk = - F_val / dFdx_val
            alphak = 1.0_rk
            xk = xk + alphak*dk
            k = k + 1
         end if
      end do
   end subroutine newton_method_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine modified_newton_method_T0(this, F, dFdx, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            real(rk), intent(in)  :: x
            real(rk) :: Fun
         end function Fun

         impure function dFun(x)
            import rk
            real(rk), intent(in)    :: x
            real(rk) :: dFun
         end function dFun
      end interface

      procedure(Fun)  :: F
      procedure(dFun) :: dFdx

      class(nlsolver), intent(inout) :: this
      real(rk),        intent(in)    :: x0
      real(rk),        intent(out)   :: x_sol
      real(rk)                       :: xk
      real(rk)                       :: F_val
      real(rk)                       :: dFdx_val
      real(rk)                       :: criteriaFun
      integer                        :: k
      logical                        :: convergenz
      real(rk)                       :: dk
      real(rk)                       :: alphak

      k          = 0
      xk         = x0
      convergenz = .false.

      do while (.not. convergenz .and. k < this%maxit)
         F_val = F(xk)
         if ((mod(k, this%nmp) == 0)) dFdx_val = dFdx(xk)

         criteriaFun = abs(F_val)

         if (this%verbosity == 1) then
            print '(g0, f12.4, 4x, e12.4, 4x, e12.4)', k, xk, F_val, dFdx_val
         end if

         if (criteriaFun <= this%TolFun) then
            convergenz = .true.
            x_sol      = xk
            return
         else
            dk = - F_val / dFdx_val
            alphak = 1.0_rk
            xk = xk + alphak*dk
            k = k + 1
         end if
      end do
   end subroutine modified_newton_method_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine quasi_fd_newton_method_T0(this, F, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            real(rk), intent(in) :: x
            real(rk)             :: Fun
         end function Fun
      end interface

      procedure(Fun) :: F

      class(nlsolver), intent(inout) :: this
      real(rk),        intent(in)    :: x0
      real(rk),        intent(out)   :: x_sol
      real(rk)                       :: xk
      real(rk)                       :: F_val
      real(rk)                       :: dFdx_val
      real(rk)                       :: criteriaFun
      integer                        :: k
      logical                        :: convergenz
      real(rk)                       :: pk
      real(rk)                       :: qk
      real(rk)                       :: dk
      real(rk)                       :: alphak

      k          = 0
      xk         = x0
      convergenz = .false.

      do while (.not. convergenz .and. k < this%maxit)
         F_val = F(xk)
         dFdx_val = derivative(f=F, x=xk, h=this%fdm_tol, method=this%fdm_method)

         criteriaFun = abs(F_val)

         if (this%verbosity == 1) then
            print '(g0, f12.4, 4x, e12.4, 4x, e12.4)', k, xk, F_val, dFdx_val
         end if

         if (criteriaFun <= this%TolFun) then
            convergenz = .true.
            x_sol      = xk
            return
         else
            dk = - F_val / dFdx_val
            alphak = 1.0_rk
            xk = xk + alphak*dk
            k = k + 1
         end if
      end do
   end subroutine quasi_fd_newton_method_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine modified_quasi_fd_newton_method_T0(this, F, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            real(rk), intent(in) :: x
            real(rk)             :: Fun
         end function Fun
      end interface

      procedure(Fun) :: F

      class(nlsolver), intent(inout) :: this
      real(rk),        intent(in)    :: x0
      real(rk),        intent(out)   :: x_sol
      real(rk)                       :: xk
      real(rk)                       :: F_val
      real(rk)                       :: dFdx_val
      real(rk)                       :: criteriaFun
      integer                        :: k
      logical                        :: convergenz
      real(rk)                       :: pk
      real(rk)                       :: qk
      real(rk)                       :: dk
      real(rk)                       :: alphak

      k          = 0
      xk         = x0
      convergenz = .false.

      do while (.not. convergenz .and. k < this%maxit)
         F_val = F(xk)
         if (mod(k, this%nmp) == 0) dFdx_val = derivative(f=F, x=xk, h=this%fdm_tol, method=this%fdm_method)

         criteriaFun = abs(F_val)

         if (this%verbosity == 1) then
            print '(g0, f12.4, 4x, e12.4, 4x, e12.4)', k, xk, F_val, dFdx_val
         end if

         if (criteriaFun <= this%TolFun) then
            convergenz = .true.
            x_sol      = xk
            return
         else
            dk = - F_val / dFdx_val
            alphak = 1.0_rk
            xk = xk + alphak*dk
            k = k + 1
         end if
      end do
   end subroutine modified_quasi_fd_newton_method_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine newton_method_T1(this, F, dFdx, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            real(rk), dimension(:), intent(in)  :: x
            real(rk), dimension(:), allocatable :: Fun
         end function Fun

         impure function dFun(x)
            import rk
            real(rk), dimension(:), intent(in)    :: x
            real(rk), dimension(:,:), allocatable :: dFun
         end function dFun
      end interface

      procedure(Fun)  :: F
      procedure(dFun) :: dFdx

      class(nlsolver),               intent(inout) :: this
      real(rk), dimension(:),        intent(in)    :: x0
      real(rk), dimension(size(x0)), intent(out)   :: x_sol
      real(rk), dimension(size(x0))                :: xk
      real(rk), dimension(:),   allocatable        :: F_val
      real(rk), dimension(:,:), allocatable        :: dFdx_val
      real(rk)                                     :: criteriaFun
      integer                                      :: k
      logical                                      :: convergenz
      real(rk), dimension(size(x0))                :: dk
      real(rk)                                     :: alphak

      k          = 0
      xk         = x0
      convergenz = .false.

      do while (.not. convergenz .and. k < this%maxit)
         F_val    = F(xk)
         dFdx_val = dFdx(xk)

         criteriaFun = norm2(F_val)

         if (this%verbosity == 1) then
            print '(g0, e12.4)', k, criteriaFun
         end if

         if (criteriaFun <= this%TolFun) then
            convergenz = .true.
            x_sol      = xk
            return
         else
            dk = - solve(dFdx_val, F_val, this%lin_method)
            alphak = 1.0_rk
            xk = xk + alphak*dk
            k = k + 1
         end if
      end do
   end subroutine newton_method_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine modified_newton_method_T1(this, F, dFdx, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            real(rk), dimension(:), intent(in)  :: x
            real(rk), dimension(:), allocatable :: Fun
         end function Fun

         impure function dFun(x)
            import rk
            real(rk), dimension(:), intent(in)    :: x
            real(rk), dimension(:,:), allocatable :: dFun
         end function dFun
      end interface

      procedure(Fun)  :: F
      procedure(dFun) :: dFdx

      class(nlsolver),               intent(inout) :: this
      real(rk), dimension(:),        intent(in)    :: x0
      real(rk), dimension(size(x0)), intent(out)   :: x_sol
      real(rk), dimension(size(x0))                :: xk
      real(rk), dimension(:),   allocatable        :: F_val
      real(rk), dimension(:,:), allocatable        :: dFdx_val
      real(rk)                                     :: criteriaFun
      integer                                      :: k
      logical                                      :: convergenz
      real(rk), dimension(size(x0))                :: dk
      real(rk)                                     :: alphak

      k          = 0
      xk         = x0
      convergenz = .false.

      do while (.not. convergenz .and. k < this%maxit)
         F_val    = F(xk)
         if (mod(k, this%nmp) == 0) dFdx_val = dFdx(xk)

         criteriaFun = norm2(F_val)

         if (this%verbosity == 1) then
            print '(g0, e12.4)', k, criteriaFun
         end if

         if (criteriaFun <= this%TolFun) then
            convergenz = .true.
            x_sol      = xk
            return
         else
            dk = - solve(dFdx_val, F_val, this%lin_method)
            alphak = 1.0_rk
            xk = xk + alphak*dk
            k = k + 1
         end if
      end do
   end subroutine modified_newton_method_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine quasi_fd_newton_method_T1(this, F, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            real(rk), dimension(:), intent(in)  :: x
            real(rk), dimension(:), allocatable :: Fun
         end function Fun
      end interface

      procedure(Fun)  :: F

      class(nlsolver),               intent(inout) :: this
      real(rk), dimension(:),        intent(in)    :: x0
      real(rk), dimension(size(x0)), intent(out)   :: x_sol
      real(rk), dimension(size(x0))                :: xk
      real(rk), dimension(:),   allocatable        :: F_val
      real(rk), dimension(:,:), allocatable        :: dFdx_val
      real(rk)                                     :: criteriaFun
      integer                                      :: k
      logical                                      :: convergenz
      real(rk), dimension(size(x0))                :: dk
      real(rk)                                     :: alphak

      k          = 0
      xk         = x0
      convergenz = .false.

      do while (.not. convergenz .and. k < this%maxit)
         F_val    = F(xk)
         dFdx_val = derivative(f=F, x=xk, h=this%fdm_tol, method=this%fdm_method)

         criteriaFun = norm2(F_val)

         if (this%verbosity == 1) then
            print '(g0, e12.4)', k, criteriaFun
         end if

         if (criteriaFun <= this%TolFun) then
            convergenz = .true.
            x_sol      = xk
            return
         else
            dk = - solve(dFdx_val, F_val, this%lin_method)
            alphak = 1.0_rk
            xk = xk + alphak*dk
            k = k + 1
         end if
      end do
   end subroutine quasi_fd_newton_method_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine modified_quasi_fd_newton_method_T1(this, F, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            real(rk), dimension(:), intent(in)  :: x
            real(rk), dimension(:), allocatable :: Fun
         end function Fun
      end interface

      procedure(Fun)  :: F

      class(nlsolver),               intent(inout) :: this
      real(rk), dimension(:),        intent(in)    :: x0
      real(rk), dimension(size(x0)), intent(out)   :: x_sol
      real(rk), dimension(size(x0))                :: xk
      real(rk), dimension(:),   allocatable        :: F_val
      real(rk), dimension(:,:), allocatable        :: dFdx_val
      real(rk)                                     :: criteriaFun
      integer                                      :: k
      logical                                      :: convergenz
      real(rk), dimension(size(x0))                :: dk
      real(rk)                                     :: alphak

      k          = 0
      xk         = x0
      convergenz = .false.

      do while (.not. convergenz .and. k < this%maxit)
         F_val = F(xk)

         if ((mod(k, this%nmp) == 0)) dFdx_val = derivative(f=F, x=xk, h=this%fdm_tol, method=this%fdm_method)

         criteriaFun = norm2(F_val)

         if (this%verbosity == 1) then
            print '(g0, e12.4)', k, criteriaFun
         end if

         if (criteriaFun <= this%TolFun) then
            convergenz = .true.
            x_sol      = xk
            return
         else
            dk = - solve(dFdx_val, F_val, this%lin_method)
            alphak = 1.0_rk
            xk = xk + alphak*dk
            k = k + 1
         end if
      end do
   end subroutine modified_quasi_fd_newton_method_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine quasi_cs_newton_method_T0(this, F, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            complex(rk), intent(in) :: x
            complex(rk)             :: Fun
         end function Fun
      end interface

      procedure(Fun)  :: F

      class(nlsolver), intent(inout) :: this
      complex(rk),     intent(in)    :: x0
      complex(rk),     intent(out)   :: x_sol
      complex(rk)                    :: xk
      complex(rk)                    :: F_val
      real(rk)                       :: dFdx_val
      real(rk)                       :: criteriaFun
      integer                        :: k
      logical                        :: convergenz
      real(rk)                       :: dk
      real(rk)                       :: alphak

      k          = 0
      xk         = x0
      convergenz = .false.

      do while (.not. convergenz .and. k < this%maxit)
         F_val    = F(xk)

         dFdx_val = derivative(f=F, x=real(xk,kind=rk), h=this%cs_tol)

         criteriaFun = abs(F_val)

         if (this%verbosity == 1) then
            print '(g0, f12.4, 4x, e12.4, 4x, e12.4)', k, real(xk, kind=rk), real(F_val, kind=rk), dFdx_val
         end if

         if (criteriaFun <= this%TolFun) then
            convergenz = .true.
            x_sol      = xk
            return
         else
            dk = - F_val / dFdx_val
            alphak = 1.0_rk
            xk = xk + alphak*dk
            k = k + 1
         end if
      end do
   end subroutine quasi_cs_newton_method_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine modified_quasi_cs_newton_method_T0(this, F, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            complex(rk), intent(in) :: x
            complex(rk)             :: Fun
         end function Fun
      end interface

      procedure(Fun)  :: F

      class(nlsolver), intent(inout) :: this
      complex(rk),     intent(in)    :: x0
      complex(rk),     intent(out)   :: x_sol
      complex(rk)                    :: xk
      complex(rk)                    :: F_val
      real(rk)                       :: dFdx_val
      real(rk)                       :: criteriaFun
      integer                        :: k
      logical                        :: convergenz
      real(rk)                       :: dk
      real(rk)                       :: alphak

      k          = 0
      xk         = x0
      convergenz = .false.

      do while (.not. convergenz .and. k < this%maxit)
         F_val    = F(xk)
         if (mod(k, this%nmp) == 0) dFdx_val = derivative(f=F, x=real(xk,kind=rk), h=this%cs_tol)

         criteriaFun = abs(F_val)

         if (this%verbosity == 1) then
            print '(g0, f12.4, 4x, e12.4, 4x, e12.4)', k, real(xk, kind=rk), real(F_val, kind=rk), dFdx_val
         end if

         if (criteriaFun <= this%TolFun) then
            convergenz = .true.
            x_sol      = xk
            return
         else
            dk = - F_val / dFdx_val
            alphak = 1.0_rk
            xk = xk + alphak*dk
            k = k + 1
         end if
      end do
   end subroutine modified_quasi_cs_newton_method_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine quasi_cs_newton_method_T1(this, F, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            complex(rk), dimension(:), intent(in)  :: x
            complex(rk), dimension(:), allocatable :: Fun
         end function Fun
      end interface

      procedure(Fun)  :: F

      class(nlsolver),                  intent(inout) :: this
      complex(rk), dimension(:),        intent(in)    :: x0
      complex(rk), dimension(size(x0)), intent(out)   :: x_sol
      complex(rk), dimension(size(x0))                :: xk
      complex(rk), dimension(:),   allocatable        :: F_val
      real(rk),    dimension(:,:), allocatable        :: dFdx_val
      real(rk)                                        :: criteriaFun
      integer                                         :: k
      logical                                         :: convergenz
      real(rk), dimension(size(x0))                   :: dk
      real(rk)                                        :: alphak

      k          = 0
      xk         = x0
      convergenz = .false.

      do while (.not. convergenz .and. k < this%maxit)
         F_val    = F(xk)
         dFdx_val = derivative(f=F, x=real(xk,kind=rk), h=this%cs_tol)

         criteriaFun = norm2(real(F_val, kind=rk))

         if (this%verbosity == 1) then
            print '(g0, e12.4)', k, criteriaFun
         end if

         if (criteriaFun <= this%TolFun) then
            convergenz = .true.
            x_sol      = xk
            return
         else
            dk = - solve(dFdx_val, real(F_val, kind=rk), this%lin_method)
            alphak = 1.0_rk
            xk = xk + alphak*dk
            k = k + 1
         end if
      end do
   end subroutine quasi_cs_newton_method_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure subroutine modified_quasi_cs_newton_method_T1(this, F, x0,  x_sol)
      interface
         impure function Fun(x)
            import rk
            complex(rk), dimension(:), intent(in)  :: x
            complex(rk), dimension(:), allocatable :: Fun
         end function Fun
      end interface

      procedure(Fun)  :: F

      class(nlsolver),                  intent(inout) :: this
      complex(rk), dimension(:),        intent(in)    :: x0
      complex(rk), dimension(size(x0)), intent(out)   :: x_sol
      complex(rk), dimension(size(x0))                :: xk
      complex(rk), dimension(:),   allocatable        :: F_val
      real(rk),    dimension(:,:), allocatable        :: dFdx_val
      real(rk)                                        :: criteriaFun
      integer                                         :: k
      logical                                         :: convergenz
      real(rk), dimension(size(x0))                   :: dk
      real(rk)                                        :: alphak

      k          = 0
      xk         = x0
      convergenz = .false.

      do while (.not. convergenz .and. k < this%maxit)
         F_val    = F(xk)
         if (mod(k, this%nmp) == 0) dFdx_val = derivative(f=F, x=real(xk,kind=rk), h=this%cs_tol)

         criteriaFun = norm2(real(F_val, kind=rk))

         if (this%verbosity == 1) then
            print '(g0, e12.4)', k, criteriaFun
         end if

         if (criteriaFun <= this%TolFun) then
            convergenz = .true.
            x_sol      = xk
            return
         else
            dk = - solve(dFdx_val, real(F_val, kind=rk), this%lin_method)
            alphak = 1.0_rk
            xk = xk + alphak*dk
            k = k + 1
         end if
      end do
   end subroutine modified_quasi_cs_newton_method_T1
   !===============================================================================

end module forsolver
