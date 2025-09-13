module external_interfaces_solver

    use forsolver_kinds, only: rk

    implicit none

    interface gesv
#if defined(REAL64)
        pure subroutine dgesv(fn, fnrhs, fa, flda, fipiv, fb, fldb, finfo)
            import rk
            implicit none
            integer,  intent(in)    :: fn, fnrhs, flda, fldb
            real(rk), intent(inout) :: fa(flda,fn), fb(fldb,fnrhs)
            integer,  intent(out)   :: finfo
            integer,  intent(out)   :: fipiv(fn)
        end subroutine dgesv
#elif defined(REAL32)
        pure subroutine sgesv(fn, fnrhs, fa, flda, fipiv, fb, fldb, finfo)
            import rk
            implicit none
            integer,  intent(in)    :: fn, fnrhs, flda, fldb
            real(rk), intent(inout) :: fa(flda,fn), fb(fldb,fnrhs)
            integer,  intent(out)   :: finfo
            integer,  intent(out)   :: fipiv(fn)
        end subroutine sgesv
#else
        pure subroutine dgesv(fn, fnrhs, fa, flda, fipiv, fb, fldb, finfo)
            import rk
            implicit none
            integer,  intent(in)    :: fn, fnrhs, flda, fldb
            real(rk), intent(inout) :: fa(flda,fn), fb(fldb,fnrhs)
            integer,  intent(out)   :: finfo
            integer,  intent(out)   :: fipiv(fn)
        end subroutine dgesv
#endif
    end interface

    interface gels
#if defined(REAL64)
        pure subroutine dgels(ftrans, fm, fn, fnrhs, fa, flda, fb, fldb, fwork, flwork, finfo)
            import :: rk
            character(len=1), intent(in)    :: ftrans
            integer,          intent(in)    :: fm, fn, fnrhs, flda, fldb, flwork
            real(rk),         intent(inout) :: fa(flda,*), fb(fldb,*)
            real(rk),         intent(in)    :: fwork(*)
            integer,          intent(out)   :: finfo
        end subroutine dgels
#elif defined(REAL32)
        pure subroutine sgels(ftrans, fm, fn, fnrhs, fa, flda, fb, fldb, fwork, flwork, finfo)
            import :: rk
            character(len=1), intent(in)    :: ftrans
            integer,          intent(in)    :: fm, fn, fnrhs, flda, fldb, flwork
            real(rk),         intent(inout) :: fa(flda,*), fb(fldb,*)
            real(rk),         intent(in)    :: fwork(*)
            integer,          intent(out)   :: finfo
        end subroutine sgels
#else
        pure subroutine dgels(ftrans, fm, fn, fnrhs, fa, flda, fb, fldb, fwork, flwork, finfo)
            import :: rk
            character(len=1), intent(in)    :: ftrans
            integer,          intent(in)    :: fm, fn, fnrhs, flda, fldb, flwork
            real(rk),         intent(inout) :: fa(flda,*), fb(fldb,*)
            real(rk),         intent(in)    :: fwork(*)
            integer,          intent(out)   :: finfo
        end subroutine dgels
#endif
    end interface

end module external_interfaces_solver
