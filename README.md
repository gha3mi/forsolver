[![GitHub](https://img.shields.io/badge/GitHub-ForSolver-blue.svg?style=social&logo=github)](https://github.com/gha3mi/forsolver)
[![Version](https://img.shields.io/github/release/gha3mi/forsolver.svg)](https://github.com/gha3mi/forsolver/releases/latest)
[![Documentation](https://img.shields.io/badge/ford-Documentation%20-blueviolet.svg)](https://gha3mi.github.io/forsolver/)
[![License](https://img.shields.io/github/license/gha3mi/forsolver?color=green)](https://github.com/gha3mi/forsolver/blob/main/LICENSE)
[![Build](https://github.com/gha3mi/forsolver/actions/workflows/CI_test.yml/badge.svg)](https://github.com/gha3mi/forsolver/actions/workflows/CI_test.yml)


**ForSolver**: A Fortran library of linear and nonlinear solvers.

## Usage

### Linear system solver

```fortran
use forsolver, only: solve

x = solve(A,b,method)
```

available methods (optional):

- ```gesv```
- ```gels```

### Nonlinear system solver

```fortran
use forsolver, only: nlsolver

call nls%set_options(&
      lin_method,&
      nl_method,&
      fdm_method,&
      fdm_tol,&
      cs_tol,&
      TolFun,&
      maxit,&
      nmp,&
      verbosity )

call nls%solve(F, dFdx, x0, x_sol)
```

available nl_methods:

- ```newton```
- ```newton-modified```
- ```newton-quasi-fd```
- ```newton-quasi-fd-modified```
- ```newton-quasi-cs```
- ```newton-quasi-cs-modified```

fd: finite difference method

cs: complex step method

## Requirements

- A Fortran Compiler
- BLAS Library
- Fortran Package Manager (fpm)

## fpm Dependency

If you want to use `ForSolver` as a dependency in your own fpm project,
you can easily include it by adding the following line to your `fpm.toml` file:

```toml
[dependencies]
forsolver = {git="https://github.com/gha3mi/forsolver.git"}
```

## Runing Tests

Execute the following commands to run tests with specific compilers:

```shell
fpm @<compiler>-test
```
`compiler: ifx, ifort, gfortran, nvfortran`


## Examples

### Example 1: Linear System Solver

```fortran
program example1

   use kinds
   use forsolver

   implicit none

   real(rk), dimension(:,:), allocatable :: A
   real(rk), dimension(:)  , allocatable :: x, b
   integer                               :: m,n, i, j

   m = 3
   n = 2

   allocate(A(m,n),b(m),x(n))

   A(1,:) = [ 1.0_rk, 5.0_rk]
   A(2,:) = [ 3.0_rk, 1.0_rk]
   A(3,:) = [-2.0_rk, 4.0_rk]

   b = [4.0_rk, -2.0_rk, 3.0_rk]

   x = solve(A, b)

end program example1
```

### Example 2: Newton's Method for Root Finding

```fortran
module my_function3
   use kinds
   implicit none
contains
   function F1(x) result(F_val)
      real(rk), intent(in) :: x
      real(rk) :: F_val
      F_val = 5.0_rk * x**3 + 8.0_rk * x - 5.0_rk
   end function F1
   function dF1dx(x) result(dFdx_val)
      real(rk), intent(in) :: x
      real(rk) :: dFdx_val
      dFdx_val = 15.0_rk * x**2 + 8.0_rk
   end function dF1dx
end module my_function3

program example2

   use forsolver
   use my_function3

   implicit none

   type(nlsolver) :: nls
   real(rk)       :: x, expected_x

   call nls%set_options(&
      nl_method   = 'newton',&
      maxit       = 100,&
      TolFun      = 1e-4_rk,&
      verbosity   = 1)

   call nls%solve(F=F1, dFdx=dF1dx, x0=10.0_rk, x_sol=x)

end program example2
```

## API Documentation

The most up-to-date API documentation for the master branch is available
[here](https://gha3mi.github.io/ForSolver/).
To generate the API documentation for `ForSolver` using
[ford](https://github.com/Fortran-FOSS-Programmers/ford) run the following
command:

```shell
ford ford.yml
```

## Contributing

Contributions to `ForSolver` are welcome!
If you find any issues or would like to suggest improvements, please open an issue.