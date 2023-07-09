# ForSolver
ForSolver is a Fortran library that provides numerical solvers for linear and non-linear systems.

-----

## Table of Contents

- [ForSolver](#forsolver)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
    - [fpm](#fpm)
  - [Linear system solver](#linear-system-solver)
  - [Nonlinear system solver](#nonlinear-system-solver)
  - [Tests](#tests)
  - [Examples](#examples)
    - [Example 1: Linear System Solver](#example-1-linear-system-solver)
    - [Example 2: Newton's Method for Root Finding](#example-2-newtons-method-for-root-finding)
  - [Documentation](#documentation)
  - [Contributing](#contributing)
-----

## Installation

### fpm
ForSolver can be cloned and then built using [fpm](https://github.com/fortran-lang/fpm), following the instructions provided in the documentation available on Fortran Package Manager.

```bash
git clone https://github.com/gha3mi/forsolver.git
cd forsolver
fpm install --perfix .
```

Or you can easily include this package as a dependency in your `fpm.toml` file.

```toml
[dependencies]
forsolver = {git="https://github.com/gha3mi/forsolver.git"}
```
-----

## Linear system solver
```fortran
use forsolver, only : solve
x = solve(A,b,method)
```
available methods (optional):
- ```gesv```
- ```gels```
-----

## Nonlinear system solver
```fortran
use forsolver, only : nlsolver

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

-----

## Tests

The `tests` directory contains test programs to verify the functionality of the `forsolver` module. To run the tests using `fpm`, you can use response files for specific compilers:

- For Intel Fortran Compiler (ifort):
```bash
fpm @ifort
```

- For Intel Fortran Compiler (ifx):
```bash
fpm @ifx
```

- For NVIDIA Compiler (nvfortran):
```bash
fpm @nvidia
```

- For GNU Fortran Compiler (gfortran):
```bash
fpm @gfortran
```
-----

## Examples

### Example 1: Linear System Solver

```fortran
program test1

   use kinds
   use forsolver, only : solve

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

   X = solve(A, b)

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
```

### Example 2: Newton's Method for Root Finding

```fortran
program test3

   use kinds
   use functions_module
   use forsolver, only : nlsolver

   implicit none

   type(nlsolver) :: nls
   real(rk)       :: x_sol

   call nls%set_options(&
      nl_method   = 'newton',&
      maxit       = 100,&
      TolFun      = 1e-15_rk,&
      verbosity   = 1)

   call nls%solve(F=F1, dFdx=dF1dx, x0=10.0_rk, x_sol=x_sol)

end program test3
```
-----

## Documentation
To generate the documentation for the `ForSolver` module using [ford](https://github.com/Fortran-FOSS-Programmers/ford) run the following command:
```bash
ford ford.yml
```
-----

## Contributing

Contributions to `ForSolver` are welcome! If you find any issues or would like to suggest improvements, please open an issue or submit a pull request.

