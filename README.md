# ForSolver
ForSolver is a Fortran library that provides numerical solvers for linear and non-linear systems.

-----

## Table of Contents

- [ForSolver](#forsolver)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
    - [fpm](#fpm)
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
module functions_module
use kinds
implicit none

contains

function F(x) result(F_val)
   real(rk), intent(in) :: x
   real(rk) :: F_val
   F_val = 5_rk * x**3 + 8_rk * x - 5_rk
end function F

function dFdx(x) result(dFdx_val)
   real(rk), intent(in) :: x
   real(rk) :: dFdx_val
   dFdx_val = 15_rk * x**2 + 8_rk
end function dFdx

end module functions_module


program test2
   use kinds
   use forsolver
   use functions_module
   implicit none
   real(rk) :: x0, tol, x_sol
   integer :: maxit

   ! Variable declaration
   x0    = 10.0_rk
   tol   = 1e-8_rk
   maxit = 100

   x_sol = solve(F, dFdx, x0, tol, maxit)

end program test2
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

