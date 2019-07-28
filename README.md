# Cascade - The C Library for Approximative Solutions to Complex Arbitrary Precision Differential Equations

Welcome to Cascade v0.7.
Cascade is a library designed to store and solve differential equations to arbitrary precision. This is accomplished through the use of [Arblib](https://arblib.org)'s acb data type, which uses ball arithmetic to store arbitrary precision floats with error bounds. Solutions are computed as power series expanions around the origin using a recursion relation between the coefficients. With the help of analytic continuation this can be turned into a solution anywhere in the complex plane.

Author: Matthias Gessinger

## Installation

This library can be build from source as a shared object library through cmake. To build the library, open a terminal in the ./build/ directory and run the commands

```bash
cmake ./
make
make install
```
Depending on your system, *make install* might have to be run with root priviliges.

## Jade

Along with Cascade, this project comes with Jade, an interface to Cascade, which can be used from the [Julia](https://julialang.org) command line. To use Jade, place the file Jade.jl somewhere where Julia can find it or add it to the LOAD_PATH via

```bash
push!(LOAD_PATH, "/path/to/jade")
```
Note that only the directory path has to be given, not the file itself!

## Examples

The easiest way to get started, is through Julia. The following example creates an *acb_ode* storing Legendre's differential equation for *n = 4*. Then a solution is computed with a precision of 256 bits (which is equivalent to 77 decimal digits) and printed to the screen.

```julia
using Jade
S = ComplexField(256);
R = AcbPolyRing(S,:z);
A = acb_ode_legendre(R,4);
@time p = powerSeries(A,S(1));
print(p)
```
The output of the above program should look something like this:

```bash
  0.000117 seconds (8 allocations: 864 bytes)
[ 0.37500000000000000000000000000000000000000000000000000000000000000000000000000 + i*0, 0 + i*0,
[-3.75000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 1e-81] + i*0, 0 + i*0, 
[ 4.37500000000000000000000000000000000000000000000000000000000000000000000000000 +/- 1e-81] + i*0 ]
```

Equivalently we can run the above program through Cascade directly. The code to perform the same computation as above, looks like this:
```C
#include <cascade.h>

int main (int argc, char **argv)
{
    acb_ode_t ODE = acb_ode_legendre(4);
    find_power_series(ODE,NULL,256);
    acb_ode_dump(ODE,NULL);
    acb_ode_clear(ODE);
    return 0;
}
```

To compile the program, run
```bash
gcc -o test test.c -lcascade
```

The output of *time ./test* should then look something like this:
```bash
Order: 2
Degree: 2
diff_eq_poly(ODE,0) = (20 + 0j)  +/-  (0, 0j)	(0 + 0j)  +/-  (0, 0j)	(0 + 0j)  +/-  (0, 0j)	
diff_eq_poly(ODE,1) = (0 + 0j)  +/-  (0, 0j)	(-2 + 0j)  +/-  (0, 0j)	(0 + 0j)  +/-  (0, 0j)	
diff_eq_poly(ODE,2) = (1 + 0j)  +/-  (0, 0j)	(0 + 0j)  +/-  (0, 0j)	(-1 + 0j)  +/-  (0, 0j)	

Solution:
[(0.375 + 0j)  +/-  (0, 0j)
(0 + 0j)  +/-  (0, 0j)
(-3.75 + 0j)  +/-  (2.8e-154, 0j)
(0 + 0j)  +/-  (0, 0j)
(4.375 + 0j)  +/-  (3.26e-154, 0j)]
real	0m0,008s
user	0m0,001s
sys	0m0,007s
```

## Memory management

Because Arb caches some constants internally, it is recommended to call *flint_cleanup()* at the end of your main program. This will clear Arb's internal cache and guarantee a clean output of *Valgrind*.

## Dependencies

Cascade uses [Arb](https://arblib.org) to store complex numbers and [Flint](http://flintlib.org) to handle memory management. Therefore both of these libraries have to be installed in order to build Cascade.

Jade uses [Nemo](https://nemocas.org), which adds a wrapper for Arblib to Julia. To use Jade, you must have Nemo installed.
