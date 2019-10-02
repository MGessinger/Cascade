# Cascade - The C Library for Approximative Solutions to Complex Arbitrary Precision Differential Equations

Welcome to Cascade v0.1.0.
Cascade is a library designed to store and solve differential equations to arbitrary precision. This is accomplished through the use of [Arblib](https://arblib.org)'s acb data type, which uses ball arithmetic to store arbitrary precision floats with error bounds. Solutions are computed as power series expanions around the origin using a recursion relation between the coefficients. With the help of analytic continuation this can be turned into a solution anywhere in the complex plane.

Author: Matthias Gessinger

## Installation

This library can be build from source as a shared object library through CMake. For a default installation, simply use the command

```bash
bash install.sh
```
to run the usual commands for executing CMake. If however any custom flags are desired, CMake will have to be run manually.
Depending on your system, the installation may have to be run with root privileges to be ablot to put all the files in the necessary places.

## Jade

Along with Cascade, this project comes with Jade, an interface to Cascade, which can be used from the [Julia](https://julialang.org) command line. To use Jade, place the file Jade.jl somewhere where Julia can find it or add it to the LOAD_PATH via

```bash
push!(LOAD_PATH, "/path/to/jade")
```
Note that only the directory path has to be given, not the file itself!

## Examples

The easiest way to get started, is through Julia. The following example creates an *acb_ode* storing Legendre's differential equation for *n = 4*. Then a solution is computed with a precision of 1024 bits (which is equivalent to 308 decimal digits) and printed to the screen.

```julia
using Jade
S = ComplexField(1024);
R = AcbPolyRing(S,:z);
A = acb_ode_legendre(R,4);
@time p = powerSeries(R(3//8),A,5)
```
The output of the above program should look something like this:

```bash
  0.007656 seconds (20 allocations: 2.281 KiB)
[ 0.37500000000000000000000000000000000000000000000000000000000000000000000000000 + i*0, 0 + i*0,
[-3.75000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 1e-313] + i*0, 0 + i*0, 
[ 4.37500000000000000000000000000000000000000000000000000000000000000000000000000 +/- 1e-313] + i*0 ]
```

Equivalently we can run the above program through Cascade directly. The code to perform the same computation as above, looks like this:
```C
#include <cascade.h>
#include <acb_poly.h>

int main (int argc, char **argv)
{
    acb_ode_t ODE = acb_ode_legendre(4);
    acb_poly_t pol;
    acb_poly_init(pol);
    acb_poly_one(pol);
    acb_set_d(acb_poly_get_coeff_ptr(pol,0),0.375);
    find_power_series(pol,ODE,5,1024);
    acb_ode_dump(ODE,NULL);
    acb_poly_printd(pol,10);
    acb_ode_clear(ODE);
    return 0;
}
```

To compile the program, run
```bash
gcc test.c -lcascade
```

The output of *time ./a.out* should then look something like this:
```bash
Order: 2
Degree: 2
diff_eq_poly(ODE,0) = (20 + 0j)  +/-  (0, 0j)	(0 + 0j)  +/-  (0, 0j)	(0 + 0j)  +/-  (0, 0j)	
diff_eq_poly(ODE,1) = (0 + 0j)  +/-  (0, 0j)	(-2 + 0j)  +/-  (0, 0j)	(0 + 0j)  +/-  (0, 0j)	
diff_eq_poly(ODE,2) = (1 + 0j)  +/-  (0, 0j)	(0 + 0j)  +/-  (0, 0j)	(-1 + 0j)  +/-  (0, 0j)	

Solution:
[(0.375 + 0j)  +/-  (0, 0j)
(0 + 0j)  +/-  (0, 0j)
(-3.75 + 0j)  +/-  (2.8e-308, 0j)
(0 + 0j)  +/-  (0, 0j)
(4.375 + 0j)  +/-  (3.26e-313, 0j)]
real	0m0,009s
user	0m0,005s
sys	0m0,005s
```

## Memory management

Because Arb caches some constants internally, it is recommended to call *flint_cleanup()* at the end of your main program. This will clear Arb's internal cache and guarantee a clean output of *Valgrind*.

## Dependencies

Cascade uses [Arb](https://arblib.org) to store complex numbers and [Flint](http://flintlib.org) to handle memory management. Therefore both of these libraries have to be installed in order to build Cascade.

Jade uses [Nemo](https://nemocas.org), which adds a wrapper for Arblib to Julia. To use Jade, you must have Nemo installed.
