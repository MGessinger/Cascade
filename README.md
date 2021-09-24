# Cascade - The C Library for Approximative Solutions to Complex Arbitrary Precision Differential Equations

Welcome to Cascade v2.1.
Cascade is a library designed to store and solve differential equations to arbitrary precision.
This is accomplished through the use of [Arblib](https://arblib.org)'s `acb_t` data type, which uses ball arithmetic to store arbitrary precision floating point numbers with error bounds.
Solutions are computed as power series expansions about the origin using a recursion relation between the coefficients.
With the help of analytic continuation this can be turned into a solution anywhere in the complex plane.

For a similar library using p-adic arithmetic, see [Implode](https://github.com/MGessinger/implode).

Author: Matthias Gessinger

## Installation

This library can be build from source as a shared object library through CMake.
For a default installation, run the following commands:

```bash
mkdir build
cd build
make
make test
```
If all tests pass, finally install the library by running
```bash
sudo make install
```

## Examples

```C
#include <acb_poly.h>
#include <cascade.h>

int main ()
{
    acb_t z;
    acb_poly_t pol;
    acb_ode_t ODE;

    acb_init(z);
    acb_poly_init(pol);
    acb_ode_legendre(ODE, 4);

    acb_set_d(z, 0.375);
    acb_poly_set_coeff_acb(pol, 0, z);
    acb_ode_solve_fuchs(pol, ODE, 5, 1024);
    acb_ode_dump(ODE, NULL);
    acb_poly_printd(pol, 10);

    acb_clear(z);
    acb_poly_clear(pol);
    acb_ode_clear(ODE);
    flint_cleanup();
    return 0;
}
```

To compile the program, run
```bash
gcc test.c -lcascade -larb -lflint
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

Because Arb caches some constants internally, it is recommended to call *flint_cleanup()* at the end of your main program.
This will clear Arb's internal cache and guarantee a clean output of *Valgrind*.

## Dependencies

Cascade uses [Arb](https://arblib.org) to store complex numbers and [Flint](http://flintlib.org) to handle memory management.
Therefore both of these libraries have to be installed in order to build Cascade, and also to build programs using Cascade!
