# Cascade - The C Library for Approximative Solutions to Complex Arbitrary Precision Differential Equations

Welcome to Cascade v0.5.3.
Cascade is a library designed to store and solve differential equations to arbitrary precision. This is accomplished through the use of Arblib's acb data type, which uses ball arithmetic to store arbitrary precision floats with error bounds. Solutions are computed as power series expanions around the origin using a recursion relation between the coefficients. With the help of analytic continuation this can be turned into a solution anywhere in the complex plane.

## Installation

This library can be build from source as a shared object library through cmake. To build the library, open a terminal in this directory and run the commands

```bash
cmake ./
make
```

## Jade

Along with Cascade, this project comes with Jade, an interface to Cascade, which can be used from the [Julia](https://julialang.org) command line. To use Jade, place the shared library file created by cmake  in 
```bash
/path/to/julia/lib/julia/libcascade.so
```
Additionally, either place the file Jade.jl somewhere where Julia can find it or add it to the LOAD_PATH via

```bash
push!(LOAD_PATH, "/path/to/jade")
```
Notice that only the directory path has to be given, not the file itself!

## Examples

The easiest way to get started, is through Julia. The following example creates an *acb_ode* storing Legendre's differential equation for $n = 4$. Then a solution is computed with a precision of 128 bits (which is equivalent to 38 decimal digits) and printed to the screen.

```julia
using Jade
S = ComplexField(128);
R = acbPolyRing(S,:z);
A = acb\_ode\_legendre(4);
@time p = powerSeries(A,S(1));
print(p*8)
```
The output of the above program should look something like this:

```bash
0.000117 seconds (8 allocations: 864 bytes)
[ 3.00000000000000000000000000000000000000 + i\*0, 0 + i\*0,
[-30.0000000000000000000000000000000000000 +/- 1e-42] + i\*0,
0 + i\*0,
[35.0000000000000000000000000000000000000 +/- 1e-42] + i\*0 ]
```
## Dependencies

Cascade uses [Arb](https://arblib.org) to store complex numbers and [Flint](http://flintlib.org) to handle memory management. Therefore both of these libraries have to be installed in order to build Cascade.

Jade uses [Nemo](https://nemocas.org), which adds a wrapper for Arblib to Julia. To use Jade, you must have Nemo installed.
