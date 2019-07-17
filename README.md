# Cascade - The C Library for Approximative Solutions to Complex Arbitrary Precision Differential Equations

Welcome to Cascade v0.5.3. Cascade is a library designed to store and solve differential equations to arbitrary precision. This is accomplished through the use of Arblib's acb data type, which uses ball arithmetic to store arbitrary precision floats with error bounds. Solutions are computed as power series expanions around the origin using a recursion relation between the coefficients. With the help of analytic continuation this can be turned into a solution anywhere in the complex plane.

## Installation

This library can be build from source as a shared object library through cmake. To build the library, open a terminal in this directory and run the commands

```bash
cmake ./
make
```

## Jade

Along with Cascade, this project comes with Jade, an interface to Cascade, which can be used from the [julia](https://julialang.org) command line. To use Jade in julia, either place the file jade.jl somewhere where julia can find it or add it to the LOAD_PATH via

```bash
push!(LOAD_PATH, "/path/to/jade")
```

Notice that only the directory path has to be given, not the file itself!

## Dependencies

Cascade uses [Arb](https://arblib.org) to store complex numbers and [Flint](http://flintlib.org) to handle memory management. Therefore both of these libraries have to be installed in order to build Cascade.

Jade uses [Nemo](https://nemocas.org), which adds a wrapper for arblib to julia. To use Jade, you must have Nemo installed.
