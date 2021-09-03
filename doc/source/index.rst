.. _Arb: http://arblib.org

**Cascade** -- the C Library for Approximative Solutions to Complex Arbitrary Precision Differential Equations
==============================================================================================================

Cascade uses Arb_ to implement complex differential equations and to find solutions in the form of power series expansions. This expansion can then be utilized to compute analytic continuation along a piecewise linear path, which provides all the necessary tools to compute monodromies of differential equations with arbitrary precision.

Cascade's documentation
-----------------------

..  toctree::
    :maxdepth: 2

    acb_ode.rst
    acb_ode_solution.rst
    cascade.rst
