.. _Arb: http://arblib.org

**CASCADE** - the C Library for Approximative Solutions to Complex Arbitrary precision Differential Equations
==============================================================================================================

Cascade_ uses Arb_ to implement complex differential equations and to find solutions in the form of power series expansions. This expansion can then be utilized to compute analytic continuation along a piece-wise linear path, which provides all the necessary tools to compute monodromies of differential equations with arbitrary precision.

Cascade's documentation
=======================

..  toctree::
    :maxdepth: 2

    diffeq.rst
    cascade.rst

Jade's documentation
====================

This project also includes Jade_, the Julia interface to Arbitrary precision Differential Equations. It provides an intuitive interface to CASCADE, which can be used from the julia command line.

..  toctree::
    :maxdepth: 2

    jade.rst
