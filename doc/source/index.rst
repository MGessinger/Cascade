.. _Arb: http://arblib.org
.. _Julia: https://julialang.org

**Cascade** -- the C Library for Approximative Solutions to Complex Arbitrary Precision Differential Equations
==============================================================================================================

Cascade uses Arb_ to implement complex differential equations and to find solutions in the form of power series expansions. This expansion can then be utilized to compute analytic continuation along a piecewise linear path, which provides all the necessary tools to compute monodromies of differential equations with arbitrary precision.

Cascade's documentation
-----------------------

..  toctree::
    :maxdepth: 2

    diffeq.rst
    cascade.rst

Jade's documentation
--------------------

This project also includes :ref:`Jade`. It provides an intuitive interface to Cascade, which can be used from the Julia_ command line.

..  toctree::
    :maxdepth: 2

    jade.rst
