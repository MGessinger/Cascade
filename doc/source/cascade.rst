.. _Cascade:

**cascade.h** -- Solving complex differential equations
==================================================================================

Cascade implements complex differential equations in a variable of type :type:`acb_ode_struct`. It provides functionality for finding power series solutions around the origin, but is also capable of performing analytic continuation along a piecewise linear path to find power series expansions anywhere in the complex plane - provided the solution is analytic in at least a small neighborhood of the origin.

Special Equations
------------------

.. function:: acb_ode_t acb_ode_legendre (ulong n)

    Returns an `acb_ode_t` initialized to Legendre's equation parametrized by a natural number *n*. The initial values are set in such a way that the computed solutions yield the Legendre functions of first kind. Because these solutions are polynomial, this is a good starting point for working with Cascade. Legendre's equation is given by:

    .. math::
        (1-x^2)y'' - 2xy' + n(n+1)y = 0.

.. function:: acb_ode_t acb_ode_bessel (acb_t nu, slong bits)

    Returns an `acb_ode_t` initialized to Legendre's equation parametrized by a complex number *nu* with a precision of *bits*. The initial values are set in such a way that for integral nu, the Bessel functions of first kind are obtained. Because these have maximally unipotent monodromy around the origin, they are a good test for computing monodromy matrices. Bessel's equation is given by:

    .. math::
        x^2y'' + xy' + (x^2 - \nu^2)y = 0

.. function:: acb_ode_t acb_ode_hypgeom (acb_t a, acb_t b, acb_t c, slong bits)

    Returns an `acb_ode_t` initialized to Euler's hypergeometric equation. The initial values are set in such a way that the series expansion yields the hypergeometric series 2F1. Euler's equation is given by:

    .. math::
        z(1-z)y'' + (c - (a + b + 1)z)y' - aby = 0

Functions
------------------

.. function:: ulong find_power_series (acb_ode_t ODE, acb_t in, slong bits)

    Finds a solution to *ODE*, by recursively computing the coefficients of the power series expansion. The iteration stops when the summands, given by *a_k\*in^k*, are small enough.  This function does not perform any tests for convergence beforehnd, but it aborts if too many coefficients are computed.

.. function:: void analytic_continuation (acb_t res, acb_ode_t ODE, acb_srcptr path, slong len, slong prec, int output_series)

    Perform analytic continuation along *path*, which stores the *len* cornes of a piecewise linear path in the complex plane. This is implemented by computing a power series expansion at each corner and then transforming the origin. The data stored inside *ODE* remains unchanged.

.. function:: void find_monodromy_matrix (acb_mat_t monodromy, acb_ode_t ODE, acb_t z0, slong bits)

    Compute the monodromy matrix of *ODE* and store it in *monodromy*. This is implemented by performing `analytic_continuation` for multiple different initial conditions. Currently the path is implemented as a polygon with 32 corners. The radius of the polygon is chosen by calling :func:`radiusOfConvergence`.

.. function:: int checkODE (acb_poly_t *polys, acb_ode_t ODE, acb_t z, slong bits)

    Plugs the solution stored in *ODE* into the actual equation and checks if the result actually vanishes. If not, *ODE* is dumped by `acb_ode_dump`.

.. function:: void radiusOfConvergence(arb_t rad, acb_ode_t ODE, slong bits)

    Sets *rad* to the radius of a ball that is garantueed not to contain any singularities of *ODE* other than zero. If zero is not a singularity of ODE, *rad* is set to zero. This is computed using Rouche's theorem applied to the leading polynomial of *ODE*. If the leading polynomial of ODE is of the form z^r*p(z) with p not vanishing at zero, then *rad* is chosen in such a way that 
    
    .. math:: |a_0| < \sum_{k=1}^{n} |a_k| rad^k

    holds, where a_k are the coefficients of p and deg(p) = n. By Rouche's theorem, this implies that *B(0,rad)* contains no roots of p.
