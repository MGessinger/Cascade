.. _Cascade:

**cascade.h** -- Solving complex differential equations
==================================================================================

Cascade implements complex differential equations in a variable of type :type:`acb_ode_struct`. It provides functionality for finding power series solutions around the origin, but is also capable of performing analytic continuation along a piecewise linear path to find power series expansions anywhere in the complex plane - provided the solution is analytic in at least a small neighborhood of the origin.

Special Equations
------------------

.. function:: void acb_ode_legendre (acb_ode_t ODE, ulong n)

	Initializes *ODE* to Legendre's equation parametrized by a natural number *n*. Because their solutions are polynomial (provided the correct initial values), this is a good starting point for working with Cascade. Legendre's equation is given by:

	.. math::
		(1-x^2)y'' - 2xy' + n(n+1)y = 0.

.. function:: void acb_ode_bessel (acb_ode_t ODE, acb_t nu, slong bits)

	Initializes *ODE* to Bessel's equation parametrized by a complex number *nu* with a precision of *bits*.  Bessel's equation is given by:

	.. math::
		x^2y'' + xy' + (x^2 - \nu^2)y = 0.

.. function:: void acb_ode_hypgeom (acb_ode_t ODE, acb_t a, acb_t b, acb_t c, slong bits)

	Initializes *ODE* to Euler's hypergeometric equation. Euler's equation is given by:

	.. math::
		z(1-z)y'' + (c - (a + b + 1)z)y' - aby = 0.

Functions
------------------

.. function:: slong truncation_bound (arb_t eta, arb_t rho, slong bits)

	Return a bound for the number of coefficients necessary to assert the Taylor expansion to have a tail of less than 2^-*bits*.

.. function:: void acb_ode_solve_fuchs (acb_poly_t res, acb_ode_t ODE, slong num_of_coeffs, slong bits)

	Computes a truncated power series solution of the differential equation defined by ODE, using the initial values provided in *res*, and truncated to degree *num_of_coeffs*. This number could (and in general should) be obtained by calling :func:`truncation_bound`.

.. function:: void analytic_continuation (acb_t res, acb_ode_t ODE, acb_srcptr path, slong len, slong num_of_coeffs, slong prec, int output_series)

	Performs analytic continuation along *path*, which stores the *len* cornes of a piecewise linear path in the complex plane. This is implemented by computing a power series expansion of degree *num_of_coeffs* at each corner and then transforming the origin. The data stored inside *ODE* remains unchanged.

.. function:: void find_monodromy_matrix (acb_mat_t monodromy, acb_ode_t ODE, slong bits)

	Compute the monodromy matrix of *ODE* and store it in *monodromy*. This is implemented by performing `analytic_continuation` for multiple different initial conditions. Currently the path is implemented as a polygon with 256 corners. The radius of the polygon is chosen by calling :func:`radius_of_convergence`.

.. function:: void radius_of_convergence (arb_t rad, acb_ode_t ODE, slong bits)

	Sets *rad* to the radius of a ball that is garantueed not to contain any singularities of *ODE* other than (possibly) zero. This is computed by bounding the inverse roots of the leading polynomial using Fujiwara's bound:

	.. math::
		R < 2 min\{\left| \frac{a_n}{a_k} \right| ^{1/k} \mid a_k \neq 0 \}.

	The inverse of this bound then yields a lower bound on the distance to the nearest singular point of *ODE*.

.. function:: void acb_poly_graeffe_transform (acb_ptr dest, acb_srcptr src, slong len, slong bits)

	Computes the Graeffe Transform of src and stores the result in dest. Aliasing is allowed.
