.. _Cascade:

Solvers for Differential Equations
==================================================================================

Cascade implements complex differential equations in a variable of type :type:`acb_ode_t`.
It provides functionality for finding power series solutions around the origin, but is also capable of performing analytic continuation along a piecewise linear path to find power series expansions anywhere in the complex plane - provided the solution is analytic in at least a small neighborhood of the origin.

Solvers
----------------------------------------------------------------------

.. function:: void acb_ode_solve_fuchs (acb_poly_t res, acb_ode_t L, slong deg, slong bits)

	Computes a truncated power series solution of the differential equation defined by *L*, using Fuchs' method.
	The initial values are provided in *res*, and the resulting series is truncated to length *deg*.
	This number could (and in general should) be obtained by calling :func:`truncation_bound`.

.. function:: void _acb_ode_solve_frobenius (acb_poly_t res, acb_ode_t L, acb_ode_solution_t rhs, slong sol_degree, slong prec)

	Compute a single solution to the linear, differential equation :math:`Ly = rhs`.
	The variable *rhs* is assumed to be initialized to multiplicity one, and to fit into one of the following two cases:

	 * The power series stored in rhs is zero, and the exponent :math:`\rho` solves the indicial equation
	 * The power series is non-zero, and none of the values :math:`\rho, \rho + 1, \rho + 2, \dots` solve the indicial equation.

	This assumption is not checked, and if neither of these cases are fulfilled, behaviour is undefined.
	When solving a homogeneous equation, it is generally recommended to call :func:`acb_ode_solve_frobenius` for an automatic algorithm choice.

.. function:: void acb_ode_solve_frobenius (acb_ode_solution_t sol, acb_ode_t L, slong deg, slong bits)

	Computes a set of generalized solutions to the homogeneous differential equation :math:`Ly = 0`.
	The variable *sol* is assumed to be initialized to correspond to a solution of the indicial equation.
	This is not confirmed, and if this is not the case, behaviour is undefined.

	The function computes all power series necessary to represent the solution(s) corresponding to the exponent stored in *sol*, truncated to length *deg*.
	This number could (and in general should) be obtained by calling :func:`truncation_bound`.

Helper Functions
----------------------------------------------------------------------

.. function:: slong truncation_bound (arb_t eta, arb_t rho, slong p)

	Return a bound for the number of coefficients necessary to assert the Taylor expansion to have a tail of less than :math:`2^{-p}`.
	The variable *rho* contains the radius of convergence of the power series, but since the error bound can't hold on the entire disk of radius *rho*, a smaller disk of radius *eta* is considered.

	It uses the folowing formula:

	.. math::
		n = \frac{\log(1 - \eta/\rho) - p \log(2)}{\log(\eta/\rho)}

.. function:: void analytic_continuation (acb_poly_t res, acb_ode_t L, acb_srcptr path, slong len, slong deg, slong bits)

	Performs analytic continuation along *path*, which stores the *len* corners of a piecewise linear path in the complex plane.
	This is implemented by computing a power series expansion of degree *deg* at each corner using the Fuchsian solver, and then transforming the origin.

.. function:: void find_monodromy_matrix (acb_mat_t mono, acb_ode_t L, slong bits)

	Compute the monodromy matrix of *L* and store it in *mono*.
	This is implemented by performing :func:`analytic_continuation` for multiple different initial conditions.
	Currently the path is implemented as a polygon with 256 corners.
	The radius of the polygon is chosen by calling :func:`radius_of_convergence`.

.. function:: void radius_of_convergence (arb_t rad, acb_ode_t L, slong bits)

	Sets *rad* to the radius of a ball that is garantueed not to contain any singularities of *L* other than (possibly) zero. This is computed by bounding the inverse roots of the leading polynomial using Fujiwara's bound:

	.. math::
		R < 2 min\{\left| \frac{a_n}{a_k} \right| ^{1/k} \mid a_k \neq 0 \}.

	The inverse of this bound then yields a lower bound on the distance to the nearest singular point of *L*.

.. function:: void indicial_polynomial (acb_poly_t f, acb_ode_t L, slong nu, slong s, slong prec);

	Compute the *nu*-th indical polynomial *f* defined by the differential operator L and shifted by *s*.
	This polynomial is defined as the coefficient of :math:`x^{\nu}` in

	.. math::
		f(\rho) = \rho (\rho-1) \dots (\rho - \lambda + 1) p_{\lambda}(x) + \dots + \rho p_1(x) + p_0(x)

	.. note::
		This polynomial is only defined for differential operators of Frobenius type.
		If *L* is not of this type, then the computation of *f* still works equally, but its value has little to no significance.

.. function:: void indicial_polynomial_evaluate (acb_t f, acb_ode_t L, slong nu, acb_t rho, slong s, slong prec);

	Compute the value of the *nu*-th indicial polynomial, defined as above, and evaluate it at the point :math:`\rho + s`.
