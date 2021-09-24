.. _solutions:

Solutions to Differential Equations
======================================================================

A variable of type :type:`acb_ode_solution_t` represents a solution to a differential equation of Frobenius type.
It stores the generalized exponent, as well as all power series necessary to define the solutions to this exponent.

Types and macros
----------------------------------------------------------------------

.. type:: acb_ode_solution_struct

.. type:: acb_ode_solution_t

	This data structure contains the defining data for a solution of a differential equation of Frobenius type.

	An `acb_ode_solution_t` is defined as an array of type `acb_ode_solution_struct` of length 1, so it can be passed by reference.

Memory management
----------------------------------------------------------------------

.. function:: void acb_ode_solution_init (acb_ode_solution_t sol, acb_t rho, slong mul, slong alpha)

	Initialize the variable *sol* to represent a solution corresponding to the generalized exponent *rho*.
	The variable *rho* is assumed to be a solution of the indicial equation of multiplicity *mul*.
	Furthermore, it is assumed that there are *alpha* roots of the indicial equation (with multiplicity), that differ from *rho* by a positive integer.

	Neither of these assumptions are confirmed.
	If either of them is untrue, behavior is undefined.

.. function:: void acb_ode_solution_clear (acb_ode_solution_t sol)

	Clears the memory allocated for *sol* by a previous call to `acb_ode_solution_init`.

	.. important::
		Because Arb caches certain values internally, it is recommended to call *flint_cleanup()* at the end of your main program.
		This will clear Arb's cache and lead to a clean output when using *Valgrind*.

Manipulation
----------------------------------------------------------------------

These functions are used internally by the solvers.
There should rarely be a need to call these directly.

.. function:: void _acb_ode_solution_update (acb_ode_solution_t sol, acb_poly_t f, slong prec)
	
	Replace the internal power series :math:`\sum g_i x^i` by :math:`\sum fg_i x^i`, and recalculate the derivatives of the :math:`g_i` with respect to *rho*.

.. function:: void _acb_ode_solution_extend (acb_ode_solution_t sol, slong nu, acb_poly_t g_nu, slong prec)

	Set the *nu*-th coefficients of the internal power series to the *rho*-derivatives of *g_nu*.

.. function:: void acb_ode_solution_evaluate (acb_t out, acb_ode_solution_t sol, acb_t a, slong prec)

	Evaluate the solution stored in *sol* at the point a.
	Because any solution may contain logarithms, a must not be zero.
