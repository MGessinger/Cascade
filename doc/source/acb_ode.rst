.. _acb-ode:

**acb_ode.h** -- Complex Differential Equations
========================================================================

An :type:`acb_ode_t` represents an ordinary, linear differential operator with complex polynomial coefficients and error bounds. The polynomials are implemented as an array of type :type:`acb_struct`. The order of the differential equation as well as the polynomial degrees are tracked internally and need no further care by the user.

Types, macros and constants
------------------------------
.. type:: acb_ode_struct

.. type:: acb_ode_t

	Contains a pointer to an array of the coefficients, the maximum degree of all the polynomials and the order of the ODE.

	An `acb_ode_t` is defined as a pointer of type `acb_ode_struct`, permitting an `acb_ode_t` to be passed by reference.

.. macro:: diff_eq_poly (ODE, i)

	Macro returning a pointer to the polynomial belonging to the *i-th* derivative.

.. macro:: diff_eq_coeff (ODE, i ,j)

	Macro returning a pointer to the *j-th* coefficient of the *i-th* polynomial. This macro is identical to *diff_eq_poly(ODE,i)->(j)*.

.. macro:: order (ODE)

	Returns the order of the differential, which is one less than the number of defining polynomials.

.. macro:: degree (ODE)

	Returns the maximum degree amongst the defining polynomials. This is computed through *acb_poly_degree* and therefore the same restrictions apply in the case of inexact polynomials.

Memory management
------------------------------------------------------------------------

.. function:: acb_ode_t acb_ode_init_blank (slong degree, slong order)

	Returns a pointer of type `acb_ode_t`, which contains space for *order+1* polynomials of length no more than *degree+1*.

.. function:: acb_ode_t acb_ode_init (acb_poly_t polys, slong order)

	Returns a pointer of type `acb_ode_t`, which has been initialized to contain *polys*.

.. function:: void acb_ode_clear (acb_ode_t ODE)

	Clears the memory allocated by a previous call to `acb_ode_init`.

	.. important::
		Because Arb caches certain values internally, it is recommended to call *flint_cleanup()* at the end of your main program. This will clear Arb's cache and lead to a clean output when using *Valgrind*.

Manipulation
------------------------------------------------------------------------

.. function:: acb_ode_t acb_ode_set (acb_ode_t dest, acb_ode_t src)

	Copies data from *src* to *dest*. If *dest* is *NULL*, a new `acb_ode_t` will be initialised and returned, otherwise only the data will be copied over.

	.. note::
		`acb_ode_set` creates a deep copy of *src* and is therefore rather slow! If *dest* is *NULL*, a pointer to a new `acb_ode_struct` is returned, otherwise *dest* itself will be returned. In either case, the return value should not be ignored but instead be stored in *dest*!

.. function:: slong acb_ode_reduce (acb_ode_t ODE)

	Finds the highest power of *z* that divides every polynomial and uses that to simplify the equation. The return value contains the exponent of z, that the equation was divided by.

.. function:: slong acb_ode_valuation (acb_ode_t ODE)

	Finds the maximum integer *v*, such that :math:`P_{ij} = 0` for all :math:`0 \leq j \leq i-v`, where :math:`P_{ij}` are the coefficients of the defining polynomials of *ODE*.

.. function:: void acb_ode_shift (acb_ode_t ODE_out, acb_ode_t ODE_in, acb_t a, slong bits)

	Transform the origin of *ODE_in* to *a* and store the result in ODE_out.

Input and Output
------------------------------------------------------------------------

.. function:: void acb_ode_dump (acb_ode_t ODE, char* file)

	Dumps the data stored in the `acb_ode_struct` into *file*. If *file == NULL*, the output will be written to *stdout*.
