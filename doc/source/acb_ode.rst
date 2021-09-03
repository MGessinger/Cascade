.. _acb-ode:

Complex Differential Operators
======================================================================

An :type:`acb_ode_t` represents an ordinary, linear differential operator with complex polynomial coefficients and error bounds.
They are implemented as an array of type :type:`acb_struct`.
The order of the differential equation as well as the polynomial degrees are tracked internally and need no further care by the user.

Types and macros
----------------------------------------------------------------------

.. type:: acb_ode_struct

.. type:: acb_ode_t

	This data structure represents a differential operator, whose coefficients are complex-valued polynomials.
	The coefficients of these polynomials should not be indexed directly.
	Instead, use the macros below.

	An `acb_ode_t` is defined as an array of type `acb_ode_struct` of length 1, so it can be passed by reference.

.. macro:: acb_ode_poly (ODE, i)

	A macro returning a pointer to the coefficient belonging to the *i-th* derivative. It is of type `acb_ptr`.
	It is of type `acb_ptr`.

.. macro:: acb_ode_coeff (ODE, i ,j)

	A macro returning a pointer to the *j-th* coefficient of the *i-th* polynomial. This macro is identical to *acb_ode_poly(ODE,i)->(j)*.
	It is of type `acb`.

.. macro:: order (ODE)

	Returns the order of the differential operator, which is one less than the number of defining polynomials.
	That means, there are :math:`order(ODE) + 1` polnomials.

.. macro:: degree (ODE)

	Returns the maximum degree amongst the defining polynomials. This is computed through *acb_poly_degree* and therefore the same restrictions apply in the case of inexact polynomials.

Memory management
----------------------------------------------------------------------

.. function:: void acb_ode_init_blank (acb_ode_t ODE, slong degree, slong order)

	Initializes the variable *ODE* to have space for :math:`order+1` polynomials of length no more than :math:`degree+1`.

.. function:: void acb_ode_init (acb_ode_t ODE, acb_poly_t polys, slong order)

	Initializes the variable *ODE* to contain *polys*.

.. function:: void acb_ode_clear (acb_ode_t ODE)

	Clears the memory allocated by a previous call to `acb_ode_init`.

	.. important::
		Because Arb caches certain values internally, it is recommended to call *flint_cleanup()* at the end of your main program.
		This will clear Arb's cache and lead to a clean output when using *Valgrind*.

Input and Output
----------------------------------------------------------------------

.. function:: void acb_ode_dump (acb_ode_t ODE, char* file)

	Dumps the data stored in the `acb_ode_struct` into *file*. If *file* is *NULL*, the output will be written to *stdout*.

Manipulation
----------------------------------------------------------------------

.. function:: void acb_ode_set (acb_ode_t dest, acb_ode_t src)

	Copies data from *src* to *dest*. The variable *dest* must be initialized by a previous call to `acb_ode_init`, but it need not be of the same size as *src* (in which case it is reallocated).

	.. note::
		`acb_ode_set` creates a deep copy of *src* and is therefore rather slow!

.. function:: slong acb_ode_reduce (acb_ode_t ODE)

	Divides the differential operator defined by *ODE* by the largest common factor of :math:`z^n`, that is shared by all coefficients. Returns the exponent *n*.

.. function:: slong acb_ode_valuation (acb_ode_t ODE)

	Finds the maximum integer *v*, such that :math:`P_{ij} = 0` for all :math:`0 \leq j \leq i-v`, where :math:`P_{ij}` are the coefficients of the defining polynomials of *ODE*.

.. function:: void acb_ode_shift (acb_ode_t ODE_out, acb_ode_t ODE_in, acb_t a, slong bits)

	Transforms the origin of *ODE_in* to *a* and stores the result in ODE_out.
	Both *ODE_in* and *ODE_out* must be initialized.
	Aliasing is permitted.

Differential Action
----------------------------------------------------------------------

.. function:: void acb_ode_apply (acb_poly_t out, acb_ode_t ODE, acb_poly_t in, slong prec)

	Apply the differential operator defined by *ODE* to the polynomial *in*, and store the result in *out*.

.. function:: int acb_ode_solves (acb_ode_t ODE, acb_poly_t res, slong deg, slong prec)

	Test if the polynomial *res* solves the differential equation defined by *ODE* up to degree *deg*.

Special Equations
----------------------------------------------------------------------

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
