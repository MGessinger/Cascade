.. _Nemo: http://nemocas.org
.. _Jade:

**Jade** - the Julia interface to Arbitrary Precision Differential Equations
====================================================================================

Jade provides a wrapper for Cascade, which can be used from Julia's REPL. It uses Nemo_ to store variables of type :type:`acb` in Julia, which will then be translated to a C-style struct whenever needed. All memory management will be performed automatically, but the C-struct can be created and destroyed if absolutely necessary, however this is strongly discouraged.

In the following, all functions are listed with their full name. However most functions are exported from Jade and can therefore be called without the *Jade.*-prefix.

Types and Aliases
--------------------

.. type:: jade_ode{P<:Poly_elem}

	Contains an array of type *P*, the order of the differential equation and a pointer of type *Nothing*. If *P == acb_poly*, then the latter is used to store a pointer to an `acb_ode_t` created by Cascade.

.. warning::

	The pointer value should never be accessed directly by the user! Pointers are created and free'd by Cascade. Any incorrect value can (and almost certainly will) lead to corruption and a Segmentation Fault.

.. type:: arb_ode

	This is an alias for jade_ode{arb_poly}.

.. type:: acb_ode

	This is an alias for jade_ode{acb_poly}.

.. type:: diff_op

	An abstract type, acting as a parent type for all jade_ode{P}-types, i.e. jade_ode{P <:Poly_elem} <: diff_op.

Memory Management
--------------------

.. function:: Jade.jade_ode{P}(polys::Vector{P})

	This is the recommended constructor for variables of type :type:`jade_ode`. The polynomials are ordered in such a way that the polynomial corresponding to the n-th derivative has the index n+1 within the array.

.. function:: Jade.translate_c(ode::acb_ode)

	An `acb_ode_t` is created by Cascade and a pointer to it is stored in *ode*. This function is not exported because it should not be called manually!

	Currently there is no method for parsing other diff_op types other than acb_ode to an acb_ode_t structure.

.. function:: Jade.delete_c(ode::acb_ode)

	The C-style struct stored in *ode* is deallocated by Cascade. This function is not exported because it usually doesn't need to be called manually!

Special Equations
-------------------

Cascade contains a few predefined differential equations. Jade can access those equations by calling to Cascade through one of the following functions. Because Jade however stores the polynomials directly, a polynomial ring has to be specified explicitly. From this ring the precision will be determined.

.. function:: Jade.acb_ode_legendre(R::Acb_poly_ring, n::Integer)

	Creates an `acb_ode` initialized to Legendre's equation. This automatically calls :func:`acb_ode_legendre` through Cascade.

.. function:: Jade.acb_ode_bessel(R::Acb_poly_ring,nu::acb)

	Creates an `acb_ode` initialized to Bessel's equation. This automatically calls :func:`acb_ode_bessel` through Cascade.

.. function:: Jade.acb_ode_hypgeom(R::Acb_poly_ring, a::acb, b::acb, c::acb)

	Creates an `acb_ode` initialized to Euler's hypergeometric equation. This automatically calls :func:`acb_ode_hypgeom` through Cascade.

Arithmetic
--------------

Simple arithmetic has been implemented to work with varibles of type *diff_op*. Addition and subtraction as well as multiplication and division by a scalar are implemented as method extensions to the +,-,*,/ functions of the Base package.
Due to limitations in Nemo's type promotion (e.g. from *fmpz* to *acb*), multiplication might fail if the types are incomptible. However all Julia base types are fully supported.

Because a variable of type :type:`jade_ode` represents a differential operator, they are callable on any polynomial, that supports Nemo's *derivative* function.

Solving ODEs
--------------------

.. function:: Jade.power_series(ode::acb_ode,p::acb_poly,n::Integer)

	Copute a power series solution of *ode*, which converges when evaluated at *target*, through Cascade. The precision is automatically determined from the polynomials in *ode*. The inital values are taken from *p*, which also stores the result.

.. function:: Jade.monodromy(ode::acb_ode)

	Compute the monodromy matrix of *ode* about the origin through Cascade.
