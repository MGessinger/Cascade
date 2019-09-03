.. _Nemo: http://nemocas.org
.. _Jade:

**Jade** - the Julia interface to Arbitrary Precision Differential Equations
====================================================================================

Jade provides an interactive interface to Cascade, which can be used from Julia's REPL. It uses Nemo_ to store variables of type :type:`acb` in Julia, which will then be translated to a C-style struct whenever needed. All memory management will be performed automatically, but you can manually create or delete the C-struct if you need to.

In the following, all functions are listed with their full name. However most functions are exported from Jade and can therefore be called without the *Jade.*-prefix.

Types and Aliases
--------------------

.. type:: jade_ode{P<:PolyElem}

    Contains an array of type *P*, the order of the differential equation and a pointer of type *nothing*. If *P == acb_poly*, then the latter is used to store a pointer to an `acb_ode_t` created by Cascade.

.. warning::

    The pointer value should never be accessed directly by the user! Pointers are created and free'd by Cascade. Any incorrect value can (and almost certainly will) lead to corruption and a Segmentation Fault.

.. type:: arb_ode

	This is an alias for jade_ode{arb_poly}.

.. type:: acb_ode

	This is an alias for jade_ode{acb_poly}.

.. type:: diffEq

	An abstract type, acting as a parent type for all jade_ode{P}-types, i.e. jade_ode{P <:PolyElem} <: diffEq.

Memory Management
--------------------

.. function:: Jade.jade_ode{P}(polys::Vector{P})

    This is the recommended constructor for variables of type :type:`jade_ode`. The polynomials are ordered in such a way that the polynomial corresponding to the n-th derivative has the index n+1 within the array.

.. function:: Jade.translateC(ode::acb_ode)

    An `acb_ode_t` is created by Cascade and a pointer to it is stored in *ode*. This function is not exported because it should not be called manually!

.. function:: Jade.deleteC(ode::acb_ode)

    The C-style struct stored in *ode* is deallocated by Cascade. This function is not exported because it usually doesn't need to be called manually!

Special Equations
-------------------

Cascade contains a few predefined differential equations. Jade can access those equations by calling to Cascade through one of the following functions. Because Jade however stores the polynomials directly, a polynomial ring has to be specified explicitly. From this ring the precision will be determined.

.. function:: Jade.acb_ode_legendre(R::AcbPolyRing, n::Integer)

    Creates an `acb_ode` initialized to Legendre's equation. This is implemented by calling :func:`acb_ode_legendre` through Cascade.

.. function:: Jade.acb_ode_bessel(R::AcbPolyRing,nu::acb)

    Creates an `acb_ode` initialized to Bessel's equation. This is implemented by calling :func:`acb_ode_bessel` through Cascade.

.. function:: Jade.acb_ode_hypgeom(R::AcbPolyRing, a::acb, b::acb, c::acb)

    Creates an `acb_ode` initialized to Euler's hypergeometric equation. This is implemented by calling :func:`acb_ode_hypgeom` through Cascade.

Solving ODEs
--------------------

.. function:: Jade.setPolynomial(ode::acb_ode, index::Integer, polynomial::acb_poly)

    Replace the polynomial at *index* by *polynomial*. If *ode* has already been translated before, the data will be cleared first. *order* will be adjusted accordingly. Remember that Julia counts from 1!

.. function:: Jade.setInitialValues(ode::acb_ode,poly::acb_poly)

    Store *poly* in the power series of the C-struct of *ode*. If the struct has not been allocated before, the function will perform that automatically.

.. function:: Jade.powerSeries(ode::acb_ode,target::acb)

    Copute a power series solution of *ode*, which converges at *target*, through Cascade. The precision is automatically determined from the polynomials in *ode*.

.. function:: Jade.monodromy(ode::acb_ode,z0=0)

    Compute the monodromy matrix of *ode* around *z0* through Cascade. The value of *z0* defaults to zero.

