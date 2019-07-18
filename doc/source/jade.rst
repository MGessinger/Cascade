.. _Nemo: http://nemocas.org
.. _Jade:

**Jade** - the Julia interface to Arbitrary precision Differential Equations
====================================================================================

Jade provides an interactive interface to Cascade, which can be used from julia's REPL. It uses Nemo_ to store variables of type :type:`acb` in julia, which will then be translated to a C-style struct whenever needed. All memory management will be performed automatically, but you can manually create or delete the C-struct if you need to.

Types
--------------------

.. type:: acb_ode

    Contains an array of type acb_poly, the order of the differential equation and a pointer of type *nothing*. The latter is used to store a pointer to an `acb_ode_t` created by Cascade.

.. warning::
    The pointer value should never be accessed directly by the user! Pointers are created and free'd by Cascade. Any incorrect value can (and almost certainly will) lead to corruption and a Segmentation Fault.

Memory Management
--------------------

.. function:: acb_ode(polys::Array{acb_poly,1})

    Creates a struct of type `acb_ode` from *polys*. The order will be determined internally.

.. function:: translateC(ode::acb_ode)

    An `acb_ode_t` is created by Cascade and a pointer to it is stored in *ode*.

.. function:: deleteC(ode::acb_ode)

    The C-style struct stored in *ode* is deallocated by Cascade.

Special Equations
-------------------

Cascade contains a few predefined differential equations. Jade can access those equations by calling to Cascade through one of the following functions. Because Jade however stores the polynomials directly, a polynomial ring has to be specified explicitly. From this ring the precision will be determined.

.. function:: acb_ode_legendre(R::AcbPolyRing, n::Integer)

    Creates an `acb_ode` initialized to Legendre's equation. This is implemented by calling :func:`acb_ode_legendre` through Cascade.

.. function:: acb_ode_bessel(R::AcbPolyRing,nu::acb)

    Creates an `acb_ode` initialized to Bessel's equation. This is implemented by calling :func:`acb_ode_bessel` through Cascade.

.. function:: acb_ode_hypgeom(R::AcbPolyRing, a::acb, b::acb, c::acb)

    Creates an `acb_ode` initialized to Euler's hypergeometric equation. This is implemented by calling :func:`acb_ode_hypgeom` through Cascade.

Solving ODEs
--------------------

.. function:: setPolynomial(ode::acb_ode, index::Integer, polynomial::acb_poly)

    Replace the polynomial at *index* by *polynomial*. If *ode* has already been translated before, the data will be cleared first. *order* will be adjusted accordingly. Remember that Julia counts from 1!

.. function:: setInitialValues(ode::acb_ode,poly::acb_poly)

    Store *poly* in the power series of the C-struct of *ode*. If the struct has not been allocated before, the function will perform that automatically.

.. function:: powerSeries(ode::acb_ode,target::acb)

    Copute a power series solution of *ode*, which converges at *target*, through Cascade. The precision is automatically determined from the polynomials in *ode*.

.. function:: monodromy(ode::acb_ode,z0=0)

    Compute the monodromy matrix of *ode* around *z0* through Cascade, which defaults to Zero.

.. important:: 

    Providing z0 explicitly does not work currently. Most likely there is an issue with the garbage collector of Julia preventing a successfull ccall, which results in a Segmentation Fault.
