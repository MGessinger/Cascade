.. _Nemo: http://nemocas.org
.. _jade:

**JADE** - the Julia interface to Arbitrary precision Differential Equations
====================================================================================

Jade provides an interactive interface to Cascade, which can be used from julia's REPL. It uses Nemo_ to store variables of type :type:`acb` in julia, which will then be translated to a C-style struct whenever needed. All memory management will be performed automatically, but you can manually create or delete the C-struct if you need to.

Types
------------

.. type:: acb_ode

    Contains an array of type acb_poly, the order of the differential equation and a pointer of type nothing. The latter is used to store a pointer to an `acb_ode_t` created by Cascade.

.. danger::
    The pointer value should never be accessed directly by the user! Pointers are created and free'd by Cascade and any incorrect value can (and almost certainly will) lead to corruption and a Segmentation Fault.

Memory Management
-----------------

.. function:: acb_ode(polys::Array{acb_poly,1})

    Creates a struct of type `acb_ode` from *polys*. The order will be determined internally.

.. function:: translateC(ode::acb_ode)

    An `acb_ode_t` is created by Cascade and a pointer to it is stored in *ode*.

.. function:: deleteC(ode::acb_ode)

    The C-style struct stored in *ode* is deallocated by Cascade.
