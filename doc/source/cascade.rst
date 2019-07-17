.. _Cascade:

**cascade.h** - Solving complex differential equations
======================================================

Cascade implements complex differential equations in a variable of type :type:`acb_ode_struct`. It provides functionality for finding power series solutions around the origin, but is also capable of performing analytic continuation along a piecewise linear path to find power series expansions anywhere in the complex plane - provided the solution is analytic in a small neighborhood of the origin.

Functions
------------------

.. function:: ulong find_power_series (acb_ode_t ODE, acb_t in, slong bits)

    Finds a solution to *ODE*, by recursively computing the coefficients of the power series expansion. The iteration stops when the summands, given by *a_k\*in^k*, are small enough.  This function does not perform any tests for convergence beforehnd, but it aborts if too many coefficients are computed.

.. function:: void analytic_continuation (acb_t res, acb_ode_t ODE, acb_srcptr path, slong len, slong prec, int output_series)

    Perform analytic continuation along *path*, which stores the corners of a polygon with *len* corners. The function computes a power series at every corner, which is then shifted to the next corner. After the end has been reached, the ODE is shifted back to the point of origin, but the power series is not changed. This allows to perform this function multiple times with varying initial conditions.

.. function:: void find_monodromy_matrix (acb_mat_struct monodromy, acb_ode_t ODE, acb_struct z0, slong bits)

    Compute the monodromy matrix of *ODE* and store it in *monodromy*. This is implemented by performing `analytic_continuation` for multiple different initial conditions.

.. function:: int checkODE (acb_poly_t *polys, acb_ode_t ODE, acb_t z, slong bits)

    Plugs the solution stored in *ODE* into the actual equation and checks if the result actually vanishes. If not, *ODE* is dumped by `acb_ode_dump`.
