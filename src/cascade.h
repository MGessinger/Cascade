#ifndef CASCADE_H_
#define CASCADE_H_

#include <acb_poly.h>
#include <acb_mat.h>
#include "acb_ode.h"

/* Compute a power series solution to ODE */
slong find_power_series (acb_poly_t res, acb_ode_t ODE, slong num_of_coeffs, slong bits);
slong truncation_order (arb_t eta, arb_t alpha, slong bits);

/* Compute analytic continuation and monodromy */
void  analytic_continuation (acb_poly_t res, acb_ode_t ODE, acb_srcptr path,
		slong len, slong num_of_coeffs, slong bits);
void  find_monodromy_matrix (acb_mat_t mono, acb_ode_t ODE, acb_t z0, slong bits);
void  radius_of_convergence (arb_t rad_of_conv, acb_ode_t ODE, slong n, slong bits);

/* Additional functions */
void  acb_poly_graeffe_transform (acb_ptr dest, acb_srcptr src, slong len, slong bits);

#endif /* CASCADE_H_ */
