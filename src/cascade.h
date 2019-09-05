#ifndef CASCADE_H_
#define CASCADE_H_

#include <acb_poly.h>
#include <acb_mat.h>
#include "acb_ode.h"

#define TRUE 1

/* Compute (and double check) a power series solution to ODE */
slong find_power_series (acb_ode_t ODE, slong numOfCoeffs, slong bits);
slong truncation_order (arb_t eta, arb_t alpha, slong bits);
int   checkODE (acb_poly_t *polys, acb_ode_t ODE, acb_t z, slong bits);

/* Compute analytic continuation and monodromy */
void  analytic_continuation (acb_t res, acb_ode_t ODE, acb_srcptr path,
                            slong len, slong numOfCoeffs, slong bits, int output_series);
void  find_monodromy_matrix (acb_mat_t monodromy, acb_ode_t ODE, acb_t z0, slong bits);
void  radiusOfConvergence(arb_t radOfConv, acb_ode_t ODE, slong n, slong bits);

/* Additional functions */
void  acb_poly_graeffe_transform(acb_ptr dest, acb_srcptr src, slong len, slong bits);

#endif /* CASCADE_H_ */
