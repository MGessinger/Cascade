#ifndef CASCADE_H_
#define CASCADE_H_

#include <acb_poly.h>
#include <acb_mat.h>
#include "acb_ode.h"
#include "examples.h"

#define TRUE 1

/* Compute (and double check) a power series solution to ODE */
ulong find_power_series (acb_ode_t ODE, acb_t in, slong bits);
int checkODE (acb_poly_t *polys, acb_ode_t ODE, acb_t z, slong bits);

void analytic_continuation (acb_t res, acb_ode_t ODE, acb_srcptr path, slong len, slong prec, int output_series);

void find_monodromy_matrix (acb_mat_t monodromy, acb_ode_t ODE, acb_t z0, slong bits);

void radiusOfConvergence(acb_ode_t ODE, arf_t radOfConv, slong bits);

#endif /* CASCADE_H_ */
