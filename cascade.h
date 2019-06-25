#ifndef MONODROMY_H
#define MONODROMY_H

#include <acb_poly.h>
#include <acb_mat.h>
#include "acb_ode.h"
#include "flint/profiler.h"

#define TRUE 1

ulong find_power_series (acb_ode_t ODE, acb_t in, slong bits);

void analytic_continuation (acb_t res, acb_ode_t ODE, acb_srcptr path, slong len, slong prec, int output_series);

void find_monodromy_matrix (acb_mat_t monodromy, acb_ode_t ODE, acb_ptr path, slong len, slong bits);

int checkODE (acb_poly_t *polys, acb_ode_t ODE, acb_t z, slong bits);

void entry_point (const ulong n, slong prec, double z, const char *file);

void acb_ode_dump(acb_ode_t ODE);

#endif
