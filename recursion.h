#include <acb_poly.h>
#include <acb_mat.h>
#include "acb_ode.h"

#define TRUE 1

#ifndef RECURSION_H
#define RECURSION_H

void analytic_continuation (acb_struct *res, acb_ode_t ODE_in, acb_srcptr path, slong len, slong prec, int output_series);

ulong find_power_series_regular(acb_t res, acb_ode_t ODE, acb_t a, slong bits);

void find_monodromy_matrix (acb_mat_t monodromy, acb_ode_t ODE, acb_ptr path, slong len, slong digits);

int checkODE (acb_poly_struct **polys, acb_ode_t ODE, slong digits);

void entry_point (ulong n, slong prec, slong z);

#endif
