#include <acb_poly.h>
#include <acb_mat.h>
#include "acb_ode.h"

#define TRUE 1

#ifndef RECURSION_H
#define RECURSION_H

void analytic_continuation (acb_struct *res, acb_ode_t ODE_in, acb_srcptr path, slong len, slong prec, int output_series);

int find_power_series_regular(acb_t res, acb_ode_t ODE, acb_t a, slong bits);

void find_monodromy_matrix (acb_mat_t monodromy, acb_ode_t ODE, acb_ptr path, slong len, slong digits);

void entry_point (slong n, slong prec, slong z);

#endif
