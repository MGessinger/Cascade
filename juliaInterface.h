#ifndef JULIAINTERFACE_H_
#define JULIAINTERFACE_H_

#include "acb_ode.h"
#include "cascade.h"

ulong set_tolerance(ulong new_tol);

acb_ode_t acb_ode_setup_empty(slong degree, slong order);

void acb_set_initial(acb_ode_t ODE, acb_poly_struct init);

acb_poly_struct find_power_series_julia(acb_ode_t ODE, acb_struct in, slong bits);

acb_mat_struct find_monodromy_julia (acb_ode_t ODE, slong bits);

#endif
