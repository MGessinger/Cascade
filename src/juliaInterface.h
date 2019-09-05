#ifndef JULIAINTERFACE_H_
#define JULIAINTERFACE_H_

#include "acb_ode.h"
#include "cascade.h"

short acb_ode_set_poly(acb_ode_t ODE, acb_poly_struct poly, slong index);

void acb_set_initial(acb_ode_t ODE, acb_poly_struct init);

acb_poly_struct find_power_series_julia(acb_ode_t ODE, acb_struct in, slong bits);

#endif
