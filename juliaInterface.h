#ifndef JULIAINTERFACE_H_
#define JULIAINTERFACE_H_

#include "acb_ode.h"
#include "cascade.h"

acb_ode_t acb_ode_setup_blank(slong degree, slong order);

void acb_set_initial(acb_ode_t ODE, acb_poly_struct init);

acb_poly_struct find_power_series_julia(acb_ode_t ODE, acb_struct in, slong bits);

#endif
