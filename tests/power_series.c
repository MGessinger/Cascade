#include "cascade.h"
#include <acb_poly.h>

int main ()
{
	acb_ode_t ODE = acb_ode_legendre(4);
	acb_poly_t pol;
	acb_poly_init(pol);

	acb_poly_set_coeff_si(pol,0,3);
	find_power_series(pol,ODE,10,1024);
	acb_ode_clear(ODE);

	int return_value = 0;
	acb_t z;
	acb_init(z);
	acb_set_si(z, 3);
	if (!acb_contains(acb_poly_get_coeff_ptr(pol, 0), z))
		return_value = 1;
	if (!acb_contains_zero(acb_poly_get_coeff_ptr(pol, 1)))
		return_value = 1;
	acb_set_si(z, -30);
	if (!acb_contains(acb_poly_get_coeff_ptr(pol, 2), z))
		return_value = 1;
	if (!acb_contains_zero(acb_poly_get_coeff_ptr(pol, 3)))
		return_value = 1;
	acb_set_si(z, 35);
	if (!acb_contains(acb_poly_get_coeff_ptr(pol, 4), z))
		return_value = 1;

	acb_clear(z);
	acb_poly_clear(pol);
	flint_cleanup();
	return return_value;
}
