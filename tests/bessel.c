#include "acb_ode.h"
#include <acb.h>
#include <flint/flint.h>

int main ()
{
	acb_t nu;
	acb_init(nu);
	acb_set_d(nu, 0.5);
	acb_ode_t ode;
	acb_ode_bessel(ode, nu, 1024);

	acb_set_si(nu, 1);
	int return_value = 0;
	if (!acb_contains(diff_eq_coeff(ode, 2, 2), nu))
		return_value = 1;
	else if (!acb_contains(diff_eq_coeff(ode, 1, 1), nu))
		return_value = 1;
	else if (!acb_contains(diff_eq_coeff(ode, 0, 2), nu))
		return_value = 1;
	acb_set_d(nu, -0.25);
	if (!acb_contains(diff_eq_coeff(ode, 0, 0), nu))
		return_value = 1;

	acb_ode_clear(ode);
	acb_clear(nu);
	flint_cleanup();
	return return_value;
}
