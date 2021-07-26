#include "acb_ode.h"
#include <acb.h>
#include <flint/flint.h>

int main ()
{
	acb_ode_t ode;
	acb_ode_legendre(ode, 5);
	acb_t exp;
	acb_init(exp);

	int return_value = 0;
	acb_set_si(exp, 30);
	if (!acb_contains(diff_eq_coeff(ode, 0, 0), exp))
		return_value = 1;
	else if (!acb_contains_zero(diff_eq_coeff(ode, 1, 0)))
		return_value = 1;
	acb_set_si(exp, -1);
	if (!acb_contains(diff_eq_coeff(ode, 2, 2), exp))
		return_value = -1;

	acb_clear(exp);
	acb_ode_clear(ode);
	flint_cleanup();
	return return_value;
}
