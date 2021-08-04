#include "acb_ode.h"
#include <acb.h>
#include <flint/flint.h>

int main ()
{
	int return_value = 0;
	acb_t a, b, c;
	acb_init(a);
	acb_init(b);
	acb_init(c);

	acb_set_si(a, 1);
	acb_set_si(b, 2);
	acb_set_si(c, 3);

	acb_ode_t ode;
	acb_ode_hypgeom(ode, a, b, c, 1024);

	acb_set_si(a, -2);
	acb_set_si(b, -4);
	if (!acb_contains(acb_ode_coeff(ode, 0, 0), a))
		return_value = 1;
	else if (!acb_contains(acb_ode_coeff(ode, 1, 1), b))
		return_value = 1;
	else if (!acb_contains(acb_ode_coeff(ode, 1, 0), c))
		return_value = 1;
	acb_set_si(a, 1);
	acb_set_si(b, -1);
	if (!acb_contains_zero(acb_ode_coeff(ode, 2, 0)))
		return_value = 1;
	else if (!acb_contains(acb_ode_coeff(ode, 2, 1), a))
		return_value = 1;
	else if (!acb_contains(acb_ode_coeff(ode, 2, 2), b))
		return_value = 1;

	acb_clear(a);
	acb_clear(b);
	acb_clear(c);
	acb_ode_clear(ode);
	flint_cleanup();
	return return_value;
}
