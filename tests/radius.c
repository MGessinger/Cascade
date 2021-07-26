#include "cascade.h"
#include <acb.h>
#include <flint/flint.h>

int main ()
{
	int return_value = 0;
	acb_ode_t ode;
	acb_ode_legendre(ode, 4);
	arb_t rad;
	arb_init(rad);
	radius_of_convergence(rad, ode, 20, 1024);

	if (!arb_contains_si(rad, 1))
		return_value = 1;

	arb_clear(rad);
	acb_ode_clear(ode);
	flint_cleanup();
	return 0;
}
