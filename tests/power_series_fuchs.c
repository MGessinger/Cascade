#include "cascade.h"
#include <acb_poly.h>

int main ()
{
	int return_value = EXIT_SUCCESS;

	slong prec, n, deg;

	flint_rand_t state;
	
	acb_poly_t result;
	acb_ode_t ODE;

	flint_randinit(state);
	acb_poly_init(result);

	for (slong iter = 0; iter < 100; iter++)
	{
		prec = 2 + n_randint(state, 126);

		acb_ode_random(ODE, state, prec);
		if (acb_contains_zero(acb_ode_coeff(ODE, order(ODE), order(ODE))))
			acb_one(acb_ode_coeff(ODE, order(ODE), order(ODE)));

		n = order(ODE) + n_randint(state, 32);
		deg = clamp(n - order(ODE), 0, n);

		acb_poly_randtest(result, state, order(ODE) + 1, prec, 8);

		acb_ode_solve_fuchs(result, ODE, n, prec);

		int solved = acb_ode_solves(ODE, result, deg, prec);

		acb_ode_clear(ODE);

		if (!solved)
		{
			return_value = EXIT_FAILURE;
			break;
		}
	}

	acb_poly_clear(result);
	flint_randclear(state);
	return return_value;
}
