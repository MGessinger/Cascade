#include "cascade.h"

int main ()
{
	int return_value = EXIT_SUCCESS;
	slong prec, degree, order, n;

	acb_ode_t ODE;
	acb_t num;
	acb_poly_t result;
	flint_rand_t state;

	/* Initialization */
	flint_randinit(state);

	acb_init(num);
	acb_poly_init(result);

	for (slong iter = 0; iter < 100; iter++)
	{
		prec = 30 + n_randint(state, 128);
		n = 2 + n_randint(state, 30);

		acb_ode_random(ODE, state, prec);

		/* Setup */
		for (slong i = 0; i <= order(ODE); i++)
			for (slong j = 0; j <= i; j++)
				acb_zero(acb_ode_coeff(ODE, i, j));
		acb_one(acb_ode_coeff(ODE, order(ODE), order(ODE)));
		acb_set_si(num, order(ODE) - 1);

		acb_poly_one(result);
		_acb_ode_solve_frobenius(result, ODE, num, n, prec);
		acb_poly_shift_left(result, result, order(ODE) - 1);

		int solved = acb_ode_solves(ODE, result, n-1, prec);

		acb_ode_clear(ODE);

		if (!solved)
		{
			return_value = EXIT_FAILURE;
			break;
		}
	}
	acb_clear(num);
	acb_poly_clear(result);
	flint_randclear(state);
	flint_cleanup();
	return return_value;
}
