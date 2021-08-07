#include "cascade.h"

int main ()
{
	int return_value = EXIT_SUCCESS;
	slong p, prec, degree, order;

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
		p = n_randprime(state, 8, 1);
		prec = 2 + n_randint(state, 62);
		degree = 2 + n_randint(state, 8);
		order = n_randint(state, degree);
		if (order <= 0)
			order = 2;

		acb_poly_init2(result, 32);

		acb_ode_init_blank(ODE, degree, order);

		/* Setup */
		for (slong i = 0; i <= order(ODE); i++)
		{
			for (slong j = i+1; j <= degree(ODE); j++)
				acb_randtest(acb_ode_coeff(ODE, i, j), state, prec, 16);
		}
		acb_one(acb_ode_coeff(ODE, order(ODE), order(ODE)));
		acb_set_si(num, order(ODE) - 1);

		acb_poly_one(result);
		_acb_ode_solve_frobenius(result, ODE, num, 32, prec);
		acb_poly_shift_left(result, result, order(ODE) - 1);

		int solved = acb_ode_solves(ODE, result, 30, prec);
		if (!solved)
			return_value = EXIT_FAILURE;

		acb_ode_clear(ODE);
	}
	acb_clear(num);
	acb_poly_clear(result);
	flint_randclear(state);
	flint_cleanup();
	return return_value;
}
