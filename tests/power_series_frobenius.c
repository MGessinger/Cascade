#include "cascade.h"

int main ()
{
	int return_value = EXIT_SUCCESS;
	slong prec, degree, order;

	flint_rand_t state;
	acb_ode_t ODE;
	acb_t rho;
	acb_ode_solution_t sol;

	/* Initialization */
	flint_randinit(state);

	acb_init(rho);

	for (slong iter = 0; iter < 100; iter++)
	{
		prec = 2 + n_randint(state, 128);
		degree = 2 + n_randint(state, 8);
		order = n_randint(state, degree);
		if (order <= 0)
			order = 2;

		acb_ode_init_blank(ODE, degree, order);
		acb_ode_solution_init(sol, rho, 2, 0);

		/* Setup */
		for (slong i = 0; i <= order(ODE); i++)
			for (slong j = i+1; j <= degree(ODE); j++)
				acb_randtest(acb_ode_coeff(ODE, i, j), state, prec, 16);

		acb_one(acb_ode_coeff(ODE, order, order));
		acb_set_si(acb_ode_coeff(ODE, order-1, order-1), order-1);

		acb_ode_solve_frobenius(sol, ODE, 32, prec);

		int solved = acb_ode_solves(ODE, sol->gens + 0, 32, prec);
		if (!solved)
			return_value = EXIT_FAILURE;

		acb_ode_clear(ODE);
		acb_ode_solution_clear(sol);
	}
	acb_clear(rho);
	flint_randclear(state);
	flint_cleanup();
	return EXIT_SUCCESS;
}
