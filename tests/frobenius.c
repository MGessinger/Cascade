#include "cascade.h"

int main ()
{
	int return_value = EXIT_SUCCESS;
	slong prec, n, s_rho;

	flint_rand_t state;
	acb_ode_t ODE;
	acb_t rho;
	acb_ode_solution_t sol;

	/* Initialization */
	flint_randinit(state);

	acb_init(rho);

	for (slong iter = 0; iter < 100; iter++)
	{
		prec = 30 + n_randint(state, 128);
		n = 2 + n_randint(state, 30);

		/* Setup */
		acb_ode_random(ODE, state, prec);
		for (slong i = 0; i <= order(ODE); i++)
			for (slong j = 0; j <= i; j++)
				acb_zero(acb_ode_coeff(ODE, i, j));
		acb_one(acb_ode_coeff(ODE, order(ODE), order(ODE)));
		acb_one(acb_ode_coeff(ODE, order(ODE)-1, order(ODE)-1));

		s_rho = order(ODE) - 2;
		acb_set_si(rho, s_rho);
		acb_ode_solution_init(sol, rho, 2, 0);

		acb_ode_solve_frobenius(sol, ODE, n, prec);
		acb_poly_shift_left(sol->gens, sol->gens, s_rho);
		int solved = acb_ode_solves(ODE, sol->gens, n, prec);

		acb_ode_solution_clear(sol);
		acb_ode_clear(ODE);

		if (!solved)
		{
			return_value = EXIT_FAILURE;
			break;
		}
	}
	acb_clear(rho);
	flint_randclear(state);
	flint_cleanup();
	return return_value;
}
