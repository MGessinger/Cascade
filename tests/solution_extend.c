#include "acb_ode.h"

int main ()
{
	/* Init */
	int return_value = EXIT_SUCCESS;
	slong prec, n;
	acb_ode_solution_t sol;
	acb_poly_t f, g;
	acb_t rho, val;

	flint_rand_t state;
	flint_randinit(state);

	acb_init(rho);
	acb_init(val);
	acb_poly_init(f);
	acb_poly_init(g);

	for (slong iter = 0; iter < 100; iter++)
	{
		prec = 30 + n_randint(state, 128);
		n = 2 + n_randint(state, 10);

		acb_randtest(rho, state, prec, 16);
		acb_ode_solution_init(sol, rho, 1 + n_randint(state, 5), 0);

		for (slong i = 0; i < n; i++)
		{
			acb_poly_randtest(f, state, 20, prec, 16);
			acb_poly_set(g, f);
			_acb_ode_solution_extend(sol, i, g, prec);
			for (slong j = 0; j < sol->multiplicity; j++)
			{
				acb_poly_evaluate(rho, f, sol->rho, prec);
				acb_poly_get_coeff_acb(val, sol->gens + j, i);
				if (!acb_equal(rho, val))
					return_value = EXIT_FAILURE;
				acb_poly_derivative(f, f, prec);
			}
		}

		acb_ode_solution_clear(sol);
	}
	acb_clear(rho);
	acb_clear(val);
	acb_poly_clear(f);
	acb_poly_clear(g);
	flint_randclear(state);
	flint_cleanup();
	return return_value;
}
