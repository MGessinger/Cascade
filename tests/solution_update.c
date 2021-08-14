#include "cascade.h"

int main ()
{
	/* Init */
	int return_value = EXIT_SUCCESS;
	acb_ode_solution_t sol;
	acb_poly_t f, g;
	acb_poly_struct H[10];
	acb_t rho, val;
	slong p, n, prec;

	flint_rand_t state;
	flint_randinit(state);

	acb_init(rho);
	acb_init(val);
	acb_poly_init(f);
	acb_poly_init(g);

	for (slong i = 0; i < 10; i++)
		acb_poly_init(H + i);

	for (slong iter = 0; iter < 100; iter++)
	{
		n = n_randint(state, 10);
		prec = 30 + n_randint(state, 128);

		acb_randtest(rho, state, prec, 16);
		acb_ode_solution_init(sol, rho, n_randint(state, 10), 0);

		/* Setup */
		acb_poly_randtest(f, state, 20, prec, 16);

		for (slong i = 0; i < n; i++) {
			acb_poly_zero(H + i);
			acb_poly_randtest(g, state, 20, prec, 16);
			acb_poly_mul(H + i, f, g, prec);
			_acb_ode_solution_extend(sol, i, g, prec);
		}

		_acb_ode_solution_update(sol, f, prec);

		/* Check */
		for (slong j = 0; j < n; j++)
		{
			for (slong i = 0; i < sol->multiplicity; i++)
			{
				acb_poly_evaluate(rho, H + j, sol->rho, prec);
				acb_poly_get_coeff_acb(val, sol->gens + i, j);
				if (!acb_overlaps(rho, val)) {
					return_value = EXIT_FAILURE;
					break;
				}
				acb_poly_derivative(H + j, H + j, prec);
			}
		}

		acb_ode_solution_clear(sol);
	}

	for (slong i = 0; i < 10; i++)
		acb_poly_clear(H + i);

	acb_clear(rho);
	acb_clear(val);
	acb_poly_clear(f);
	acb_poly_clear(g);
	flint_randclear(state);
	flint_cleanup();
	return return_value;
}
