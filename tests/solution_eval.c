#include "acb_ode.h"

int main ()
{
	int return_value = EXIT_SUCCESS;
	slong prec;
	fmpz_t binom;
	acb_t rho, exp, num;
	acb_t temp1, temp2;
	acb_ode_solution_t sol;
	flint_rand_t state;

	/* Initialization */
	flint_randinit(state);

	fmpz_init(binom);
	acb_init(rho);
	acb_init(exp);
	acb_init(num);
	acb_init(temp1);
	acb_init(temp2);

	for (slong iter = 0; iter < 100; iter++)
	{
		prec = 30 + n_randint(state, 128);

		acb_randtest(rho, state, prec, 10);
		acb_ode_solution_init(sol, rho, 1 + n_randint(state, 10), 0);

		for (slong i = 0; i < sol->M; i++)
			acb_poly_randtest(sol->gens, state, 5, prec, 10);

		acb_randtest(rho, state, prec, 10);
		acb_ode_solution_evaluate(exp, sol, rho, prec);

		acb_zero(num);
		for (slong i = 0; i < sol->M; i++)
		{
			acb_poly_evaluate(temp1, sol->gens + i, rho, prec);

			fmpz_bin_uiui(binom, sol->M - 1, i);
			acb_mul_fmpz(temp2, temp1, binom, prec);

			acb_log(temp1, rho, prec);
			acb_pow_si(temp1, temp1, sol->M - i - 1, prec);
			acb_mul(temp2, temp2, temp1, prec);

			acb_add(num, num, temp2, prec);
		}
		acb_pow(temp1, rho, sol->rho, prec);
		acb_mul(num, num, temp1, prec);

		acb_ode_solution_clear(sol);

		if (!acb_overlaps(exp, num))
		{
			return_value = EXIT_FAILURE;
			break;
		}
	}

	acb_clear(rho);
	acb_clear(num);
	acb_clear(exp);
	acb_clear(temp1);
	acb_clear(temp2);
	fmpz_clear(binom);
	flint_randclear(state);
	flint_cleanup();
	return return_value;
}
