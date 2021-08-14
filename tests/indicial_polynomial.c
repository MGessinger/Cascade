#include "acb_ode.h"
#include "cascade.h"

int main ()
{
	int return_value = EXIT_SUCCESS;
	slong prec, rho;

	acb_t num, exp;
	acb_ode_t ODE;
	acb_poly_t poly, indicial;
	flint_rand_t state;

	/* Initialization */
	flint_randinit(state);

	acb_init(num);
	acb_init(exp);

	acb_poly_init(poly);
	acb_poly_init(indicial);

	for (slong iter = 0; iter < 100; iter++)
	{
		rho = n_randint(state, 10);
		prec = 30 + n_randint(state, 128);

		/* Setup */
		acb_ode_random(ODE, state, prec);
		for (slong i = 1; i <= order(ODE); i++)
			for (slong j = 0; j < i; j++)
				acb_zero(acb_ode_coeff(ODE, i, j));

		acb_poly_zero(poly);
		acb_one(exp);
		acb_poly_set_coeff_acb(poly, rho, exp);

		acb_ode_apply(poly, ODE, poly, prec * 3);

		for (slong i = 0; i <= degree(ODE)+1; i++)
		{
			acb_poly_get_coeff_acb(exp, poly, rho+i);

			/* Directly compute the "indicial coefficient" */
			acb_set_si(num, rho);
			indicial_polynomial_evaluate(num, ODE, i, num, 0, prec);
			if (!acb_overlaps(exp, num))
			{
				return_value = EXIT_FAILURE;
				break;
			}

			/* Compare to indicial polynomial */
			acb_set_si(num, rho);
			indicial_polynomial(indicial, ODE, i, 0, prec);
			acb_poly_evaluate(num, indicial, num, prec);
			if (!acb_overlaps(exp, num))
			{
				return_value = EXIT_FAILURE;
				break;
			}

			/* Compute f(0 + rho) instead of f(rho + 0) */
			acb_zero(num);
			indicial_polynomial_evaluate(num, ODE, i, num, rho, prec);
			if (!acb_overlaps(exp, num))
			{
				return_value = EXIT_FAILURE;
				break;
			}

			/* Again, compare to indicial polynomial */
			acb_zero(num);
			indicial_polynomial(indicial, ODE, i, rho, prec);
			acb_poly_evaluate(num, indicial, num, prec);
			if (!acb_overlaps(exp, num))
			{
				return_value = EXIT_FAILURE;
				break;
			}
		}

		acb_ode_clear(ODE);
	}

	/* Memory Cleanup */
	acb_poly_clear(poly);
	acb_poly_clear(indicial);
	acb_clear(num);
	acb_clear(exp);
	flint_randclear(state);
	flint_cleanup();
	return return_value;
}
