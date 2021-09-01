#include "acb_ode.h"

void acb_ode_solution_init (acb_ode_solution_t sol, acb_t rho, slong mul, slong alpha)
{
	sol->multiplicity = mul;
	sol->alpha = alpha;

	acb_init(sol->rho);
	acb_set(sol->rho, rho);

	slong mu = sol->multiplicity + sol->alpha;
	sol->gens = flint_malloc(mu * sizeof(acb_poly_struct));
	for (slong i = 0; i < mu; i++)
		acb_poly_init(sol->gens + i);
}

void acb_ode_solution_clear (acb_ode_solution_t sol)
{
	acb_clear(sol->rho);
	slong mu = sol->multiplicity + sol->alpha;
	for (slong i = 0; i < mu; i++)
		acb_poly_clear(sol->gens + i);
	flint_free(sol->gens);
}

void _acb_ode_solution_update (acb_ode_solution_t sol, acb_poly_t f, slong prec)
{
	acb_struct *F;
	acb_t temp1, temp2;

	slong mu = sol->multiplicity + sol->alpha;
	F = flint_malloc(mu * sizeof(acb_struct));
	if (F == NULL)
		return;

	acb_init(temp1);
	acb_init(temp2);

	acb_one(temp2);

	for (slong k = 0; k < mu; k++)
	{
		acb_init(F + k);
		acb_poly_evaluate(F + k, f, sol->rho, prec);
		acb_poly_derivative(f, f, prec);

		acb_mul(F + k, F + k, temp2, prec);
		acb_set_si(temp1, mu - 1 - k);
		acb_mul(temp2, temp2, temp1, prec);
		acb_set_si(temp1, k + 1);
		acb_div(temp2, temp2, temp1, prec);
	}

	for (slong n = mu - 1; n >= 0; n--)
	{
		acb_poly_scalar_mul(sol->gens + n, sol->gens + n, F + 0, prec);
		acb_set_si(temp2, n);
		for (slong k = 1; k <= n; k++)
		{
			acb_poly_scalar_mul(f, sol->gens + (n - k), F + k, prec);
			acb_poly_add(sol->gens + n, sol->gens + n, f, prec);

			acb_set_si(temp1, n - k);
			acb_mul(F + k, F + k, temp1, prec);
			acb_div(F + k, F + k, temp2, prec);
		}
		acb_clear(F + n);
	}

	acb_poly_zero(f);

	acb_clear(temp1);
	acb_clear(temp2);
	flint_free(F);
}

void _acb_ode_solution_extend (acb_ode_solution_t sol, slong nu, acb_poly_t g_nu, slong prec)
{
	acb_t temp;
	acb_init(temp);
	slong mu = sol->multiplicity + sol->alpha;
	for (slong i = 0; i < mu; i++)
	{
		acb_poly_evaluate(temp, g_nu, sol->rho, prec);
		acb_poly_set_coeff_acb(sol->gens + i, nu, temp);
		acb_poly_derivative(g_nu, g_nu, prec);
	}
	acb_poly_zero(g_nu);
	acb_clear(temp);
}
