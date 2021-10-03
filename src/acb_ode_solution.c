#include "acb_ode.h"

void acb_ode_solution_init (acb_ode_solution_t sol, acb_t rho, slong mul, slong alpha)
{
	sol->mul = mul;
	sol->M = mul + alpha;

	acb_init(sol->rho);
	acb_set(sol->rho, rho);

	sol->gens = flint_malloc(sol->M * sizeof(acb_poly_struct));
	for (slong i = 0; i < sol->M; i++)
		acb_poly_init(sol->gens + i);
}

void acb_ode_solution_clear (acb_ode_solution_t sol)
{
	acb_clear(sol->rho);
	for (slong i = 0; i < sol->M; i++)
		acb_poly_clear(sol->gens + i);
	flint_free(sol->gens);
}

void acb_ode_solution_evaluate (acb_t out, acb_ode_solution_t sol, acb_t a, slong prec)
{
	acb_t l, p, res;
	slong binom = 1;

	if (acb_is_zero(a))
	{
		acb_indeterminate(out);
		return;
	}

	acb_init(l);
	acb_init(p);
	acb_init(res);

	acb_log(l, a, prec);
	acb_poly_evaluate(res, sol->gens, a, prec);
	for (slong i = 1; i < sol->M; i++)
	{
		acb_mul(res, res, l, prec);

		acb_poly_evaluate(p, sol->gens + i, a, prec);
		binom = (binom * (sol->M - i + 1)) / i;
		acb_mul_si(p, p, binom, prec);

		acb_add(res, res, p, prec);
	}

	acb_pow(p, a, sol->rho, prec);
	acb_mul(out, res, p, prec);

	acb_clear(res);
	acb_clear(l);
	acb_clear(p);
}

void _acb_ode_solution_update (acb_ode_solution_t sol, acb_poly_t f, slong prec)
{
	acb_struct *F;
	slong binom = 1;

	F = flint_malloc(sol->M * sizeof(acb_struct));
	if (F == NULL)
		return;

	for (slong k = 0; k < sol->M; k++)
	{
		acb_init(F + k);
		acb_poly_evaluate(F + k, f, sol->rho, prec);
		acb_poly_derivative(f, f, prec);

		acb_mul_si(F + k, F + k, binom, prec);
		binom = (binom * (sol->M - 1 - k)) / (k + 1);
	}

	for (slong n = sol->M - 1; n >= 0; n--)
	{
		acb_poly_scalar_mul(sol->gens + n, sol->gens + n, F + 0, prec);
		for (slong k = 1; k <= n; k++)
		{
			acb_poly_scalar_mul(f, sol->gens + (n - k), F + k, prec);
			acb_poly_add(sol->gens + n, sol->gens + n, f, prec);

			acb_mul_si(F + k, F + k, n - k, prec);
			acb_div_si(F + k, F + k, n, prec);
		}
		acb_clear(F + n);
	}

	acb_poly_zero(f);
	flint_free(F);
}

void _acb_ode_solution_extend (acb_ode_solution_t sol, slong nu, acb_poly_t g_nu, slong prec)
{
	acb_t temp;
	acb_init(temp);
	for (slong i = 0; i < sol->M; i++)
	{
		acb_poly_evaluate(temp, g_nu, sol->rho, prec);
		acb_poly_set_coeff_acb(sol->gens + i, nu, temp);
		acb_poly_derivative(g_nu, g_nu, prec);
	}
	acb_poly_zero(g_nu);
	acb_clear(temp);
}

void _acb_ode_solution_normalize (acb_ode_solution_t sol, slong prec)
{
	acb_t t;
	acb_init(t);

	slong alpha = sol->M - sol->mul;
	acb_poly_get_coeff_acb(t, sol->gens + alpha, 0);

	for (slong i = 0; i < sol->M; i++)
		acb_poly_scalar_div(sol->gens + i, sol->gens + i, t, prec);

	acb_clear(t);
}
