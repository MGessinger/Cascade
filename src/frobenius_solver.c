#include "acb_ode.h"
#include "cascade.h"

void indicial_polynomial (acb_poly_t result, acb_ode_t ODE, slong nu, slong shift, slong prec)
{
	/* Compute f_ν(ρ-ς) [for definition of f, see Frobenius Paper, equation 3] */
	acb_poly_zero(result);
	nu += acb_ode_valuation(ODE);
	if (nu > degree(ODE))
		return;

	acb_t temp1, temp2;

	acb_init(temp1);
	acb_init(temp2);

	acb_poly_fit_length(result, order(ODE) + 1);

	slong lambda = clamp(degree(ODE) - nu, 0, order(ODE));
	for (; lambda >= 0; lambda--)
	{
		for (slong j = acb_poly_length(result); j > 0; j--)
		{
			acb_poly_get_coeff_acb(temp1, result, j);
			acb_poly_get_coeff_acb(temp2, result, j - 1);
			acb_mul_si(temp1, temp1, shift - lambda, prec);
			acb_add(temp1, temp2, temp1, prec);
			acb_poly_set_coeff_acb(result, j, temp1);
		}
		acb_poly_get_coeff_acb(temp1, result, 0);
		acb_mul_si(temp1, temp1, shift - lambda, prec);
		if (lambda + nu >= 0)
			acb_add(temp1, temp1, acb_ode_coeff(ODE, lambda, lambda + nu), prec);
		acb_poly_set_coeff_acb(result, 0, temp1);
	}

	acb_clear(temp1);
	acb_clear(temp2);
}

void indicial_polynomial_evaluate (acb_t result, acb_ode_t ODE, slong nu, acb_t rho, slong shift, slong prec)
{
	nu += acb_ode_valuation(ODE);
	if (nu > degree(ODE))
	{
		acb_zero(result);
		return;
	}

	acb_t temp1, out;

	acb_init(temp1);
	acb_init(out);

	slong lambda = clamp(degree(ODE)-nu, 0, order(ODE));
	for (; lambda >= 0; lambda--)
	{
		acb_add_si(temp1, rho, shift - lambda, prec);
		acb_mul(out, out, temp1, prec);

		if (lambda + nu < 0)
			continue;

		acb_add(out, out, acb_ode_coeff(ODE, lambda, lambda + nu), prec);
	}
	acb_set(result, out);

	acb_clear(temp1);
	acb_clear(out);
}

void _acb_ode_solve_frobenius (acb_poly_t res, acb_ode_t ODE, acb_ode_solution_t rhs, slong sol_degree, slong prec)
{
	acb_t g_new, indicial, g_i;
	acb_init(g_new);
	acb_init(indicial);
	acb_init(g_i);

	acb_t rho;
	acb_init(rho);
	acb_set(rho, rhs->rho);

	acb_poly_fit_length(res, sol_degree + 1);
	if (acb_poly_is_zero(rhs->gens))
		acb_poly_one(res);
	else
	{
		slong v = acb_ode_valuation(ODE);
		acb_sub_si(rho, rho, v, prec);

		acb_poly_get_coeff_acb(g_new, rhs->gens, 0);
		indicial_polynomial_evaluate(indicial, ODE, 0, rho, 0, prec);
		acb_div(g_new, g_new, indicial, prec);

		acb_poly_set_coeff_acb(res, 0, g_new);
	}

	for (slong nu = 1; nu <= sol_degree; nu++)
	{
		acb_poly_get_coeff_acb(g_new, rhs->gens, nu);

		slong i = clamp(nu, 1, degree(ODE));
		indicial_polynomial_evaluate(indicial, ODE, i, rho, nu - i, prec);
		do
		{
			acb_poly_get_coeff_acb(g_i, res, nu - i);
			acb_mul(indicial, indicial, g_i, prec);
			acb_sub(g_new, g_new, indicial, prec);

			i--;
			indicial_polynomial_evaluate(indicial, ODE, i, rho, nu - i, prec);
		} while (i > 0);
		acb_div(g_new, g_new, indicial, prec);
		acb_poly_set_coeff_acb(res, nu, g_new);

	}
	acb_clear(g_new);
	acb_clear(indicial);
	acb_clear(g_i);
	acb_clear(rho);
}

void acb_ode_solve_frobenius (acb_ode_solution_t sol, acb_ode_t ODE, slong sol_degree, slong prec)
{
	if (sol->M == 1)
	{
		_acb_ode_solve_frobenius(sol->gens, ODE, sol, sol_degree, prec);
		return;
	}

	acb_t temp;
	acb_poly_t indicial;
	acb_poly_t g_new;
	acb_poly_struct *g_rho;

	g_rho = flint_malloc( degree(ODE) * sizeof(acb_poly_struct) );
	if (g_rho == NULL)
		return;

	acb_init(temp);
	acb_poly_init(indicial);
	acb_poly_init(g_new);

	for (slong i = 0; i < degree(ODE); i++)
		acb_poly_init(g_rho + i);
	acb_poly_one(g_rho);

	for (slong i = 0; i < sol->M; i++)
	{
		acb_poly_zero(sol->gens + i);
		acb_poly_fit_length(sol->gens + i, sol_degree + 1);
	}
	acb_poly_one(sol->gens);

	for (slong nu = 1; nu <= sol_degree; nu++)
	{
		/* Compute the new coefficient (as a function of rho) */
		slong i = clamp(nu, 1, degree(ODE));
		indicial_polynomial(indicial, ODE, i, nu - i, prec);

		do
		{
			acb_poly_mul(indicial, indicial, g_rho + (i - 1), prec);
			acb_poly_sub(g_new, g_new, indicial, prec);

			i--;
			indicial_polynomial(indicial, ODE, i, nu - i, prec);
		} while (i > 0);

		/* Rescale the indicial polynomial, to keep coefficients small */
		indicial_polynomial_evaluate(temp, ODE, 0, sol->rho, nu, prec);
		if (!acb_contains_zero(temp))
		{
			acb_poly_scalar_div(indicial, indicial, temp, prec);
			acb_poly_scalar_div(g_new, g_new, temp, prec);
		}

		/* Multiply all relevant g_nu(rho) by f(rho + nu) */
		int all_zero = acb_poly_is_zero(g_new);
		for (i = degree(ODE) - 1; i > 0; i--)
		{
			acb_poly_mul(g_rho + i, g_rho + (i - 1), indicial, prec);
			all_zero &= acb_poly_is_zero(g_rho + i);
		}
		acb_poly_set(g_rho, g_new);

		/* Update the G^(i) */
		_acb_ode_solution_update(sol, indicial, prec);

		/* Extend the g^(k)_nu */
		_acb_ode_solution_extend(sol, nu, g_new, prec);

		if (all_zero)
			break;
	}

	_acb_ode_solution_normalize(sol, prec);

	acb_poly_clear(g_new);
	acb_poly_clear(indicial);
	for (slong i = 0; i < degree(ODE); i++)
		acb_poly_clear(g_rho + i);
	flint_free(g_rho);
	acb_clear(temp);
}
