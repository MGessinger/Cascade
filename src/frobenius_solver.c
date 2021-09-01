#include "acb_ode.h"
#include "cascade.h"

void indicial_polynomial (acb_poly_t result, acb_ode_t ODE, slong nu, slong shift, slong prec)
{
	/* Compute f_ν(ρ-ς) [for definition of f, see Frobenius Paper, equation 3] */
	acb_poly_zero(result);
	if (nu > degree(ODE))
		return;

	acb_t temp1;
	acb_poly_t horner;

	acb_init(temp1);
	acb_poly_init(horner);

	acb_one(temp1);
	acb_poly_set_coeff_acb(horner, 1, temp1);

	acb_poly_fit_length(result, order(ODE) + 1);

	slong lambda = clamp(degree(ODE) - nu, 0, order(ODE));
	for (; lambda >= 0; lambda--)
	{
		acb_set_si(temp1, shift - lambda);
		acb_poly_set_coeff_acb(horner, 0, temp1);
		acb_poly_mul(result, result, horner, prec);

		acb_poly_get_coeff_acb(temp1, result, 0);
		acb_add(temp1, temp1, acb_ode_coeff(ODE, lambda, lambda + nu), prec);
		acb_poly_set_coeff_acb(result, 0, temp1);
	}

	acb_clear(temp1);
	acb_poly_clear(horner);
}

void indicial_polynomial_evaluate (acb_t result, acb_ode_t ODE, slong nu, acb_t rho, slong shift, slong prec)
{
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
		acb_add(out, out, acb_ode_coeff(ODE, lambda, lambda + nu), prec);
	}
	acb_set(result, out);

	acb_clear(temp1);
	acb_clear(out);
}

void _acb_ode_solve_frobenius (acb_poly_t res, acb_ode_t ODE, acb_t rho, slong sol_degree, slong prec)
{
	acb_t g_new, indicial, g_i;
	acb_init(g_new);
	acb_init(indicial);
	acb_init(g_i);

	acb_poly_fit_length(res, sol_degree + 1);
	acb_poly_one(res);
	for (slong nu = 1; nu <= sol_degree; nu++)
	{
		acb_zero(g_new);

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
}

void acb_ode_solve_frobenius (acb_ode_solution_t sol, acb_ode_t ODE, slong sol_degree, slong prec)
{
	if (sol->multiplicity == 1 && sol->alpha == 0)
	{
		_acb_ode_solve_frobenius(sol->gens, ODE, sol->rho, sol_degree, prec);
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

	for (slong i = 0; i < sol->multiplicity; i++)
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

	acb_poly_clear(g_new);
	acb_poly_clear(indicial);
	for (slong i = 0; i < degree(ODE); i++)
		acb_poly_clear(g_rho + i);
	flint_free(g_rho);
	acb_clear(temp);
}
