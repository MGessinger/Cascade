#include "acb_ode.h"

#define UNDEFINED -0xFFFF

/* Static function */

static inline int max_degree (acb_poly_t *polys, slong order)
{
	slong deg, poly_max_degree = 0;
	while (polys[order] == NULL || acb_poly_is_zero(polys[order]))
		order--;
	if (order <= 0)
		return -1;

	for (slong i = 0; i <= order; i++)
	{
		if (polys[i] == NULL)
			continue;
		deg = acb_poly_degree(polys[i]);
		if (poly_max_degree < deg)
			poly_max_degree = deg;
	}
	return poly_max_degree;
}

/* Setup and memory management*/

void acb_ode_init_blank (acb_ode_t ODE, slong degree, slong order)
{
	ODE->order = order;
	ODE->degree = degree;
	ODE->polys = NULL;
	ODE->alloc = 0;
	ODE->valuation = UNDEFINED;

	if (degree < 0 || order <= 0)
		return;

	ODE->alloc = (order + 1) * (degree + 1);
	ODE->polys = _acb_vec_init((order(ODE)+1)*(degree(ODE)+1));
}

void acb_ode_init (acb_ode_t ODE, acb_poly_t *polys, slong order)
{
	slong degree = max_degree(polys, order);

	acb_ode_init_blank (ODE, degree, order);
	for (slong i = 0; i <= order(ODE); i++)
	{
		if (polys[i] == NULL)
			continue;
		for (slong j = 0; j < acb_poly_length(polys[i]); j++)
		{
			acb_poly_get_coeff_acb(acb_ode_coeff(ODE, i, j), polys[i], j);
		}
	}
}

void acb_ode_clear (acb_ode_t ODE)
{
	if (ODE->alloc <= 0)
		return;

	for (slong i = 0; i < ODE->alloc; i++)
		acb_clear(ODE->polys + i);

	flint_free(ODE->polys);
}

void acb_ode_set (acb_ode_t ODE_out, acb_ode_t ODE_in)
{
	if (ODE_out->alloc < ODE_in->alloc)
	{
		acb_ode_clear(ODE_out);
		acb_ode_init_blank(ODE_out, order(ODE_in), degree(ODE_in));
	}
	_acb_vec_set(ODE_out->polys, ODE_in->polys, ODE_in->alloc);
	degree(ODE_out) = degree(ODE_in);
	order(ODE_out) = order(ODE_in);
	ODE_out->valuation = acb_ode_valuation(ODE_in);
}

void acb_ode_random (acb_ode_t ode, flint_rand_t state, slong prec)
{
	slong degree, order;

	degree = 3 + n_randint(state, 7);
	order = 2 + n_randint(state, degree-2);

	acb_ode_init_blank(ode, degree, order);

	for (slong i = 0; i <= order(ode); i++)
		for (slong j = 0; j <= degree(ode); j++)
			acb_randtest(acb_ode_coeff(ode, i, j), state, prec, 16);
}

/* I/O */

void acb_ode_dump (acb_ode_t ODE, char *file)
{
	/* Dumps the ODE to file. If file is NULL, dump to stdout */
	FILE *out = stdout;
	if (file != NULL)
		out = fopen(file, "w");
	if (out == NULL)
		return;
	flint_fprintf(out, "Order: %w\nDegree: %w\n", order(ODE), degree(ODE));
	for (slong i = 0; i <= order(ODE); i++)
	{
		flint_fprintf(out, "acb_ode_poly(ODE, %w) = ", i);
		for (slong j = 0; j <= degree(ODE); j++)
		{
			acb_fprintd(out, acb_ode_coeff(ODE, i, j), 20);
			flint_fprintf(out, "\t");
		}
		flint_fprintf(out, "\n");
	}
	if (out != stdout)
		fclose(out);
}

/* Transformations */

void acb_ode_shift (acb_ode_t ODE_out, acb_ode_t ODE_in, acb_srcptr a, slong bits)
{
	/* Shifts the origin to a */
	acb_ode_set(ODE_out, ODE_in);
	if (acb_is_zero(a))
		return;
	if (degree(ODE_in) == 0)
		return;
	for (slong j = 0; j <= order(ODE_out); j++)
		_acb_poly_taylor_shift(acb_ode_poly(ODE_out, j), a, degree(ODE_out)+1, bits);
	ODE_out->valuation = UNDEFINED;
}

slong acb_ode_reduce (acb_ode_t ODE)
{
	/* Divides all polynomials by z^n, if they share such a factor */
	slong reduced = 0;
	int all_zero = 1;
	do {
		for (slong i = 0; i < order(ODE); i++)
			all_zero &= acb_is_zero(acb_ode_coeff(ODE, i, reduced));
		reduced++;
	} while (all_zero);
	reduced -= 1; /* From the loop body */
	if (reduced <= 0)
		return 0;

	slong new_deg = degree(ODE)-reduced;
	for (slong i = 0; i<= order(ODE); i++)
		_acb_poly_shift_right(ODE->polys + i*(new_deg+1), ODE->polys + i*(degree(ODE)+1), degree(ODE)+1, reduced);
	degree(ODE) = new_deg;
	ODE->valuation = UNDEFINED;
	return reduced;
}

slong acb_ode_valuation (acb_ode_t ODE)
{
	if (ODE->valuation != UNDEFINED)
		return ODE->valuation;

	slong val = degree(ODE);
	for (int i = 0; i <= order(ODE); i++)
	{
		slong v = 0;
		while (acb_is_zero(acb_ode_coeff(ODE, i, v)))
		{
			v++;
			if (v > degree(ODE))
				break;
		}

		v = v - i;
		if (v < val)
			val = v;
	}
	ODE->valuation = val;
	return val;
}

/* Differential Action */

void acb_ode_apply (acb_poly_t out, acb_ode_t ODE, acb_poly_t in, slong prec)
{
	acb_poly_t deriv, acc;
	acb_poly_init(deriv);
	acb_poly_init(acc);

	acb_poly_set(deriv, in);
	acb_poly_zero(out);
	for (slong i = 0; i <= order(ODE); i++)
	{
		acb_poly_zero(acc);
		for (slong j = 0; j <= degree(ODE); j++)
			acb_poly_set_coeff_acb(acc, j, acb_ode_coeff(ODE, i, j));

		acb_poly_mul(acc, acc, deriv, prec);
		acb_poly_add(out, out, acc, prec);

		acb_poly_derivative(deriv, deriv, prec);
	}
	acb_poly_clear(deriv);
	acb_poly_clear(acc);
}

int acb_ode_solves (acb_ode_t ODE, acb_poly_t res, slong deg, slong prec)
{
	int solved = 1;
	acb_t coeff;
	acb_poly_t out;

	acb_init(coeff);
	acb_poly_init(out);
	acb_poly_fit_length(out, deg);
	acb_ode_apply(out, ODE, res, prec);
	for (slong i = 0; i < deg; i++)
	{
		acb_poly_get_coeff_acb(coeff, out, i);
		if (!acb_is_finite(coeff))
		{
			flint_printf("Coefficient %w/%w is infinite.\n", i, deg);
			solved = 0;
			break;
		}
		else if (!acb_contains_zero(coeff))
		{
			flint_printf("Coefficient %w/%w is non-zero.\n", i, deg);
			acb_printd(coeff, 20);
			flint_printf("\nPrecision %w.\n", prec);
			solved = 0;
			break;
		}
	}

	acb_poly_clear(out);
	acb_clear(coeff);
	return solved;
}
