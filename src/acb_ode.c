#include "acb_ode.h"

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
	if (  (degree(ODE_out) != degree(ODE_in))
	   || (order(ODE_out) != order(ODE_in))   )
	{
		acb_ode_clear(ODE_out);
		acb_ode_init_blank(ODE_out, order(ODE_in), degree(ODE_in));
	}
	_acb_vec_set(ODE_out->polys, ODE_in->polys, ODE_in->alloc);
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
	return reduced;
}

slong acb_ode_valuation (acb_ode_t ODE)
{
	slong val = 0;
	for (int i = 0; i <= order(ODE); i++)
	{
		slong v = i;
		while (acb_is_zero(acb_ode_coeff(ODE, i, i-v)))
			v--;

		if (v > val)
			val = v;
	}
	return val;
}
