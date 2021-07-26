#include "acb_ode.h"

/* Static function */

static short interpret (acb_poly_t *polys, acb_ode_t ODE)
{
	/* Find the characteristic "size" of the ODE */
	slong poly_max_deg = 0, ord = order(ODE);
	while (polys[ord] == NULL || acb_poly_is_zero(polys[ord]))
		ord--;
	if (ord <= 0)
		return -1;
	for (slong i = 0; i <= ord; i++)
	{
		if (polys[i] == NULL)
			continue;
		if (poly_max_deg < acb_poly_degree(polys[i]))
			poly_max_deg = acb_poly_degree(polys[i]);
	}
	order(ODE) = ord;
	degree(ODE) = poly_max_deg;
	return 1;
}

/* Setup and memory management*/

void acb_ode_init_blank (acb_ode_t ODE, slong degree, slong order)
{
	/* Prepare the Differential equation for later use */
	order(ODE) = order;
	degree(ODE) = degree;
	ODE->polys = NULL;
	if (degree < 0 || order <= 0)
		return;
	ODE->polys = _acb_vec_init((order(ODE)+1)*(degree(ODE)+1));
}

void acb_ode_init (acb_ode_t ODE, acb_poly_t *polys, slong order)
{
	/* Create a differential operator defined by *polys* */
	order(ODE) = order;
	if (interpret(polys, ODE) != 1)
		return;

	ODE->polys = NULL;
	if (degree(ODE) < 0 || order(ODE) <= 0)
		return;
	ODE->polys = _acb_vec_init((order(ODE)+1)*(degree(ODE)+1));
	for (slong i = 0; i <= order(ODE); i++)
	{
		if (polys[i] == NULL)
			continue;
		for (slong j = 0; j <= acb_poly_degree(polys[i]); j++)
		{
			acb_poly_get_coeff_acb(diff_eq_coeff(ODE, i, j), polys[i], j);
		}
	}
}

void acb_ode_clear (acb_ode_t ODE)
{
	if (degree(ODE) < 0 || order(ODE) <= 0)
		return;
	/* Free memory allocated for ODE*/
	_acb_vec_clear(diff_eq_poly(ODE, 0), (order(ODE)+1)*(degree(ODE)+1));
}

void acb_ode_set (acb_ode_t ODE_out, acb_ode_t ODE_in)
{
	if (  (degree(ODE_out) != degree(ODE_in))
	   || (order(ODE_out) != order(ODE_in))   )
	{
		acb_ode_clear(ODE_out);
		acb_ode_init_blank(ODE_out, order(ODE_in), degree(ODE_in));
	}
	_acb_vec_set(ODE_out->polys, ODE_in->polys, (order(ODE_in)+1)*(degree(ODE_in)+1));
}

void acb_ode_set_poly (acb_ode_t ODE, acb_poly_t poly, slong index)
{
	if (index < 0 || index > order(ODE)+1)
		return;
	slong deg = acb_poly_degree(poly);
	if (deg > degree(ODE))
		deg = degree(ODE);    /* Automatically truncates */
	_acb_vec_set(diff_eq_poly(ODE, index), acb_poly_get_coeff_ptr(poly, 0), deg+1);
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
		flint_fprintf(out, "diff_eq_poly(ODE, %w) = ", i);
		for (slong j = 0; j <= degree(ODE); j++)
		{
			acb_fprintd(out, diff_eq_coeff(ODE, i, j), 20);
			flint_fprintf(out, "\t");
		}
		flint_fprintf(out, "\n");
	}
	if (out != stdout)
		fclose(out);
	return;
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
		_acb_poly_taylor_shift(diff_eq_poly(ODE_out, j), a, degree(ODE_out)+1, bits);
}

slong acb_ode_reduce (acb_ode_t ODE)
{
	/* Divides all polynomials by z^n, if they share such a factor */
	slong reduced = 0;
	int all_zero = 1;
	do {
		for (slong i = 0; i < order(ODE); i++)
			all_zero &= acb_is_zero(diff_eq_coeff(ODE, i, reduced));
		reduced++;
	} while (all_zero);
	reduced -= 1; /* From the loop body */
	if (reduced <= 0)
		return 0;
	acb_ode_t newODE;
	acb_ode_init_blank(newODE, degree(ODE)-reduced, order(ODE));
	for (slong i = 0; i<= order(ODE); i++)
		_acb_poly_shift_right(diff_eq_poly(newODE, i), diff_eq_poly(ODE, i), degree(ODE)+1, reduced);
	free(ODE->polys);
	ODE->polys = newODE->polys;
	degree(ODE) -= reduced;
	return reduced;
}

slong acb_ode_valuation (acb_ode_t ODE)
{
	slong val = 0;
	for (int i = 0; i <= order(ODE); i++)
	{
		slong v = i;
		while (acb_is_zero(diff_eq_coeff(ODE, i, i-v)))
			v--;

		if (v > val)
			val = v;
	}
	return val;
}
