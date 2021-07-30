#include "cascade.h"

static inline slong clamp (slong in, slong min, slong max)
{
	if (in < min)
		return min;
	else if (in > max)
		return max;
	else
		return in;
}

slong truncation_order (arb_t eta, arb_t alpha, slong bits)
{
	/* Compute the number of coefficients necessary to obtain a truncation precision of 2^-bits */
	if (!arb_is_finite(eta) || !arb_is_finite(alpha))
		return 2*bits;
	if (!arb_lt(eta, alpha))
		return 0;
	slong n = 0;
	arb_t r;
	arb_init(r);
	/* Prepare r as a specific convex combination */
	arb_add(r, eta, alpha, bits);
	arb_mul_2exp_si(r, r, -1);
	/* The value of r is uninteresting. Only the quotient eta/r is relevant. */
	arb_div(r, eta, r, bits);

	/* Compute the formula for n */
	arb_t N, temp;
	arb_init(N);
	arb_init(temp);

	arb_one(temp);
	arb_sub(temp, temp, r, bits);
	arb_log(N, temp, bits);

	arb_const_log2(temp, bits);
	arb_submul_si(N, temp, bits, bits);

	arb_log(temp, r, bits);
	arb_div(N, N, temp, bits);

	arb_ceil(N, N, bits);
	arb_get_unique_fmpz(&n, N);

	arb_clear(N);
	arb_clear(temp);
	arb_clear(r);
	return n;
}

void acb_ode_solve_fuchs (acb_poly_t res, acb_ode_t ODE, slong num_of_coeffs, slong bits)
{
	/* Iteratively compute the first num_of_coeffs coefficients of the power series solution of the ODE around zero */
	acb_t temp1; acb_init(temp1);
	acb_t temp2; acb_init(temp2);
	acb_t new_coeff; acb_init(new_coeff);

	fmpz_t fac; fmpz_init(fac);
	slong i_min, i_max;
	slong v = acb_ode_valuation(ODE);

	acb_poly_fit_length(res, num_of_coeffs);
	for (slong b_max = v; b_max < num_of_coeffs; b_max++)
	{
		acb_zero(new_coeff);
		slong exp = b_max - v;
		/* Loop through the known coefficients of the power series */
		slong b_min = clamp(exp - degree(ODE), 0, b_max);
		acb_set(temp2, diff_eq_coeff(ODE, 0, exp-b_min));
		slong b = b_min;
		do {
			acb_poly_get_coeff_acb(temp1, res, b);
			acb_mul(temp2, temp1, temp2, bits);
			acb_sub(new_coeff, new_coeff, temp2, bits);

			acb_zero(temp2);
			fmpz_one(fac);
			/* Loop through the polynomials */
			b++;
			i_min = clamp(b - exp, 0, v);
			i_max = clamp(b - b_min, 0, order(ODE));
			for (slong i = i_min; i <= i_max; i++)
			{
				acb_set(temp1, diff_eq_coeff(ODE, i, i + exp - b));
				acb_mul_fmpz(temp1, temp1, fac, bits);
				acb_add(temp2, temp2, temp1, bits);
				fmpz_mul_si(fac, fac, b-i);
			}
			fmpz_rfac_uiui(fac, b-i_min+1, i_min);
			acb_mul_fmpz(temp2, temp2, fac, bits);
		} while (b != b_max);
		acb_div(new_coeff, new_coeff, temp2, bits);
		acb_poly_set_coeff_acb(res, b_max, new_coeff);
	}

	acb_clear(new_coeff);
	acb_clear(temp1);
	acb_clear(temp2);
	fmpz_clear(fac);
}

void analytic_continuation (acb_poly_t res, acb_ode_t ODE, acb_srcptr path,
		slong len, slong num_of_coeffs, slong bits)
{
	/* Evaluate a solution along the given piecewise linear path */
	acb_t a; acb_init(a);
	acb_ode_t ODE_shift;
	acb_ode_init_blank(ODE_shift, degree(ODE), order(ODE));
	for (slong time = 0; time+1 < len; time++)
	{
		acb_ode_shift(ODE_shift, ODE, path+time, bits);
		acb_ode_solve_fuchs(res, ODE_shift, num_of_coeffs, bits);
		acb_sub(a, path+time+1, path+time, bits);
		acb_poly_taylor_shift(res, res, a, bits);
	}
	acb_ode_clear(ODE_shift);
	acb_clear(a);
}

void find_monodromy_matrix (acb_mat_t mono, acb_ode_t ODE, slong bits)
{
	/* Choose a path for the analytic continuation */
	arb_t rad_of_conv;
	arb_init(rad_of_conv);
	radius_of_convergence(rad_of_conv, ODE, 40, bits);
	if (arb_is_zero(rad_of_conv))
		return;
	if (!arb_is_finite(rad_of_conv))
		arb_one(rad_of_conv);
	else
		arb_div_si(rad_of_conv, rad_of_conv, 2, bits);
	arb_get_mid_arb(rad_of_conv, rad_of_conv);

	/* Now initialize the path to move along and a polynomial to store the power series */
	slong steps = 256;
	acb_ptr path = _acb_vec_init(steps+1);
	acb_poly_t res;
	acb_poly_init(res);
	_acb_vec_unit_roots(path, steps, steps, bits);
	_acb_vec_scalar_mul_arb(path, path, steps, rad_of_conv, bits);

	/* Find the number of coefficients necessary every time
	 * (I am abusing path+steps for this because I can) */
	acb_sub(path+steps, path, path+1, bits);
	acb_abs(acb_realref(path+steps), path+steps, bits);
	slong num_of_coeffs = truncation_order(acb_realref(path+steps), rad_of_conv, bits);
	acb_set(path+steps, path);
	arb_clear(rad_of_conv);

	/* Compute the function along the chosen path */
	for (slong i = 0; i < order(ODE); i++)
	{
		acb_poly_zero(res);
		acb_poly_set_coeff_si(res, i, 1);
		analytic_continuation(res, ODE, path, steps+1, num_of_coeffs, bits);
		_acb_vec_set(acb_mat_entry_ptr(mono, i, 0), acb_poly_get_coeff_ptr(res, 0), order(ODE));
	}
	acb_mat_transpose(mono, mono);
	acb_poly_clear(res);
	_acb_vec_clear(path, steps+1);
}

void radius_of_convergence (arb_t rad_of_conv, acb_ode_t ODE, slong n, slong bits)
{
	/* Find the radius of convergence of the power series expansion */
	slong deg = 0;
	acb_ptr P = _acb_vec_init(degree(ODE)+1);
	_acb_vec_set(P, diff_eq_poly(ODE, order(ODE)), degree(ODE)+1);
	deg = degree(ODE);
	for (slong i = 0; i < deg/2; i++)
		acb_swap(P+i, P+deg-i);
	while (acb_contains_zero(P+deg))
		deg--;
	if (deg == 0)
	{
		/* There are no singularities (outside zero, possibly) */
		arb_indeterminate(rad_of_conv);
		_acb_vec_clear(P, degree(ODE)+1);
		return;
	}

	/* Graeffe Transform */
	for (slong counter = 0; counter < n; counter++)
		_acb_poly_graeffe_transform(P, P, deg+1, bits);

	/* Evaluation */
	fmpq_t root;
	fmpq_init(root);
	fmpq_one(root);
	fmpq_mul_2exp(root, root, n);
	fmpq_inv(root, root);
	_acb_poly_root_bound_fujiwara(arb_radref(rad_of_conv), P, deg+1);
	arb_get_rad_arb(rad_of_conv, rad_of_conv);
	arb_pow_fmpq(rad_of_conv, rad_of_conv, root, bits);
	arb_inv(rad_of_conv, rad_of_conv, bits);

	/* Error bound */
	arb_t radius;
	arb_init(radius);
	arb_set_si(radius, 2*deg);
	arb_pow_fmpq(radius, radius, root, bits);
	arb_sub_si(radius, radius, 1, bits);
	arb_mul(radius, rad_of_conv, radius, bits);
	arb_add_error(rad_of_conv, radius);

	arb_clear(radius);
	fmpq_clear(root);
	_acb_vec_clear(P, degree(ODE)+1);
}
