#include "cascade.h"

slong truncation_order (arb_t eta, arb_t alpha, slong bits)
{
	/* Compute the number of coefficients necessary to obtain a truncation precision of 2^-bits */
	if (!arb_is_finite(eta) || !arb_is_finite(alpha))
		return 2*bits;
	if (!arb_lt(eta,alpha))
		return 0;
	slong n = 0;
	arb_t r;
	arb_init(r);
	/* Prepare r as a specific convex combination */
	arb_add(r,eta,alpha,bits);
	arb_mul_2exp_si(r,r,-1);
	/* The value of r is uninteresting. Only the quotient eta/r is relevant. */
	arb_div(r,eta,r,bits);

	/* Compute the formula for n */
	arb_t N, temp;
	arb_init(N);
	arb_init(temp);

	arb_one(temp);
	arb_sub(temp,temp,r,bits);
	arb_log(N,temp,bits);

	arb_const_log2(temp,bits);
	arb_submul_si(N,temp,bits,bits);

	arb_log(temp,r,bits);
	arb_div(N,N,temp,bits);

	arb_ceil(N,N,bits);
	arb_get_unique_fmpz(&n,N);

	arb_clear(N);
	arb_clear(temp);
	arb_clear(r);
	return n;
}

void acb_poly_graeffe_transform(acb_ptr dest, acb_srcptr src, slong len, slong bits)
{
	/* Computes the Graeffe transform of src. In- and output can be aliased. */
	slong q = (len-1)/2;
	acb_ptr pe = _acb_vec_init(q+1);
	acb_ptr po = _acb_vec_init(len);
	for (slong i = 0; i < len; i++)
	{
		if (i%2 == 0)
			acb_set(pe+(i/2),src+i);
		else
			acb_set(po+(i/2),src+i);
	}
	_acb_poly_mul(dest,po,q+1,po,q+1,bits);
	_acb_poly_shift_left(dest,dest,len-1,1);
	_acb_vec_neg(dest,dest,len);
	_acb_poly_mul(po,pe,q+1,pe,q+1,bits);
	_acb_vec_add(dest,dest,po,len,bits);

	_acb_vec_clear(pe,q+1);
	_acb_vec_clear(po,len);
	return;
}

slong find_power_series(acb_poly_t res, acb_ode_t ODE, slong numOfCoeffs, slong bits)
{
	/* Iteratively compute the first numOfCoeffs coefficients of the power series solution of the ODE around zero */
	if (!ODE || !res)
		return 0;
	if (numOfCoeffs <= 0 || acb_poly_is_zero(res))
		return 1;

	/* Only now does it make sense to initialise variables */
	acb_t temp; acb_init(temp);
	acb_t newCoeff; acb_init(newCoeff); /* Stores the latest coefficient which is being computed */

	/* Now compute the recursion */
	fmpz_t fac;
	fmpz_init(fac);
	slong minIndex, newIndex;
	slong polyMin, polyMax;
	slong offset, facStart;

	acb_poly_fit_length(res,numOfCoeffs);
	/* Negating the constant coefficient of the leading polynomial saves a few cycles later on */
	acb_neg(diff_eq_coeff(ODE,order(ODE),0),diff_eq_coeff(ODE,order(ODE),0));

	for (newIndex = order(ODE); newIndex < numOfCoeffs; newIndex++)
	{
		acb_zero(newCoeff);
		minIndex = newIndex - degree(ODE) - order(ODE);
		if (minIndex < 0)
			minIndex = 0;
		facStart = newIndex - order(ODE);

		for (slong oldIndex = newIndex-1; oldIndex >= minIndex; oldIndex--)
		{
			if (acb_poly_get_coeff_ptr(res,oldIndex) == NULL)
				continue;
			/* Loop through the polynomials */
			polyMax = oldIndex - minIndex;
			/* No more than order(ODE) terms can contribute: */
			if (polyMax > order(ODE))
				polyMax = order(ODE);

			offset = oldIndex-facStart;
			if (offset <= 0)
			{
				polyMin = 0;
				fmpz_one(fac);
			}
			else
			{
				polyMin = offset;
				fmpz_rfac_uiui(fac,facStart+1,polyMin);
			}
			acb_mul_fmpz(temp,diff_eq_coeff(ODE,polyMin,polyMin-offset),fac,bits);
			for (slong polyIndex = polyMin; polyIndex < polyMax; polyIndex++)
			{
				fmpz_mul_si(fac,fac,oldIndex - polyIndex);
				acb_addmul_fmpz(temp,diff_eq_coeff(ODE,polyIndex+1,polyIndex+1-offset),fac,bits);
			}
			if (acb_is_zero(temp))
				continue;
			acb_addmul(newCoeff,acb_poly_get_coeff_ptr(res,oldIndex),temp,bits);
		}
		/* Divide by the coefficient of a_b where b = newIndex + order(ODE) */
		fmpz_rfac_uiui(fac,facStart+1,order(ODE));
		acb_mul_fmpz(temp,diff_eq_coeff(ODE,order(ODE),0),fac,bits);
		acb_div(newCoeff,newCoeff,temp,bits);

		if (!acb_is_finite(newCoeff))
		{
			flint_printf("A coefficient was evaluated to be NaN. Aborting.\n");
			newIndex = 0;
			break;
		}
		acb_poly_set_coeff_acb(res,newIndex,newCoeff);
	}
	if (newIndex == 0)
		acb_ode_dump(ODE,"odedump.txt");

	/* Undo the negation performed earlier */
	acb_neg(diff_eq_coeff(ODE,order(ODE),0),diff_eq_coeff(ODE,order(ODE),0));

	acb_clear(newCoeff);
	acb_clear(temp);
	fmpz_clear(fac);
	return newIndex;
}

void analytic_continuation (acb_poly_t res, acb_ode_t ODE, acb_srcptr path,
		slong len, slong numOfCoeffs, slong bits)
{
	/* Evaluate a solution along the given piecewise linear path */
	acb_t a; acb_init(a);
	acb_ode_t ODE_shift = acb_ode_init_blank(degree(ODE),order(ODE));
	slong time = 0;
	for (; time+1 < len; time++)
	{
		acb_ode_shift(ODE_shift,ODE,path+time,bits);
		acb_sub(a,path+time+1,path+time,bits);
		if (find_power_series(res,ODE_shift,numOfCoeffs,bits) == 0)
		{
			flint_printf("The power series expansion did not converge from ");
			acb_printd(path+time,10);
			flint_printf(" to ");
			acb_printd(path+time+1,10);
			flint_printf(" where t = %w.\n",time);
			break;
		}
		acb_poly_taylor_shift(res,res,a,bits);
	}
	acb_ode_clear(ODE_shift);
	acb_clear(a);
	return;
}

void find_monodromy_matrix (acb_mat_t mono, acb_ode_t ODE, acb_t z0, slong bits)
{
	if (ODE == NULL)
	{
		flint_printf("No differential operator was provided. Please confirm your input.\n");
		return;
	}
	arb_t radOfConv;
	arb_init(radOfConv);
	/* Move to the given singularity */
	if (z0 != NULL && acb_is_finite(z0))
		acb_ode_shift(ODE,ODE,z0,bits);

	/* Choose a path for the analytic continuation */
	radiusOfConvergence(radOfConv,ODE,40,bits);
	if (arb_is_zero(radOfConv))
		return;
	if (!arb_is_finite(radOfConv))
		arb_one(radOfConv);
	else
		arb_div_si(radOfConv,radOfConv,2,bits);
	arb_get_mid_arb(radOfConv,radOfConv);

	/* Now initialize the path to move along and a polynomial to store the power series */
	slong steps = 256;
	acb_ptr path = _acb_vec_init(steps+1);
	acb_poly_t res;
	acb_poly_init(res);
	_acb_vec_unit_roots(path, steps, steps, bits);
	_acb_vec_scalar_mul_arb(path,path,steps,radOfConv,bits);

	/* Find the number of coefficients necessary every time (I am abusing path+steps for this because I can) */
	acb_sub(path+steps,path,path+1,bits);
	acb_abs(acb_realref(path+steps),path+steps,bits);
	slong numOfCoeffs = truncation_order(acb_realref(path+steps),radOfConv,bits);
	acb_set(path+steps,path);
	arb_clear(radOfConv);

	/* Compute the function along the chosen path */
	for (slong i = 0; i < order(ODE); i++)
	{
		acb_poly_zero(res);
		acb_poly_set_coeff_si(res,i,1);
		analytic_continuation(res,ODE,path,steps+1,numOfCoeffs,bits);
		_acb_vec_set(acb_mat_entry_ptr(mono,i,0),acb_poly_get_coeff_ptr(res,0),order(ODE));
	}
	acb_mat_transpose(mono,mono);
	/* Move back to the start */
	if (z0 != NULL && acb_is_finite(z0))
	{
		acb_neg(z0,z0);
		acb_ode_shift(ODE,ODE,z0,bits);
		acb_neg(z0,z0);
	}
	acb_poly_clear(res);
	_acb_vec_clear(path,steps+1);
	return;
}

void radiusOfConvergence(arb_t radOfConv, acb_ode_t ODE, slong n, slong bits)
{
	/* Find the radius of convergence of the power series expansion */
	if (ODE == NULL)
	{
		arb_zero(radOfConv);
		return;
	}
	slong deg = 0;
	acb_ptr P = _acb_vec_init(degree(ODE)+1);
	_acb_vec_set(P,diff_eq_poly(ODE,order(ODE)),degree(ODE)+1);
	deg = degree(ODE);
	for (slong i = 0; i < deg/2; i++)
		acb_swap(P+i,P+deg-i);
	while (acb_contains_zero(P+deg))
		deg--;
	if (deg == 0)
	{
		/* There are no singularities (outside zero, possibly) */
		arb_indeterminate(radOfConv);
		_acb_vec_clear(P,degree(ODE)+1);
		return;
	}

	/* Graeffe Transform */
	for (slong counter = 0; counter < n; counter++)
		acb_poly_graeffe_transform(P,P,deg+1,bits);

	/* Evaluation */
	fmpq_t root;
	fmpq_init(root);
	fmpq_one(root);
	fmpq_mul_2exp(root,root,n);
	fmpq_inv(root,root);
	_acb_poly_root_bound_fujiwara(arb_radref(radOfConv),P,deg+1);
	arb_get_rad_arb(radOfConv,radOfConv);
	arb_pow_fmpq(radOfConv,radOfConv,root,bits);
	arb_inv(radOfConv,radOfConv,bits);

	/* Error bound */
	arb_t radius;
	arb_init(radius);
	arb_set_si(radius,2*deg);
	arb_pow_fmpq(radius,radius,root,bits);
	arb_sub_si(radius,radius,1,bits);
	arb_mul(radius,radOfConv,radius,bits);
	arb_add_error(radOfConv,radius);

	arb_clear(radius);
	fmpq_clear(root);
	_acb_vec_clear(P,degree(ODE)+1);
	return;
}
