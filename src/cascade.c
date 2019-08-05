#include "cascade.h"

slong convergence_tolerance = 2;

slong truncation_order (arb_t eta, arb_t alpha, slong bits)
{
    /* Compute the number of coefficients necessary to obtain a truncation precision of 2^-bits */
    if (!arb_is_finite(eta) || !arb_is_finite(alpha))
        return convergence_tolerance*bits;
    if (!arb_lt(eta,alpha))
        return 0;
    slong n = 0;
    arb_t r;
    arb_init(r);
    /* Prepare r as a specific convex combination */
    arb_mul_si(r,eta,801,bits);
    arb_addmul_si(r,alpha,199,bits);
    arb_div_si(r,r,1000,bits);
    /* The value of r is uninteresting. Only the quotient eta/r is relevant. */
    arb_div(r,eta,r,bits);

    /* Compute the formula for n */
    arb_t N, temp;
    arb_init(N);
    arb_init(temp);

    arb_sub_si(temp,r,1,bits);
    arb_neg(temp,temp);
    arb_log(N,temp,bits);

    arb_const_log2(temp,bits);
    arb_submul_si(N,temp,bits,bits);

    arb_log(temp,r,bits);
    arb_div(N,N,temp,bits);

    arb_ceil(N,N,bits);
    if (arb_get_unique_fmpz(&n,N) == 0)
        n = convergence_tolerance*bits;

    arb_clear(N);
    arb_clear(temp);
    arb_clear(r);
    return n;
}

ulong find_power_series(acb_ode_t ODE, acb_t in, arb_t rad, slong bits)
{
    /* Iteratively compute the summands of the power series solution near z=in (where z0 = 0).
     * Only converges when z0=0 is an ordinary point and and B(0,|in|) contains no singularities. */
    if (ODE == NULL)
        return 0;
    if (acb_poly_is_zero(ODE->solution))
        return 1;

    /* First bound the number of coefficientes that are computed.
     * This assures a proper enclosure on the one hand as well as a terminating do-while-loop on the other. */
    arb_t eta;
    arb_init(eta);
    if (in != NULL)
        acb_abs(eta,in,bits);
    else
        arb_indeterminate(eta);
    slong num_of_coeffs = truncation_order(eta,rad,bits);
    arb_clear(eta);
    if (num_of_coeffs <= 0)
        return 0;

    /* Only now does it make sense to initialise variables */
    acb_t temp; acb_init(temp);
    acb_t newCoeff; acb_init(newCoeff); /* Stores the latest coefficient which is being computed */
    acb_t oldCoeff; acb_init(oldCoeff); /* Stores the old coefficient which is looped over */
    mag_t ubound;   mag_init(ubound);   /* Stores an upper bound for the new coefficient */

    /* Now compute the recursion */
    slong realError = 0, imagError = 0;
    slong minIndex, polyIndex, newIndex = 0;
    for (int oldIndex = 0; oldIndex < order(ODE); oldIndex++)
    {
        acb_poly_get_coeff_acb(oldCoeff,ODE->solution,oldIndex);
        if (realError < arb_rel_error_bits(acb_realref(oldCoeff)))
            realError = arb_rel_error_bits(acb_realref(oldCoeff));
        if (realError < convergence_tolerance*bits)
            realError = convergence_tolerance*bits;

        if (imagError < arb_rel_error_bits(acb_imagref(oldCoeff)))
            imagError = arb_rel_error_bits(acb_imagref(oldCoeff));
        if (imagError < convergence_tolerance*bits)
            imagError = convergence_tolerance*bits;
    }
    do {
        acb_zero(newCoeff);
        minIndex = (newIndex > degree(ODE)) ? newIndex - degree(ODE) : 0;

        for (slong oldIndex = newIndex+order(ODE)-1; oldIndex >= minIndex; oldIndex--)
        {
            /* Loop through the polynomials */
            polyIndex = order(ODE);
            /* No more than degree(ODE) terms can contribute: */
            if (polyIndex > degree(ODE) + oldIndex - newIndex)
                polyIndex = degree(ODE) + oldIndex - newIndex;
            if (polyIndex > oldIndex)
                polyIndex = oldIndex;

            acb_poly_get_coeff_acb(oldCoeff,ODE->solution,oldIndex);
            acb_set(temp,diff_eq_coeff(ODE,polyIndex,polyIndex+newIndex-oldIndex));
            for (polyIndex--; polyIndex >= 0; polyIndex--)
            {
                acb_mul_si(temp,temp,oldIndex - polyIndex,bits);
                if (polyIndex + newIndex - oldIndex < 0)
                    continue;
                acb_add(temp,temp,diff_eq_coeff(ODE,polyIndex,polyIndex+newIndex-oldIndex),bits);
            }
            acb_addmul(newCoeff,oldCoeff,temp,bits);
        }
        /* Divide by the coefficient of a_b where b = newIndex + order(ODE) */
        fmpz_rfac_uiui(&minIndex,newIndex+1,order(ODE)); /* Notice that minIndex is no longer needed at this point
                                                            and can therefore act as a temporary variable */
        acb_mul_si(temp,diff_eq_coeff(ODE,order(ODE),0),-1*minIndex,bits);
        acb_div(newCoeff,newCoeff,temp,bits);

        arb_get_mag(ubound,acb_realref(newCoeff));
        mag_mul_2exp_si(arb_radref(acb_realref(newCoeff)),ubound,-realError);
        arb_get_mag(ubound,acb_imagref(newCoeff));
        mag_mul_2exp_si(arb_radref(acb_imagref(newCoeff)),ubound,-imagError);

        if (acb_is_zero(newCoeff))
            continue;
        if (!acb_is_finite(newCoeff))
        {
            flint_printf("A coefficient was evaluated to be NaN. Aborting.\n");
            newIndex = NON_CONVERGENT;
            break;
        }
        acb_poly_set_coeff_acb(ODE->solution,newIndex+order(ODE),newCoeff);
    } while (++newIndex < num_of_coeffs);
    if (newIndex == NON_CONVERGENT)
        acb_ode_dump(ODE,"odedump.txt");
    else
        acb_poly_truncate(ODE->solution,order(ODE)+newIndex+1);

    mag_clear(ubound);
    acb_clear(newCoeff); acb_clear(oldCoeff);
    acb_clear(temp);
    return newIndex;
}

int checkODE (acb_poly_t *polys, acb_ode_t ODE, acb_t z, slong bits)
{
    acb_poly_t result, polyder, summand;
    acb_poly_init(polyder);
    acb_poly_init(summand);
    acb_poly_init(result);

    acb_t res; acb_init(res);
    mag_t absValue; mag_init(absValue);
    int incorrect = 0;

    acb_poly_set(polyder,ODE->solution);
    for (slong n = 0; n <= order(ODE); n++)
    {
        if (polys[n] == NULL)
            continue;
        acb_poly_mul(summand,polyder,polys[n],bits);
        acb_poly_add(result,result,summand,bits);
        acb_poly_derivative(polyder,polyder,bits);
    }
    acb_poly_evaluate(res,result,z,bits);
    acb_get_mag(absValue,res);
    if (mag_cmp_2exp_si(absValue,-bits*0.90) >= 0)
    {
        incorrect = 1;
        acb_ode_dump(ODE,"odedump.txt");
    }

    acb_poly_clear(summand);
    acb_poly_clear(polyder);
    acb_poly_clear(result);

    mag_clear(absValue);
    acb_clear(res);
    return incorrect;
}

void analytic_continuation (acb_t res, acb_ode_t ODE, acb_srcptr path, slong len, slong bits, int output_solution)
{
    /* Evaluate a solution along the given piecewise linear path */
    acb_t a; acb_init(a);
    arb_t rad; arb_init(rad);
    acb_set(a,path);

    slong time = 0;
    for (; time+1 < len; time++)
    {
        acb_ode_shift(ODE,a,bits);
        acb_sub(a,path+time+1,path+time,bits);
        acb_abs(rad,path+time,bits);
        if (find_power_series(ODE,a,rad,bits) == 0)
        {
            flint_printf("The power series expansion did not converge from ");
            acb_printd(path+time,10);
            flint_printf(" to ");
            acb_printd(path+time+1,10);
            flint_printf(" where t = %w.\n",time);
            time++;
            break;
        }
        acb_poly_taylor_shift(ODE->solution,ODE->solution,a,bits);
    }
    /* Move it back to the point of origin */
    acb_neg(a,path+time-1);
    acb_ode_shift(ODE,a,bits);
    if (output_solution == TRUE)
    {
        if (acb_poly_length(ODE->solution) < order(ODE))
            _acb_vec_set(res,acb_poly_get_coeff_ptr(ODE->solution,0),acb_poly_length(ODE->solution));
        else
            _acb_vec_set(res,acb_poly_get_coeff_ptr(ODE->solution,0),order(ODE));
    }

    acb_clear(a);
    return;
}

void find_monodromy_matrix (acb_mat_t monodromy, acb_ode_t ODE, acb_t z0, slong bits)
{
    if (ODE == NULL)
    {
        flint_printf("The ODE is the NULL-pointer. Please confirm input.\n");
        return;
    }
    slong steps = 32;
    acb_ptr path = _acb_vec_init(steps+1);
    acb_t radOfConv;
    acb_init(radOfConv);
    /* Move to the given singularity */
    if (z0 != NULL && acb_is_finite(z0))
        acb_ode_shift(ODE,z0,bits);

    /* Choose a path for the analytic continuation */
    radiusOfConvergence(acb_realref(radOfConv),ODE,bits);
    if (acb_is_zero(radOfConv))
        return;
    acb_div_si(radOfConv,radOfConv,convergence_tolerance,bits);

    _acb_vec_unit_roots(path, steps, steps, bits);
    acb_one(path+steps);
    _acb_vec_scalar_mul(path,path,steps+1,radOfConv,bits);

    /* Compute the function along the chosen path */
    for (slong i = 0; i < order(ODE); i++)
    {
        acb_poly_zero(ODE->solution);
        acb_poly_set_coeff_si(ODE->solution,i,1);
        analytic_continuation(acb_mat_entry(monodromy,i,0),ODE,path,steps+1,bits,TRUE);
    }
    acb_mat_transpose(monodromy,monodromy);
    /* Move back to the start */
    if (z0 != NULL && acb_is_finite(z0))
    {
        acb_neg(z0,z0);
        acb_ode_shift(ODE,z0,bits);
        acb_neg(z0,z0);
    }
    acb_clear(radOfConv);
    _acb_vec_clear(path,steps+1);
    return;
}

void graeffe_transform(acb_ptr dest, acb_srcptr src, slong len, slong bits)
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

void radiusOfConvergence(arb_t radOfConv, acb_ode_t ODE, slong bits)
{
    /* Find the radius of cervegence of the power series expansion */
    if (ODE == NULL)
    {
        arb_indeterminate(radOfConv);
        return;
    }
    slong valuation = 0;
    arb_pos_inf(radOfConv);
    arb_t R;
    arb_init(R);
    acb_ptr P = _acb_vec_init(degree(ODE)+1);
    arb_ptr Q = _arb_vec_init(degree(ODE)+1);
    _acb_vec_set(P,diff_eq_poly(ODE,order(ODE)),degree(ODE)+1);
    while (acb_contains_zero(diff_eq_coeff(ODE,order(ODE),valuation)))
        valuation++;
    _acb_poly_shift_right(P,P,degree(ODE)+1,valuation);
    slong it = 1; /* A counter to avoid infinite loops */
    while (it <= 1048576)
    {
        _acb_poly_majorant(Q,P,degree(ODE)+1,bits);
        for (slong i = 0; i <= degree(ODE); i++)
        {
            if (arb_is_zero(Q+i))
                continue;
            arb_div(R,Q,Q+i,bits);
            arb_root_ui(R,R,i,bits);
            if (arb_lt(R,radOfConv))
                arb_set(radOfConv,R);
        }
        arb_div_ui(radOfConv,radOfConv,2,bits);
        arb_root_ui(radOfConv,radOfConv,it,bits);
        graeffe_transform(P,P,degree(ODE)+1,bits);
        it*=2;
    }
    arb_printd(radOfConv,20);
    flint_printf("\n");
    arb_clear(R);
    _arb_vec_clear(Q,degree(ODE)+1);
    _acb_vec_clear(P,degree(ODE)+1);
    return;
}
