#include "cascade.h"

slong convergence_tolerance = 2;

ulong find_power_series(acb_ode_t ODE, acb_t in, slong bits)
{
    /* Iteratively compute the summands of the power series solution near z=in (where z0 = 0).
     * Only converges when z0=0 is an ordinary point and and B(0,|in|) contains no singularities.
     * This is not tested here! */
    if (ODE == NULL)
        return 0;
    if (acb_poly_is_zero(ODE->solution))
        return 1;
    acb_t temp; acb_init(temp);
    acb_t newCoeff; acb_init(newCoeff); /* Stores the latest coefficient which is being computed */
    acb_t oldCoeff; acb_init(oldCoeff); /* Stores the old coefficient which is looped over */
    acb_t power;    acb_init(power);    /* Stores the current power of in */
    mag_t ubound;   mag_init(ubound);   /* Stores an upper bound for each old coefficient */
    acb_one(power);

    slong num_of_nonzero;
    slong realError = 0, imagError = 0;
    slong minIndex, maxPoly, newIndex = 0;
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
        num_of_nonzero = 0;
        acb_zero(newCoeff);
        minIndex = (newIndex > degree(ODE)) ? newIndex - degree(ODE) : 0;

        for (slong oldIndex = newIndex+order(ODE)-1; oldIndex >= minIndex; oldIndex--)
        {
            acb_poly_get_coeff_acb(oldCoeff,ODE->solution,oldIndex);

            /* Ignore zero terms immediately: */
            if (acb_is_zero(oldCoeff))
                continue;

            /* Check for proper convergence of the summands */
            acb_mul(temp,oldCoeff,power,bits);
            acb_get_mag(ubound,temp);
            if (mag_cmp_2exp_si(ubound,-bits) >= 0)
                num_of_nonzero++;
            acb_zero(temp);

            /* Loop through the polynomials */
            maxPoly = order(ODE);
            /* No more than degree(ODE) terms can contribute: */
            if (maxPoly > degree(ODE) + oldIndex - newIndex)
                maxPoly = degree(ODE) + oldIndex - newIndex;
            if (maxPoly > oldIndex)
                maxPoly = oldIndex;
            for (slong polyIndex = maxPoly; polyIndex >= 0; polyIndex--)
            {
                if (polyIndex + newIndex - oldIndex >= 0)
                    acb_add(temp,temp,diff_eq_coeff(ODE,polyIndex,polyIndex+newIndex-oldIndex),bits);

                if (polyIndex != 0)
                    acb_mul_si(temp,temp,oldIndex - polyIndex + 1,bits);
            }
            acb_addmul(newCoeff,oldCoeff,temp,bits);
        }
        /* Divide by the coefficient of a_b where b = newIndex + order(ODE) */
        acb_set_ui(temp,newIndex+1);
        acb_rising_ui(temp,temp,order(ODE),bits);
        acb_mul(temp,temp,diff_eq_coeff(ODE,order(ODE),0),bits);
        acb_neg(temp,temp);
        acb_div(newCoeff,newCoeff,temp,bits);

        arb_get_mag(ubound,acb_realref(newCoeff));
        mag_mul_2exp_si(arb_radref(acb_realref(newCoeff)),ubound,-realError);
        arb_get_mag(ubound,acb_imagref(newCoeff));
        mag_mul_2exp_si(arb_radref(acb_imagref(newCoeff)),ubound,-imagError);

        if (!acb_is_finite(newCoeff))
        {
            flint_printf("A coefficient was evaluated to be NaN. Aborting.\n");
            newIndex = NON_CONVERGENT;
            break;
        }
        acb_poly_set_coeff_acb(ODE->solution,newIndex+order(ODE),newCoeff);
        if (in != NULL)
            acb_mul(power,power,in,bits);

        if (++newIndex >= convergence_tolerance*bits)
        {
            flint_printf("%w summands have been computed and no convergence was achieved. Aborting.\n",newIndex);
            flint_printf("Target was %w bits.\n",bits);
            newIndex = NON_CONVERGENT;
            break;
        }
    } while (num_of_nonzero != 0);
    if (newIndex == NON_CONVERGENT)
        acb_ode_dump(ODE,"odedump.txt");
    else
        acb_poly_truncate(ODE->solution,order(ODE)+newIndex+1);

    mag_clear(ubound);
    acb_clear(newCoeff); acb_clear(oldCoeff);
    acb_clear(power); acb_clear(temp);
    return newIndex;
}

void analytic_continuation (acb_t res, acb_ode_t ODE, acb_srcptr path, slong len, slong bits, int output_solution)
{
    /* Evaluate a solution along the given piecewise linear path */
    acb_t a; acb_init(a);
    acb_set(a,path);

    slong time = 0;
    for (; time+1 < len; time++)
    {
        acb_ode_shift(ODE,a,bits);
        acb_sub(a,path+time+1,path+time,bits);
        if (find_power_series(ODE,a,bits) == 0)
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
    if (acb_is_finite(z0))
        acb_ode_shift(ODE,z0,bits);

    /* Choose a path for the analytic continuation */
    radiusOfConvergence(ODE,arb_midref(acb_realref(radOfConv)),bits);
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
    acb_clear(radOfConv);
    _acb_vec_clear(path,steps+1);
    /* Move back to the start */
    if (acb_is_finite(z0))
    {
        acb_neg(z0,z0);
        acb_ode_shift(ODE,z0,bits);
    }
    return;
}
