#include "cascade.h"

slong convergence_tolerance = 2;

ulong find_power_series(acb_ode_t ODE, acb_t in, slong bits) {
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
    slong realError, imagError;
    slong newIndex = 0;
    do {
        realError = imagError = 0;
        num_of_nonzero = 0;
        acb_zero(newCoeff);
        slong minIndex = (newIndex > degree(ODE)) ? newIndex - degree(ODE) : 0;

        for (slong oldIndex = newIndex+order(ODE)-1; oldIndex >= minIndex; oldIndex--)
        {
            acb_poly_get_coeff_acb(oldCoeff,ODE->solution,oldIndex);

            /* Ignore zero terms immediately: */
            if (acb_is_zero(oldCoeff))
                continue;
            if (realError < arb_rel_error_bits(acb_realref(oldCoeff)))
                realError = arb_rel_error_bits(acb_realref(oldCoeff));
            if (realError < convergence_tolerance*bits)
                realError = convergence_tolerance*bits;

            if (imagError < arb_rel_error_bits(acb_imagref(oldCoeff)))
                imagError = arb_rel_error_bits(acb_imagref(oldCoeff));
            if (imagError < convergence_tolerance*bits)
                imagError = convergence_tolerance*bits;

            /* Check for proper convergence of the summands */
            acb_mul(temp,oldCoeff,power,bits);
            acb_get_mag(ubound,temp);
            if (mag_cmp_2exp_si(ubound,-bits) >= 0)
                num_of_nonzero++;
            acb_zero(temp);

            /* Loop through the polynomials */
            for (slong polyIndex = order(ODE); polyIndex >= 0; polyIndex--)
            {
                /* No more than degree(ODE) terms can contribute: */
                if (polyIndex > degree(ODE) + oldIndex - newIndex)
                    polyIndex = degree(ODE) + oldIndex - newIndex;
                if (polyIndex > oldIndex)
                    polyIndex = oldIndex;
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
        acb_ode_dump(ODE);

    mag_clear(ubound);
    acb_clear(newCoeff); acb_clear(oldCoeff);
    acb_clear(power); acb_clear(temp);
    return newIndex;
}

void analytic_continuation (acb_t res, acb_ode_t ODE, acb_srcptr path, slong len, slong bits, int output_solution)
{
    acb_t a; acb_init(a);
    slong time = 0;

    /* Evaluate a solution along the given piecewise linear path */
    acb_set(a,path);
    for (; time+1 < len; time++)
    {
        acb_ode_shift(ODE,a,bits);
        acb_sub(a,path+time+1,path+time,bits);
        acb_poly_truncate(ODE->solution,order(ODE));
        if (find_power_series(ODE,a,bits) == 0)
        {
            flint_printf("The power series expansion did not converge from ");
            acb_printd(path+time,10);
            flint_printf(" to ");
            acb_printd(path+time+1,10);
            flint_printf(" where t = %w.\n",time);
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

int checkODE (acb_poly_t *polys, acb_ode_t ODE, acb_t z, slong bits) {
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
        flint_printf("The differential equation has not been solved correctly. These are the coefficients:\n");
        acb_poly_printd(ODE->solution,10);
        flint_printf("\n, which evaluates to ");
        mag_print(absValue);
        flint_printf(" at ");
        acb_printd(z,10);
        flint_printf(".\n");
    }

    acb_poly_clear(summand);
    acb_poly_clear(polyder);
    acb_poly_clear(result);

    mag_clear(absValue);
    acb_clear(res);
    return incorrect;
}

void find_monodromy_matrix (acb_mat_t monodromy, acb_ode_t ODE, acb_ptr path, slong len, slong bits) {
    acb_mat_init(monodromy,order(ODE),order(ODE));
    if (monodromy == NULL) /* Further Error checking elsewhere */
        return;
    for (slong i = 0; i < order(ODE); i++)
    {
        acb_poly_zero(ODE->solution);
        acb_poly_set_coeff_si(ODE->solution,i,1);
        analytic_continuation(acb_mat_entry(monodromy,i,0),ODE,path,len,bits,TRUE);
    }
    acb_t determinant; acb_init(determinant);
    acb_mat_det(determinant,monodromy,bits);
    if (order(ODE) > 1)
    {
        flint_printf("The monodromy matrix has the determinant ");
        acb_printn(determinant,bits,ARB_STR_CONDENSE * 50);
        flint_printf("\n");
    }
    acb_clear(determinant);
    return;
}

void acb_ode_dump(acb_ode_t ODE)
{
    FILE *out = fopen("data/odedump.txt","w");
    if (out == NULL)
        return;
    flint_fprintf(out,"Order: %w\nDegree: %w\n",order(ODE),degree(ODE));
    for (slong i = 0; i <= order(ODE); i++)
    {
        for (slong j = 0; j <= degree(ODE); j++)
        {
            acb_fprintd(out,diff_eq_coeff(ODE,i,j),10);
            flint_fprintf(out,"\t");
        }
        flint_fprintf(out,"\n");
    }
    flint_fprintf(out,"\n");
    acb_poly_fprintd(out,ODE->solution,10);
    fclose(out);
    return;
}

void entry_point (const ulong maxOrder, slong bits, double z_val, const char *file) {
    bits = bits*3.32193 + 5;
    slong steps = 128;

    ulong numOfPols;
    acb_poly_t *polys = acb_ode_fread(&numOfPols,file,maxOrder,10*bits);
    if (polys == NULL)
        return;

    acb_t res, z;
    acb_init(res); acb_init(z);
    acb_set_d(z,z_val);
    acb_mat_t monodromy;

    acb_ptr path = _acb_vec_init(steps+1);
    for (slong i = 0; i < steps; i++) {
        acb_set_si(path+i,i);
        acb_div_ui(path+i,path+i,steps/2,bits);
        acb_exp_pi_i(path+i,path+i,bits);
    }
    acb_set(path+steps,path);
    _acb_vec_scalar_div_ui(path,path,steps+1,8,steps/2*bits);
    flint_printf("Path has been initialised!\n");

    acb_ode_t ODE = acb_ode_init(polys,NULL,numOfPols);
    if (ODE != NULL)
    {
        //~ acb_poly_set_coeff_si(ODE->solution,0,1);
        //~ find_power_series(ODE,z,bits);
        //~ checkODE(polys,ODE,z,bits);
        find_monodromy_matrix(monodromy,ODE,path,steps+1,bits);
        if (monodromy == NULL)
            flint_printf("Could not allocate memory. Please try again.\n");
        else
        {
            acb_mat_printd(monodromy,bits/33.2193);
            acb_mat_clear(monodromy);
        }
    }

    for (ulong i = 0; i <= numOfPols; i++)
        if (polys[i] != NULL)
            acb_poly_clear(polys[i]);
    flint_free(polys);
    acb_ode_clear(ODE);
    acb_clear(res); acb_clear(z);
    _acb_vec_clear(path,steps+1);
    flint_cleanup();
    return;
}

int main (int argc, char **argv) {
    entry_point((argc >= 3) ? atol(argv[2]) : 1,
                  (argc >= 4) ? atol(argv[3]) : 50,
                  (argc >= 5) ? atof(argv[4]) : 1,
                  argv[1]);
    flint_cleanup();
    return 0;
}
