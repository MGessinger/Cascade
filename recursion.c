#include "recursion.h"

ulong find_power_series_regular(acb_ode_t ODE, acb_t in, slong bits) {
    /* Iteratively compute the summands of the power series solution near in (where z0 = 0)
     * Only converges when z0=0 is an ordinary point and and B(0,|in|) contains no singularities.
     * This is not tested here! */
    acb_t temp; acb_init(temp);
    acb_t newCoeff; acb_init(newCoeff); /* Stores the latest coefficient which is being computed */
    acb_t oldCoeff; acb_init(oldCoeff); /* Stores the old coefficient which is looped over */
    acb_t power;    acb_init(power);    /* Stores the current power of in */
    mag_t ubound;   mag_init(ubound);   /* Stores an upper bound for each old coefficient */

    if (in == NULL) /* No input means no output. */
    {
        flint_printf("No input was provided so no output will be computed.\n");
    }
    else
        acb_one(power);

    slong num_of_nonzero;
    slong newIndex = 0;
    slong currentError, maxRealError, maxImagError;
    do {
        num_of_nonzero = 0;
        maxRealError = maxImagError = -bits;
        acb_zero(newCoeff);
        slong minIndex = (newIndex > degree(ODE)) ? newIndex - degree(ODE) : 0;

        for (slong oldIndex = newIndex+order(ODE)-1; oldIndex >= minIndex; oldIndex--)
        {
            acb_poly_get_coeff_acb(oldCoeff,ODE->series,oldIndex);

            /* Ignore zero terms immediately: */
            if (acb_is_zero(oldCoeff))
                continue;

            /* Check for proper convergence of the summands */
            acb_mul(temp,oldCoeff,power,bits);
            acb_get_mag(ubound,temp);
            if (mag_cmp_2exp_si(ubound,-bits) >= 0)
                num_of_nonzero++;
            acb_zero(temp);

            /* Check error bounds */
            currentError = arb_rel_error_bits(acb_realref(oldCoeff));
            if (currentError >= maxRealError)
                maxRealError = currentError;
            currentError = arb_rel_error_bits(acb_imagref(oldCoeff));
            if (currentError >= maxImagError)
                maxImagError = currentError;

            /* Loop through the polynomials */
            for (slong polyIndex = order(ODE); polyIndex >= 0; polyIndex--)
            {
                /* No more than degree(ODE) terms can contribute: */
                if (polyIndex > degree(ODE) + oldIndex - newIndex)
                    polyIndex = degree(ODE) + oldIndex - newIndex;
                if (polyIndex > oldIndex)
                    polyIndex = oldIndex;
                if (polyIndex + newIndex - oldIndex < 0)
                    break;

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

        if (!acb_is_finite(newCoeff))
        {
            flint_printf("A coefficient was evaluated to be NaN. Aborting.\n");
            newIndex = NON_CONVERGENT;
            break;
        }
        acb_poly_set_coeff_acb(ODE->series,newIndex+order(ODE),newCoeff);
        acb_mul(power,power,in,bits);

        if (++newIndex >= 2*bits)
        {
            flint_printf("%w summands have been computed and no convergence was achieved. Aborting.\n",newIndex);
            newIndex = NON_CONVERGENT;
            break;
        }
    } while (num_of_nonzero != 0);
    if (newIndex == NON_CONVERGENT)
        acb_ode_dump(ODE);

    flint_printf("The maximum relative Error amongst all coefficients was (%w,%w).\n",maxRealError,maxImagError);

    mag_clear(ubound);
    acb_clear(newCoeff); acb_clear(oldCoeff);
    acb_clear(power); acb_clear(temp);

    return newIndex;
}

void analytic_continuation (acb_t res, acb_ode_t ODE_in, acb_srcptr path, slong len, slong digits, int output_series)
{
    acb_ode_t ODE = acb_ode_set(NULL,ODE_in);
    if (ODE == NULL) return;
    acb_t a; acb_init(a);
    slong accuracy, bits = digits/1.414; /* Start with a much smaller value to get a good estimate very quickly */

    acb_set(a,path);
    for (slong i = 1; i < len; i++)
    {
        acb_ode_shift(ODE,a,bits);
        acb_sub(a,path+i,path+(i-1),bits);
        if (find_power_series_regular(ODE,a,bits) == 0)
            break;
        acb_poly_taylor_shift(ODE->series,ODE->series,a,bits);
        accuracy = acb_rel_accuracy_bits(res);
        /* If the precision of the result is not good enough, start over */
        if (accuracy < digits && accuracy >= 0)
        {
            flint_printf("Ran out of precision at %w bits. ", bits);
            flint_printf("Result dropped to a decimal-accuracy of %w at t = %w.\n",accuracy,i);
            ODE = acb_ode_set(ODE,ODE_in);
            i = 0;
            bits += (digits-accuracy)+len; /* This has proven to be a good choice */
        }
        if (bits >= 1.5*len*digits)
        {
            flint_printf("The internal precision exceeded a reasonable bound. Aborting.\n");
            break;
        }
    }
    if (output_series == TRUE)
        _acb_vec_set(res,acb_poly_get_coeff_ptr(ODE->series,0),order(ODE));
    else
        flint_printf("It required a working precision of %w bits to reach the requested %w digits.\n\n",bits,digits/3.32193);

    acb_ode_clear(ODE);
    acb_clear(a);
    return;
}

int checkODE (acb_poly_t *polys, acb_ode_t ODE, slong digits) {
    acb_poly_t result, polyder, summand;
    acb_poly_init(polyder);
    acb_poly_init(summand);
    acb_poly_init(result);

    arf_t absValue;
    arf_init(absValue);
    int printed = 0;

    acb_poly_set(polyder,ODE->series);
    for (slong n = 0; n <= order(ODE); n++)
    {
        if (polys[n] == NULL)
            continue;
        acb_poly_mul(summand,polyder,polys[n],digits);
        acb_poly_add(result,result,summand,digits);
        acb_poly_derivative(polyder,polyder,digits);
    }
    for (slong n = acb_poly_degree(result); n >= 0; n--)
    {
        acb_get_abs_ubound_arf(absValue,acb_poly_get_coeff_ptr(result,n),digits);
        if (arf_cmpabs_2exp_si(absValue,-digits*0.95) >= 0)
        {
            if (printed == 0)
            {
                flint_printf("The differential equation was not solved correctly. These are the coefficients:\n");
                printed = 1;
                acb_poly_printd(ODE->series,10);
                flint_printf("\n\n");
            }
            acb_printn(acb_poly_get_coeff_ptr(result,n),digits,ARB_STR_CONDENSE * 10);
            flint_printf(" = a_%w\n",n);
        }
    }
    acb_poly_clear(summand);
    acb_poly_clear(polyder);
    acb_poly_clear(result);
    arf_clear(absValue);
    return printed;
}

void find_monodromy_matrix (acb_mat_t monodromy, acb_ode_t ODE, acb_ptr path, slong len, slong digits) {
    acb_mat_init(monodromy,order(ODE),order(ODE));
    if (monodromy == NULL) return; /* Further Error checking elsewhere */
    for (slong i = 0; i < order(ODE); i++)
    {
        acb_poly_zero(ODE->series);
        acb_poly_set_coeff_si(ODE->series,i,1);
        analytic_continuation(acb_mat_entry(monodromy,i,0),ODE,path,len,digits,TRUE);
    }
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
        for (slong j = 0; j < degree(ODE); j++)
        {
            acb_fprintd(out,diff_eq_coeff(ODE,i,j),10);
            flint_fprintf(out,"\t");
        }
        flint_fprintf(out,"\n");
    }
    flint_fprintf(out,"\n");
    acb_poly_fprintd(out,ODE->series,10);
    flint_fprintf(out,"\n\n\n");
    fclose(out);
    return;
}

void entry_point (ulong maxOrder, slong digits, slong z_val, const char *file) {
    digits = digits*3.32193; /* Multiply by log_2(10) to obtain the number of binary digits */
    slong steps = 256;

    /* Polynomials */
    ulong numOfPols;
    acb_poly_t *polys = acb_ode_fread(&numOfPols,file,maxOrder,10*digits);
    if (polys == NULL)
        return;

    /* Input/Output */
    acb_t res, z;
    acb_init(res); acb_init(z);
    acb_set_si(z,z_val);
    acb_mat_t monodromy;

    /* Path */
    acb_ptr path = _acb_vec_init(steps+1);
    for (slong i = 0; i < steps; i++) {
        acb_set_si(path+i,i);
        acb_div_ui(path+i,path+i,steps/2,digits);
        acb_exp_pi_i(path+i,path+i,steps/2*digits);
    }
    acb_set(path+steps,path);
    _acb_vec_scalar_div_ui(path,path,steps+1,4,steps/2*digits);

    acb_ode_t ODE = acb_ode_init(polys,NULL,z,numOfPols,digits);
    if (ODE != NULL)
    {
        find_monodromy_matrix(monodromy,ODE,path,steps+1,digits);
        if (monodromy == NULL)
            flint_printf("Could not allocate memory. Please try again.\n");
        else
        {
            acb_mat_printd(monodromy,7);
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
    ulong n =       (argc >= 3) ? atol(argv[2]) : 1;
    slong digits =  (argc >= 4) ? atol(argv[3]) : 100;
    slong z =       (argc >= 5) ? atol(argv[4]) : 1;

    entry_point(n,digits,z,argv[1]);
    flint_cleanup();
    return 0;
}
