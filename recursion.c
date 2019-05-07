#include "recursion.h"

void analytic_continuation (acb_t res, acb_ode_t ODE_in, acb_srcptr path, slong len, slong digits, int output_series)
{
    acb_ode_t ODE = acb_ode_copy(NULL,ODE_in);
    if (ODE == NULL) return;
    acb_t a; acb_init(a);
    slong bits = digits/1.41, accuracy;

    acb_set(a,path);
    for (slong i = 1; i <= len; i++)
    {
        if (degree(ODE) >= 1)
            for (slong j = 0; j <= order(ODE); j++)
            {
                if (ODE->polys[j] == NULL) continue;
                _acb_poly_taylor_shift(ODE->polys[j],a,degree(ODE)+1,bits);
            }
        acb_sub(a,path+(i%len),path+(i-1),bits); /* When i==len, wrap back around to the beginning */
        if (find_power_series_regular(res,ODE,a,bits) == 0)
            break;
        acb_poly_taylor_shift(ODE->series,ODE->series,a,bits);
        accuracy = acb_rel_accuracy_bits(res);
        if (accuracy < digits && accuracy >= 0) /* If the precision of the result is not good enough, start over */
        {
            flint_printf("Ran out of precision at %w bits. ", bits);
            flint_printf("Result dropped to a decimal-accuracy of %w at t = %w.\n",accuracy,i);
            ODE = acb_ode_copy(ODE,ODE_in);
            i = 0;
            bits += (digits-accuracy)+len; /* This has proven to be a good choice */
        }
        if (bits >= len*digits)
        {
            flint_printf("Aborted. This does in no way make sense.\n");
            break;
        }
    }
    if (output_series == TRUE)
        _acb_vec_set(res,acb_poly_get_coeff_ptr(ODE->series,0),order(ODE));
    else
        flint_printf("It required a working precision of %w bits to reach the requested %w bits.\n\n",bits,digits);

    acb_ode_clear(ODE);
    acb_clear(a);
    return;
}

ulong find_power_series_regular(acb_t out, acb_ode_t ODE, acb_t in, slong bits) {
    /* Iteratively compute the summands of the power series solution near in (where z0 = 0)
     * Only converges when z0=0 is an ordinary point and and B(0,|in|) contains no singularities. */
    acb_t temp; acb_init(temp);
    acb_t new_coeff; acb_init(new_coeff);
    acb_t result; acb_init(result); /* Stores the result to allow for aliasing of out and in */
    arf_t err; arf_init(err); /* Measures how well the series has converged yet */

    if (in == NULL) /* No input means no output. */
    {
        flint_printf("No input was provided so no output will be computed.\n");
        out = NULL;
    }
    if (out != NULL)
    {
        acb_poly_truncate(ODE->series,order(ODE));
        acb_poly_evaluate(result,ODE->series,in,bits);
    }
    slong num_of_nonzero;
    slong newIndex = 0;
    do {
        num_of_nonzero = 0;
        acb_zero(new_coeff);
        for (slong oldIndex = newIndex+order(ODE)-1; oldIndex >= 0; oldIndex--)
        {
            /* Ignore zero terms immediately: */
            if (acb_poly_get_coeff_ptr(ODE->series,oldIndex) == NULL)
                continue;
            if (acb_is_zero(acb_poly_get_coeff_ptr(ODE->series,oldIndex)))
                continue;
            flint_printf("%w, ",oldIndex);
            num_of_nonzero++;
            acb_zero(temp);
            /* Loop through the polynomials */
            for (slong polyIndex = order(ODE); polyIndex >= 0; polyIndex--)
            {
                /* No more than degree(ODE) terms can contribute: */
                if (polyIndex > degree(ODE) + oldIndex - newIndex)
                {
                    polyIndex = degree(ODE) + oldIndex - newIndex + 1;
                    continue;
                }
                if (polyIndex > oldIndex)
                    polyIndex = oldIndex;
                if (polyIndex + newIndex - oldIndex < 0)
                    break;
                acb_add(temp,temp,diff_eq_coeff(ODE,polyIndex,polyIndex+newIndex-oldIndex),bits);
                if (polyIndex != 0)
                    acb_mul_si(temp,temp,oldIndex - polyIndex + 1,bits);
            }
            acb_addmul(new_coeff,acb_poly_get_coeff_ptr(ODE->series,oldIndex),temp,bits);
        }
        flint_printf(" to find %w.\n",newIndex);
        if (num_of_nonzero == 0)
        {
            printf("Zero!");
            break;
        }
        /* Divide by the appropriate coefficients */
        acb_set_ui(temp,newIndex+1);
        acb_rising_ui(temp,temp,order(ODE),bits);
        acb_mul(temp,temp,diff_eq_coeff(ODE,order(ODE),0),bits);
        acb_neg(temp,temp);

        acb_div(new_coeff,new_coeff,temp,bits);
        acb_poly_set_coeff_acb(ODE->series,newIndex+order(ODE),new_coeff);

        /* Compute the whole summand as a_b*z^b where b = newIndex + order(ODE) */
        if (in != NULL)
        {
            acb_pow_ui(temp,in,newIndex+order(ODE),bits);
            acb_mul(temp,new_coeff,temp,bits);
        }
        else
            acb_set(temp,new_coeff);
        /* If necessary, do the summation */
        if (out != NULL)
            acb_add(result,result,temp,bits);

        if (++newIndex >= 5000)
        {
            flint_printf("5000 summands have been computed and no convergence was achieved. Aborting.\n");
            newIndex = 0;
            break;
        }
        if (acb_is_zero(temp)) /* Due to unprecise imaginary parts, this might not catch all zeroes. Try (x+2)y'' - y' - y = 0 */
            continue;
        acb_get_abs_ubound_arf(err,temp,bits);
    } while (arf_cmpabs_2exp_si(err,-bits) >= 0);

    if (out != NULL)
        acb_set(out,result);

    acb_clear(new_coeff);
    acb_clear(result); acb_clear(temp);
    arf_clear(err);
    return newIndex;
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

void entry_point (ulong n, slong digits, slong z_val) {
    digits = digits*3.32193; /* Multiply by log_2(10) to obtain the number of binary digits */
    slong steps = 128;
    if (n == 0 || n >= 100) /* If n >= 100, that must be an error! */
        return;

    /* Polynomials */
    acb_poly_struct **polys = flint_malloc((n+1)*sizeof(acb_poly_t));
    for (ulong i = 0; i <= n; i++) polys[i] = NULL;

    acb_poly_t pol0,pol1,pol2;
    acb_poly_init(pol0); acb_poly_init(pol1); acb_poly_init(pol2);
    polys[0] = pol0; polys[1] = pol1; polys[2] = pol2;
    acb_poly_set_coeff_si(polys[2],1,1);
    acb_poly_set_coeff_si(polys[2],0,2);
    acb_poly_set_coeff_si(polys[1],0,-1);
    acb_poly_set_coeff_si(polys[0],0,-1);

    /* Input/Output */
    acb_t res, z;
    acb_init(res); acb_init(z);
    acb_set_si(z,z_val);
    acb_mat_t monodromy;

    /* Path */
    acb_ptr path = _acb_vec_init(steps);
    for (slong i = 0; i < steps; i++) {
        acb_set_si(path+i,i);
        acb_div_ui(path+i,path+i,steps/2,digits);
        acb_exp_pi_i(path+i,path+i,steps/2*digits);
    }

    acb_ode_t ODE = acb_ode_init(polys,NULL,z,n,digits);
    //~ if (ODE != NULL)
    //~ {
        //~ find_monodromy_matrix(monodromy,ODE,path,steps,digits);
        //~ if (monodromy == NULL)
            //~ flint_printf("Could not allocate memory. Please try again.\n");
        //~ else
        //~ {
            //~ acb_mat_printd(monodromy,7);
            //~ acb_mat_clear(monodromy);
        //~ }
    //~ }
    acb_poly_set_coeff_si(ODE->series,1,1);
    acb_poly_set_coeff_si(ODE->series,0,1);
    find_power_series_regular(res,ODE,z,digits);
    checkODE(polys,ODE,digits);

    for (ulong i = 0; i <= n; i++)
        if (polys[i] != NULL)
            acb_poly_clear(polys[i]);
    flint_free(polys);
    acb_ode_clear(ODE);
    acb_clear(res); acb_clear(z);
    _acb_vec_clear(path,steps);
    flint_cleanup();
    return;
}

int checkODE (acb_poly_struct **polys, acb_ode_t ODE, slong digits) {
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
    return 1;
}

int main (int argc, char **argv) {
    ulong n =       (argc >= 2) ? atol(argv[1]) : 1;
    slong digits =  (argc >= 3) ? atol(argv[2]) : 100;
    slong z =       (argc >= 4) ? atol(argv[3]) : 1;

    entry_point(n,digits,z);
    return 0;
}
