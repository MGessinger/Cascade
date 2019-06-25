#include "juliaInterface.h"

/*ulong set_tolerance(ulong new_tol)
{
    return convergence_tolerance = new_tol;
}*/

acb_ode_t acb_ode_setup_blank(slong degree, slong order)
{
    if (degree < 0 || order <= 0)
        return NULL;
    /* Prepare the Differential equation for later use */
    acb_ode_t ODE = flint_malloc(sizeof(acb_ode_struct));
    if (ODE == NULL) {
        flint_printf("Initalisation of the differential equation failed. Please try again.\n");
        return NULL;
    }
    order(ODE) = order;
    degree(ODE) = degree;
    ODE->polys = flint_malloc((order(ODE)+1)*sizeof(acb_ptr));
    for (slong i = 0; i <= order(ODE); i++) {
        (ODE->polys)[i] = _acb_vec_init(degree(ODE)+1);
    }
    acb_poly_init(ODE->solution);
    return ODE;
}

short acb_ode_set_poly(acb_ode_t ODE, acb_poly_struct poly, slong index)
{
    if (ODE == NULL)
        return INVALID_DATA;
    if (index > order(ODE))
        return INVALID_DATA;
    if (acb_poly_degree(&poly) > degree(ODE))
        return INVALID_DATA;
    for (slong j = 0; j <= acb_poly_degree(&poly); j++) {
        acb_poly_get_coeff_acb(diff_eq_coeff(ODE,index,j),&poly,j);
    }
    return ORDINARY;
}

void acb_set_initial(acb_ode_t ODE, acb_poly_struct init)
{
    if (ODE == NULL)
        return;
    acb_poly_set(ODE->solution,&init);
    return;
}
    
acb_poly_struct find_power_series_julia(acb_ode_t ODE, acb_struct in, slong bits)
{
    if (acb_poly_is_zero(ODE->solution))
    {
        flint_printf("No initial values have been specified. Please run setInitalValues first.\n");
        return ODE->solution[0];
    }
    find_power_series(ODE, &in, bits);
    flint_cleanup();
    return ODE->solution[0];
}

acb_mat_struct find_monodromy_julia (acb_ode_t ODE, slong bits)
{
    acb_mat_t monodromy;
    if (ODE == NULL)
    {
        printf("The ODE is the NULL-pointer. Please confirm input.\n");
        return monodromy[0];
    }
    slong steps = 16;
    acb_ptr path = _acb_vec_init(steps+1);
    
    for (slong i = 0; i < steps; i++)
    {
        acb_set_si(path+i,i);
        acb_div_ui(path+i,path+i,steps/2,bits);
        acb_exp_pi_i(path+i,path+i,bits);
    }
    acb_set(path+steps,path);
    _acb_vec_scalar_div_ui(path,path,steps+1,128,steps/2*bits);
    find_monodromy_matrix(monodromy,ODE,path,steps+1,bits);
    _acb_vec_clear(path,steps+1);
    flint_cleanup();
    return monodromy[0];
}
