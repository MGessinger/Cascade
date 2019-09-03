#include "juliaInterface.h"

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
    arb_t rad, eta;
    arb_init(rad);
    arb_init(eta);
    acb_abs(eta,&in,bits);
    radiusOfConvergence(rad,ODE,20,bits);
    find_power_series(ODE, truncation_order(eta,rad,bits), bits);
    arb_clear(rad);
    arb_clear(eta);
    return ODE->solution[0];
}
