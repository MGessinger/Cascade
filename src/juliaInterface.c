#include "juliaInterface.h"

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

acb_poly_struct find_power_series_julia(acb_poly_t res, acb_ode_t ODE, acb_struct in, slong bits)
{
    acb_poly_t poly;
    arb_t rad, eta;
    arb_init(rad);
    arb_init(eta);
    acb_poly_init(poly);
    acb_poly_set(poly,res);
    acb_abs(eta,&in,bits);
    radiusOfConvergence(rad,ODE,20,bits);
    find_power_series(poly,ODE, truncation_order(eta,rad,bits), bits);
    arb_clear(rad);
    arb_clear(eta);
    return poly[0];
}
