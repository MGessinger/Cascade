#include "examples.h"

acb_ode_t acb_ode_legendre(ulong n)
{
    acb_poly_t* polys = flint_malloc(3*sizeof(acb_poly_t));
    if (polys == NULL)
        return NULL;
    for (int i = 0; i < 3; i++)
        acb_poly_init(polys[i]);
    acb_poly_set_coeff_si(polys[0],0,n*(n+1));
    acb_poly_set_coeff_si(polys[1],1,-2);
    acb_poly_set_coeff_si(polys[2],0,1);
    acb_poly_set_coeff_si(polys[2],2,-1);

    acb_poly_t initial;
    acb_poly_init(initial);
    acb_poly_one(initial);

    acb_ode_t ODE = acb_ode_init(polys,initial,2);

    for (int i = 0; i < 3; i++)
        acb_poly_clear(polys[i]);
    flint_free(polys);
    acb_poly_clear(initial);
    return ODE;
}
