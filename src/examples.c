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
    fmpz_t Int;
    acb_t coeff;
    acb_poly_init(initial);
    fmpz_init(Int);
    acb_init(coeff);

    fmpz_bin_uiui(Int,n,(n+(n%2))/2);
    if (n % 2 == 1)
        fmpz_mul_ui(Int,Int,n+1);
    acb_set_fmpz(coeff,Int);
    acb_mul_2exp_si(coeff,coeff,-n);
    acb_poly_set_coeff_acb(initial,n%2,coeff);

    acb_ode_t ODE = acb_ode_init(polys,initial,2);

    for (int i = 0; i < 3; i++)
        acb_poly_clear(polys[i]);
    flint_free(polys);
    acb_poly_clear(initial);
    acb_clear(coeff);
    fmpz_clear(Int);
    return ODE;
}

acb_ode_t acb_ode_bessel(acb_struct nu, slong bits)
{
    acb_printd(&nu,10);
    acb_poly_t* polys = flint_malloc(3*sizeof(acb_poly_t));
    if (polys == NULL)
        return NULL;
    for (int i = 0; i < 3; i++)
        acb_poly_init(polys[i]);
    acb_poly_set_coeff_si(polys[2],2,1);
    acb_poly_set_coeff_si(polys[1],1,1);
    acb_poly_set_coeff_si(polys[0],2,1);
    acb_mul(&nu,&nu,&nu,bits);
    acb_neg(acb_poly_get_coeff_ptr(polys[0],0),&nu);

    acb_poly_t initial;
    acb_poly_init(initial);
    if (acb_is_zero(&nu))
        acb_poly_one(initial);
    else
        acb_poly_set_coeff_si(initial,1,1);

    acb_ode_t ODE = acb_ode_init(polys,initial,2);

    for (int i = 0; i < 3; i++)
        acb_poly_clear(polys[i]);
    flint_free(polys);
    acb_poly_clear(initial);
    return ODE;
}
