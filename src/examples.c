#include "acb_ode.h"

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

	acb_ode_t ODE = acb_ode_init(polys,2);

	for (int i = 0; i < 3; i++)
		acb_poly_clear(polys[i]);
	flint_free(polys);
	return ODE;
}

acb_ode_t acb_ode_bessel(acb_t nu, slong bits)
{
	acb_poly_t* polys = flint_malloc(3*sizeof(acb_poly_t));
	if (polys == NULL)
		return NULL;
	for (int i = 0; i < 3; i++)
		acb_poly_init(polys[i]);
	acb_poly_set_coeff_si(polys[2],2,1);
	acb_poly_set_coeff_si(polys[1],1,1);
	acb_poly_set_coeff_si(polys[0],2,1);
	acb_mul(nu,nu,nu,bits);
	acb_neg(acb_poly_get_coeff_ptr(polys[0],0),nu);

	acb_ode_t ODE = acb_ode_init(polys,2);

	for (int i = 0; i < 3; i++)
		acb_poly_clear(polys[i]);
	flint_free(polys);
	return ODE;
}

acb_ode_t acb_ode_hypgeom(acb_t a, acb_t b, acb_t c, slong bits)
{
	acb_t temp;
	acb_init(temp);
	acb_poly_t* polys = flint_malloc(3*sizeof(acb_poly_t));
	if (polys == NULL)
		return NULL;
	for (int i = 0; i < 3; i++)
		acb_poly_init(polys[i]);
	/* z*(1-z) */
	acb_poly_set_coeff_si(polys[2],1,1);
	acb_poly_set_coeff_si(polys[2],2,-1);
	/* -ab */
	acb_mul(temp,a,b,bits);
	acb_neg(temp,temp);
	acb_poly_set_coeff_acb(polys[0],0,temp);
	/* c - (a+b+1)z */
	acb_poly_set_coeff_acb(polys[1],0,c);
	acb_one(temp);
	acb_add(temp,temp,a,bits);
	acb_add(temp,temp,b,bits);
	acb_neg(temp,temp);
	acb_poly_set_coeff_acb(polys[1],1,temp);

	acb_ode_t ODE = acb_ode_init(polys,2);

	for (int i = 0; i < 3; i++)
		acb_poly_clear(polys[i]);
	flint_free(polys);
	acb_clear(temp);
	return ODE;
}
