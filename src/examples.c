#include "acb_ode.h"

void acb_ode_legendre (acb_ode_t ODE, ulong n)
{
	acb_ode_init_blank(ODE, 2, 2);

	acb_set_si(diff_eq_coeff(ODE, 0, 0), n*(n+1));
	acb_set_si(diff_eq_coeff(ODE, 1, 1), -2);
	acb_set_si(diff_eq_coeff(ODE, 2, 0), 1);
	acb_set_si(diff_eq_coeff(ODE, 2, 2), -1);
}

void acb_ode_bessel (acb_ode_t ODE, acb_t nu, slong bits)
{
	acb_ode_init_blank(ODE, 2, 2);

	acb_set_si(diff_eq_coeff(ODE, 2, 2), 1);
	acb_set_si(diff_eq_coeff(ODE, 1, 1), 1);
	acb_set_si(diff_eq_coeff(ODE, 0, 2), 1);
	acb_mul(nu, nu, nu, bits);
	acb_neg(diff_eq_coeff(ODE, 0, 0), nu);
}

void acb_ode_hypgeom (acb_ode_t ODE, acb_t a, acb_t b, acb_t c, slong bits)
{
	acb_ode_init_blank(ODE, 2,2);

	acb_t temp;
	acb_init(temp);
	/* z*(1-z) */
	acb_set_si(diff_eq_coeff(ODE, 2, 1), 1);
	acb_set_si(diff_eq_coeff(ODE, 2, 2), -1);
	/* c - (a+b+1)z */
	acb_set(diff_eq_coeff(ODE, 1, 0), c);
	acb_one(temp);
	acb_add(temp, temp, a, bits);
	acb_add(temp, temp, b, bits);
	acb_neg(diff_eq_coeff(ODE, 1, 1), temp);
	/* -ab */
	acb_mul(temp, a, b, bits);
	acb_neg(diff_eq_coeff(ODE, 0, 0), temp);

	acb_clear(temp);
}
