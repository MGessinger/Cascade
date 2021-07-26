#include "cascade.h"
#include <acb_poly.h>

int main ()
{
	slong prec = 1024;
	acb_t z; acb_init(z);
	acb_zero(z);
	acb_ode_t ODE;
	acb_ode_bessel(ODE, z, prec);
	acb_ode_reduce(ODE);

	acb_poly_t poly; acb_poly_init(poly);
	acb_poly_one(poly);

	acb_ode_solve_fuchs(poly, ODE, 150, prec);
	acb_ode_clear(ODE);

	int return_value = 0;
	acb_t c; acb_init(c);
	acb_one(z);
	for (int i = 2; i < 150; i++) {
		acb_poly_get_coeff_acb(c, poly, i);
		if (i%2 != 0)
		{
			if (!acb_is_zero(c))
			{
				return_value = 1;
				break;
			}
		}
		else
		{
			acb_div_si(z, z, -i*i, prec);
			if (!acb_contains(c, z))
			{
				return_value = 1;
				break;
			}
		}
	}

	acb_clear(z);
	acb_clear(c);
	acb_poly_clear(poly);
	flint_cleanup();
	return return_value;
}
