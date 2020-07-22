#include "acb_ode.h"
#include <acb.h>
#include <flint/flint.h>

int main ()
{
	acb_t z;
	acb_init(z);
	acb_zero(z);

	acb_ode_t ode = acb_ode_bessel(z, 1024);
	slong red = acb_ode_reduce(ode);

	acb_ode_clear(ode);
	acb_clear(z);
	flint_cleanup();
	return (red - 1);
}
