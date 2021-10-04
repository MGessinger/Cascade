#include "acb_ode.h"
#include <acb.h>
#include <flint/flint.h>

int main ()
{
	int return_value = EXIT_SUCCESS;
	slong prec, n;

	acb_ode_t ode;
	flint_rand_t state;

	flint_randinit(state);

	for (slong iter = 0; iter < 100; iter++)
	{
		prec = 30 + n_randint(state, 128);

		acb_ode_random(ode, state, prec);

		n = n_randint(state, degree(ode));
		for (slong i = 0; i <= order(ode); i++)
			for (slong j = 0; j < n; j++)
				acb_zero(acb_ode_coeff(ode, i, j));

		int red = acb_ode_reduce(ode);

		acb_ode_clear(ode);

		if (red < n)
		{
			return_value = EXIT_FAILURE;
			break;
		}
	}

	flint_randclear(state);
	flint_cleanup();
	return return_value;
}
