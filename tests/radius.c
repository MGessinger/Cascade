#include "cascade.h"
#include <acb.h>
#include <flint/flint.h>

int main ()
{
	int return_value = EXIT_SUCCESS;
	int close;
	slong prec, n;

	arb_t rad, temp;
	acb_ode_t ode;
	flint_rand_t state;

	arb_init(rad);
	arb_init(temp);
	acb_ode_init_blank(ode, 15, 1);
	flint_randinit(state);

	for (slong iter = 0; iter < 100; iter++)
	{
		prec = 30 + n_randint(state, 128);
		n = 1 + n_randint(state, 14);

		_acb_vec_zero(acb_ode_poly(ode, 1), degree(ode) + 1);

		for (slong i = 0; i < n; i++)
			acb_randtest(acb_ode_coeff(ode, 0, i), state, prec, 15);
		_acb_poly_product_roots(acb_ode_poly(ode, 1), acb_ode_poly(ode, 0), n, prec);

		radius_of_convergence(rad, ode, 20, prec);

		for (slong i = 0; i < n; i++)
		{
			acb_abs(temp, acb_ode_coeff(ode, 0, i), prec);
			if (arb_contains_zero(temp))
				continue;
			close = (arb_lt(rad, temp) || arb_overlaps(rad, temp));
			if (!close)
			{
				return_value = EXIT_FAILURE;
				break;
			}
		}
		if (return_value == EXIT_FAILURE)
			break;
	}

	arb_clear(rad);
	arb_clear(temp);
	acb_ode_clear(ode);
	flint_randclear(state);
	flint_cleanup();
	return return_value;
}
