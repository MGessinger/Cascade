#include "cascade.h"
#include <acb_poly.h>
#include <flint/flint.h>

int main ()
{
	int return_value = 0; /* 0 means success. */
	acb_ptr vec = _acb_vec_init(5);
	acb_set_si(vec+0, 1);
	acb_set_si(vec+1, 2);
	acb_set_si(vec+2, 3);
	acb_set_si(vec+3, 4);
	acb_set_si(vec+4, 5);

	acb_poly_graeffe_transform(vec, vec, 5, 1024);

	/* Check the results */
	acb_t expected;
	acb_init(expected);
	acb_set_si(expected, 1);
	if (!acb_contains(vec+0, expected))
		return_value |= 1 << 0;
	acb_set_si(expected, 2);
	if (!acb_contains(vec+1, expected))
		return_value |= 1 << 1;
	acb_set_si(expected, 3);
	if (!acb_contains(vec+2, expected))
		return_value |= 1 << 2;
	acb_set_si(expected, 14);
	if (!acb_contains(vec+3, expected))
		return_value |= 1 << 3;
	acb_set_si(expected, 25);
	if (!acb_contains(vec+4, expected))
		return_value |= 1 << 4;
	acb_clear(expected);
	_acb_vec_clear(vec, 5);
	flint_cleanup();
	return return_value;
}

// Expected output-coefficients:
// 1 2 3 14 25
