#ifndef EXAMPLES_H_
#define EXAMPLES_H_

#include "acb_ode.h"

acb_ode_t acb_ode_legendre(ulong n);

acb_ode_t acb_ode_bessel(acb_t nu, slong bits);

acb_ode_t acb_ode_hypgeom(acb_t a, acb_t b, acb_t c, slong bits);

#endif /* EXAMPLES _H_ */
