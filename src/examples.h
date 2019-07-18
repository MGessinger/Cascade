#ifndef EXAMPLES_H_
#define EXAMPLES_H_

#include "acb_ode.h"

acb_ode_t acb_ode_legendre(ulong n);

acb_ode_t acb_ode_bessel(acb_struct nu, slong bits);

acb_ode_t acb_ode_hypgeom(acb_struct a, acb_struct b, acb_struct c, slong bits);

#endif /* EXAMPLES _H_ */
