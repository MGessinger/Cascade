#ifndef EXAMPLES_H_
#define EXAMPLES_H_

#include "acb_ode.h"

acb_ode_t acb_ode_legendre(ulong n);

acb_ode_t acb_ode_bessel(acb_struct nu, slong bits);

#endif /* EXAMPLES _H_ */
