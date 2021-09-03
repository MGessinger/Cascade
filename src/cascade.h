#ifndef CASCADE_H_
#define CASCADE_H_

#include <acb_poly.h>
#include <acb_mat.h>
#include "acb_ode.h"

void	radius_of_convergence (arb_t rad_of_conv, acb_ode_t ODE, slong n, slong bits);
slong	truncation_order (arb_t eta, arb_t alpha, slong bits);

/* ============================== Fuchs Solver ============================== */

void	acb_ode_solve_fuchs (acb_poly_t res, acb_ode_t ODE, slong deg, slong bits);

/* Compute analytic continuation and monodromy */
void	analytic_continuation (acb_poly_t res, acb_ode_t ODE, acb_srcptr path,
		slong len, slong deg, slong bits);
void	find_monodromy_matrix (acb_mat_t mono, acb_ode_t ODE, slong bits);

/* ============================ Frobenius Solver ============================ */

void	indicial_polynomial (acb_poly_t result, acb_ode_t ODE, slong nu, slong shift, slong prec);
void	indicial_polynomial_evaluate (acb_t result, acb_ode_t ODE, slong nu, acb_t rho, slong shift, slong prec);

void	_acb_ode_solve_frobenius (acb_poly_t res, acb_ode_t ODE, acb_t rho, slong sol_degree, slong prec);
void	acb_ode_solve_frobenius (acb_ode_solution_t sol, acb_ode_t ODE, slong sol_degree, slong prec);

/* Inlines */

static inline slong clamp (slong in, slong min, slong max)
{
	if (in < min)
		return min;
	else if (in > max)
		return max;
	else
		return in;
}

#endif /* CASCADE_H_ */
