#ifndef ACB_ODE_H
#define ACB_ODE_H

#include <acb.h>
#include <acb_poly.h>

#define degree(ODE) ((ODE)->degree)
#define order(ODE) ((ODE)->order)
#define acb_ode_poly(ODE, i) (((ODE)->polys) + (i)*(degree(ODE)+1))
#define acb_ode_coeff(ODE, i, j) (acb_ode_poly(ODE, i) + (j))

typedef struct acb_ode_struct {
	slong order;
	slong degree;
	slong alloc;
	acb_ptr polys;
} acb_ode_struct;

typedef	acb_ode_struct acb_ode_t[1];

/* Setup and memory management */
void	acb_ode_init_blank (acb_ode_t ODE, slong degree, slong order);
void	acb_ode_init (acb_ode_t ODE, acb_poly_t *polys, slong order);
void	acb_ode_clear (acb_ode_t ODE);
void	acb_ode_set (acb_ode_t ODE_out, acb_ode_t ODE_in);

/* I/O */
void	acb_ode_dump (acb_ode_t ODE, char *file);

/* Transformations */
void	acb_ode_shift (acb_ode_t ODE_out, acb_ode_t ODE_in, acb_srcptr a, slong bits);
slong	acb_ode_reduce (acb_ode_t ODE);
slong	acb_ode_valuation (acb_ode_t ODE);

/* Predefined Equations (in examples.c) */
void	acb_ode_legendre (acb_ode_t ODE, ulong n);
void	acb_ode_bessel (acb_ode_t ODE, acb_t nu, slong bits);
void	acb_ode_hypgeom (acb_ode_t ODE, acb_t a, acb_t b, acb_t c, slong bits);

#endif
