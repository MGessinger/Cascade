#ifndef ACB_ODE_H
#define ACB_ODE_H

#include <acb.h>
#include <acb_poly.h>

#define INVALID_DATA (-1)
#define ORDINARY (1)

#define degree(ODE) ((ODE)->degree)
#define order(ODE) ((ODE)->order)
#define diff_eq_coeff(ODE,i,j) ((ODE)->polys[i] + (j))
#define diff_eq_poly(ODE,i) ((ODE)->polys[i])

typedef struct acb_ode_struct {
    slong order;
    slong degree;
    acb_struct **polys;
    acb_poly_t solution;
} acb_ode_struct;

typedef acb_ode_struct* acb_ode_t;

/* Predefined Equations (in examples.c) */
acb_ode_t acb_ode_legendre(ulong n);
acb_ode_t acb_ode_bessel(acb_t nu, slong bits);
acb_ode_t acb_ode_hypgeom(acb_t a, acb_t b, acb_t c, slong bits);

/* Setup and memory management */

acb_ode_t acb_ode_init (acb_poly_t *polys, acb_poly_t initial, slong order);
void acb_ode_clear (acb_ode_t ODE);
acb_ode_t acb_ode_set (acb_ode_t ODE_out, acb_ode_t ODE_in);

/* I/O */

acb_ode_t acb_ode_fread(const char *fileName, ulong maxOrder, slong bits);
void acb_ode_dump(acb_ode_t ODE, char *file);

/* Transformations */

void acb_ode_shift (acb_ode_t ODE, acb_t a, slong bits);
slong acb_ode_reduce (acb_ode_t ODE);

/* Non-documented function (because terrible) */

void parsePoly(acb_poly_t polyOut, const char *polyString, slong strLength, slong bits);

#endif
