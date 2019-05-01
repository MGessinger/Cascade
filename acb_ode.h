#include <acb.h>
#include <acb_poly.h>
#include <acb_calc.h>

#ifndef ACB_ODE_H
#define ACB_ODE_H

#define INVALID_DATA (-1)
#define ORDINARY (0)
#define NON_CONVERGENT (1)
#define SINGULAR (2)

#define degree(ODE) ((ODE)->degree)
#define order(ODE) ((ODE)->order)
#define diff_eq_coeff(ODE,i,j) \
    ((j) <= (degree(ODE)) ? (ODE)->polys[i] + (j) : NULL)

typedef struct acb_ode_struct {
    slong order;
    slong degree;
    acb_struct **polys;
    acb_poly_t series;
} acb_ode_struct;

typedef acb_ode_struct* acb_ode_t;

short precondition(acb_poly_struct **polys, acb_ode_t ODE, acb_t z, slong prec);

acb_ode_t acb_ode_init (acb_poly_struct **polys, acb_poly_t initial, acb_t z, slong order, slong prec);

void acb_ode_clear(acb_ode_t ODE);

acb_ode_t acb_ode_copy(acb_ode_t ODE_out, acb_ode_t ODE_in);

#endif
