#ifndef ACB_ODE_H
#define ACB_ODE_H

#include <acb.h>
#include <acb_poly.h>

#define INVALID_DATA (-1)
#define NON_CONVERGENT (0)
#define SINGULAR (1)
#define ORDINARY (2)

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

short precondition (acb_poly_t *polys, acb_ode_t ODE);

void radiusOfConvergence(acb_ode_t ODE, arf_t radOfConv, slong bits);

acb_ode_t acb_ode_init (acb_poly_t *polys, acb_poly_t initial, slong order);

void acb_ode_clear (acb_ode_t ODE);

acb_ode_t acb_ode_set (acb_ode_t ODE_out, acb_ode_t ODE_in);

void acb_ode_shift (acb_ode_t ODE, acb_t a, slong bits);

acb_poly_t* acb_ode_fread(ulong *numberOfPols, const char *fileName, ulong maxOrder, slong bits);

void parsePoly(acb_poly_t polyOut, const char *polyString, slong strLength, slong bits);

slong acb_ode_reduce (acb_ode_t ODE);

#endif
