#include "acb_ode.h"

short precondition (acb_poly_struct **polys, acb_ode_t ODE, acb_t z, slong prec) {
    /* Error Checking */
    if (polys == NULL) {
        flint_printf("Please provide defining polynomials.\n");
        return INVALID_DATA;
    }
    /* Find the highest degree amongst the given polynomials */
    slong polyMaxLength = 0, ord = ODE->order;
    for (slong i = ord; i >= 0; i--)
    {
        if (polys[i] == NULL)
        {
            if (polyMaxLength == 0) ord--;
        }
        else if (polyMaxLength < acb_poly_length(polys[i]))
            polyMaxLength = acb_poly_length(polys[i]);
    }
    order(ODE) = ord;
    degree(ODE) = polyMaxLength;
    if (ord <= 0)
    {
        flint_printf("The order of the differential equation has to be positive.\n");
        return INVALID_DATA;
    }
    flint_printf("The polynomials have degree at most %wd .\n",polyMaxLength-1);
    if (z == NULL)
        return ORDINARY;

    /* Check if z is a singular point (regularity is elsewhere) */
    acb_t r; acb_init(r);
    acb_t m; acb_init(m);
    acb_poly_evaluate(r,polys[ord],m,prec);
    if (acb_is_zero(r))
    {
        flint_printf("Warning. This is a singular point. Proceed? (y/n)\n");
        if (getchar() != 'y')
            return SINGULAR;
    }

    /* Check if the series expansion converges */
    acb_ptr polyder = _acb_vec_init(acb_poly_degree(polys[ord]));
    _acb_poly_derivative(polyder,polys[ord]->coeffs,acb_poly_degree(polys[ord])+1,prec);
    _acb_poly_root_inclusion(r, m, acb_poly_get_coeff_ptr(polys[ord],0), polyder, acb_poly_degree(polys[ord])+1, prec);
    int convergent = ORDINARY;
    if (!acb_is_exact(r) && !acb_contains(r,z))
    {
        convergent = NON_CONVERGENT;
        flint_printf("There is definitely a root in the way. We have an enclosure of \n");
        acb_printn(r,prec,0); flint_printf(" and we are evaluating at "); acb_printn(z,5,0); flint_printf("\n");
    }
    _acb_vec_clear(polyder,acb_poly_degree(polys[ord]));
    acb_clear(r);
    acb_clear(m);
    return convergent;
}

acb_ode_t acb_ode_init (acb_poly_struct **polys, acb_poly_t initial, acb_t z, slong order, slong prec) {
    /* Prepare the Differential equation for later use */
    acb_ode_t ODE = flint_malloc(sizeof(acb_ode_struct));
    if (ODE == NULL) {
        flint_printf("Initalisation of the differential equation failed. Please try again.\n");
        return NULL;
    }
    order(ODE) = order;
    int ode_case = precondition(polys,ODE,z,prec);
    if (ode_case != ORDINARY) {
        flint_free(ODE);
        return NULL;
    }
    ODE->polys = flint_malloc((order(ODE)+1)*sizeof(acb_ptr));
    for (slong i = 0; i <= order(ODE); i++) {
        (ODE->polys)[i] = _acb_vec_init(degree(ODE)+1);
        if (polys[i] == NULL) continue;
        for (slong j = 0; j <= acb_poly_degree(polys[i]); j++) {
            acb_poly_get_coeff_acb(diff_eq_coeff(ODE,i,j),polys[i],j);
        }
    }
    acb_poly_init(ODE->series);
    if (initial != NULL)
        acb_poly_set(ODE->series,initial);
    return ODE;
}

void acb_ode_clear (acb_ode_t ODE) {
    if (ODE == NULL) return;
    for (int i = 0; i <= order(ODE); i++) {
        _acb_vec_clear((ODE->polys)[i],degree(ODE)+1);
    }
    flint_free(ODE->polys);
    acb_poly_clear(ODE->series);
    flint_free(ODE);
    return;
}

acb_ode_t acb_ode_copy (acb_ode_t ODE_out, acb_ode_t ODE_in) {
    /* Copy data from ODE_in to an existing ODE structure or create a new one */
    if (ODE_out == NULL) {
        ODE_out = flint_malloc(sizeof(acb_ode_struct));
        if (ODE_out == NULL) {
            flint_printf("Initalisation of the differential equation failed. Please try again.\n");
            return NULL;
        }
        order(ODE_out) = order(ODE_in);
        degree(ODE_out) = degree(ODE_in);
        ODE_out->polys = flint_malloc((order(ODE_out)+1)*sizeof(acb_ptr));
        for (slong i = 0; i <= order(ODE_out); i++) {
            (ODE_out->polys)[i] = _acb_vec_init(degree(ODE_out)+1);
        }
        acb_poly_init(ODE_out->series);
    }
    for (slong i = 0; i <= order(ODE_out); i++) {
        for (slong j = 0; j <= degree(ODE_out); j++) {
            acb_set(diff_eq_coeff(ODE_out,i,j),diff_eq_coeff(ODE_in,i,j));
        }
    }
    acb_poly_set(ODE_out->series,ODE_in->series);
    return ODE_out;
}

void acb_ode_shift (acb_ode_t ODE, acb_t a, slong bits)
{
    if (degree(ODE) >= 1)
        for (slong j = 0; j <= order(ODE); j++)
        {
            if (ODE->polys[j] == NULL)
                continue;
            _acb_poly_taylor_shift(ODE->polys[j],a,degree(ODE)+1,bits);
        }
    return;
}
