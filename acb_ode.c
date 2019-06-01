#include "acb_ode.h"

short precondition (acb_poly_t *polys, acb_ode_t ODE) {
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
    if (ord <= 0 || polyMaxLength <= 0)
    {
        flint_printf("The order of the differential equation has to be positive.\n");
        return INVALID_DATA;
    }
    flint_printf("The polynomials have degree at most %wd .\n",polyMaxLength-1);
    return ORDINARY;
}

acb_ode_t acb_ode_init (acb_poly_t *polys, acb_poly_t initial, slong order) {
    /* Prepare the Differential equation for later use */
    acb_ode_t ODE = flint_malloc(sizeof(acb_ode_struct));
    if (ODE == NULL) {
        flint_printf("Initalisation of the differential equation failed. Please try again.\n");
        return NULL;
    }
    order(ODE) = order;
    int ode_case = precondition(polys,ODE);
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

acb_ode_t acb_ode_set (acb_ode_t ODE_out, acb_ode_t ODE_in) {
    /* Copy data from ODE_in to an existing ODE structure or create a new one */
    if (ODE_out == NULL)
    {
        ODE_out = flint_malloc(sizeof(acb_ode_struct));
        if (ODE_out == NULL)
        {
            flint_printf("Initalisation of the differential equation failed. Please try again.\n");
            return NULL;
        }
        order(ODE_out) = order(ODE_in);
        degree(ODE_out) = degree(ODE_in);
        ODE_out->polys = flint_malloc((order(ODE_out)+1)*sizeof(acb_ptr));
        for (slong i = 0; i <= order(ODE_out); i++)
        {
            (ODE_out->polys)[i] = _acb_vec_init(degree(ODE_out)+1);
        }
        acb_poly_init(ODE_out->series);
    }
    for (slong i = 0; i <= order(ODE_out); i++)
    {
        _acb_vec_set(diff_eq_poly(ODE_out,i),diff_eq_poly(ODE_in,i),degree(ODE_out)+1);
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

acb_poly_t* acb_ode_fread(ulong *numberOfPols, const char *fileName, ulong maxOrder, slong bits)
{
    if (maxOrder == 0)
        maxOrder = UWORD_MAX;
    FILE *input = fopen(fileName,"r");
    if (input == NULL)
    {
        flint_printf("Could not open file %s. Please confirm input!\n",fileName);
        return NULL;
    }
    char poly[512];
    long unsigned derivative = 0;
    int length = 0;
    if (fscanf(input,"%*c%lu*(",&derivative) == 0)
    {
        flint_printf("The file format is wrong. Please make sure to declare the degree of derivation first.\n");
        fclose(input);
        return NULL;
    }
    *numberOfPols = derivative;
    if (derivative > maxOrder)
    {
        flint_printf("The order of the ODE was larger than allowed. Aborted.\n");
        fclose(input);
        return NULL;
    }
    acb_poly_t *polys = malloc((derivative+1)*sizeof(acb_poly_t));
    if (polys == NULL)
    {
        flint_printf("Could not allocate memory. Please try again.\n");
        fclose(input);
        return NULL;
    }
    for (ulong i = 0; i <= *numberOfPols; i++)
        acb_poly_init(polys[i]);
    do {
        if (fscanf(input,"%[^)*]%n",poly,&length) != 0)
            parsePoly(polys[derivative],poly,length,bits);
    } while (fscanf(input,"%*[^a-z]%*c%lu*(",&derivative) != EOF);
    fclose(input);
    return polys;
}

void parsePoly(acb_poly_t polyOut, const char *polyString, const slong strLength, slong bits)
{
    if (strLength == 0 || polyString[0] == '\0')
    {
        acb_poly_zero(polyOut);
        return;
    }
    if (polyOut == NULL)
    {
        flint_printf("The output polynomial is the NULL pointer. Please check your input.\n");
        return;
    }
    flint_printf("Parsing %s into a polynomial... ",polyString,strLength);
    acb_t coeff; acb_init(coeff);

    int lengthOfReal = 0, lengthOfImag;
    char realPart[128];
    char imagPart[128];
    slong totalLength = 0;
    slong index = 0;

    while (totalLength < strLength)
    {
        acb_zero(coeff);
        lengthOfReal = lengthOfImag = 0;
        if (sscanf(polyString+totalLength,"%[^, +]%n",realPart,&lengthOfReal) != 0)
        {
            arb_set_str(acb_realref(coeff),realPart,bits);
            totalLength += lengthOfReal;
            acb_poly_set_coeff_acb(polyOut,index,coeff);
        }
        if (totalLength >= strLength)
            break;
        if (sscanf(polyString+totalLength,"%*[ +]%[^j,]j%n",imagPart,&lengthOfImag) != 0)
        {
            if (lengthOfImag != 0)
            {
                arb_set_str(acb_imagref(coeff),imagPart,bits);
                totalLength += lengthOfImag;
                acb_poly_set_coeff_acb(polyOut,index,coeff);
            }
        }
        sscanf(polyString+totalLength,"%*[, +]%n",&lengthOfReal);
        totalLength += lengthOfReal;
        index++;
    }
    flint_printf("Done!\n");
    acb_clear(coeff);
    return;
}

acb_ode_t acb_ode_simplify(acb_ode_t ODE)
{
    if (ODE == NULL)
        return NULL;
    slong shift = 0;
    acb_ode_t ODEfixed = NULL;
    for (; shift <= order(ODE); shift++)
    {
        /* shift finds the number of leading Zero polynomials */
        if (!_acb_vec_is_zero(diff_eq_poly(ODE,shift),degree(ODE)+1))
        {
            break;
        }
    }
    if (shift != 0)
    {
        ODEfixed = acb_ode_set(NULL,ODE);
        order(ODEfixed) -= shift;
        for (slong i = 0; i <= order(ODEfixed); i++)
            _acb_vec_set(diff_eq_poly(ODEfixed,i),diff_eq_poly(ODE,i+shift),degree(ODE)+1);
        for (slong j = 1; j <= shift; j++)
        {
            _acb_vec_clear(diff_eq_poly(ODEfixed,j+order(ODEfixed)),degree(ODE)+1);
        }
    }
    return ODEfixed;
}

slong acb_ode_reduce (acb_ode_t ODE)
{
    /* Divides all polynomials by z^n, if they share such a factor */
    if (ODE == NULL)
        return 0;
    slong reduced = 0;
    while (acb_is_zero(diff_eq_coeff(ODE,order(ODE),reduced)))
    {
        reduced++;
    }
    for (slong i = 0; i < order(ODE); i++)
    {
        /* reduced finds the number of leading Zero coefficients */
        for (slong j = 0; j < reduced; j++)
            if(!acb_is_zero(diff_eq_coeff(ODE,i,j)))
            {
                reduced = j;
                break;
            }
    }
    if (reduced != 0)
    {
        for (slong i = 0; i<= order(ODE); i++)
            _acb_poly_shift_right(diff_eq_poly(ODE,i),diff_eq_poly(ODE,i),degree(ODE)+1,reduced);
    }
    return reduced;
}
