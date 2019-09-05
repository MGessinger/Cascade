#include "acb_ode.h"

/* Static function */

static short interpret (acb_poly_t *polys, acb_ode_t ODE)
{
    /* Find the characteristic "size" of the ODE */
    slong polyMaxDeg = 0, ord = order(ODE);
    while (polys[ord] == NULL || acb_poly_is_zero(polys[ord]))
        ord--;
    if (ord <= 0)
        return INVALID_DATA;
    for (slong i = 0; i <= ord; i++)
    {
        if (polys[i] == NULL)
            continue;
        if (polyMaxDeg < acb_poly_degree(polys[i]))
            polyMaxDeg = acb_poly_degree(polys[i]);
    }
    order(ODE) = ord;
    degree(ODE) = polyMaxDeg;
    return ORDINARY;
}

/* Setup and memory management*/

acb_ode_t acb_ode_init_blank(slong degree, slong order)
{
    if (degree < 0 || order <= 0)
        return NULL;
    /* Prepare the Differential equation for later use */
    acb_ode_t ODE = flint_malloc(sizeof(acb_ode_struct));
    if (ODE == NULL) {
        flint_printf("Initalisation of the differential equation failed. Please try again.\n");
        return NULL;
    }
    order(ODE) = order;
    degree(ODE) = degree;
    ODE->polys = flint_malloc((order(ODE)+1)*sizeof(acb_ptr));
    for (slong i = 0; i <= order(ODE); i++) {
        (ODE->polys)[i] = _acb_vec_init(degree(ODE)+1);
    }
    acb_poly_init(ODE->solution);
    return ODE;
}

acb_ode_t acb_ode_init (acb_poly_t *polys, acb_poly_t initial, slong order)
{
    /* Prepare the Differential equation for later use */
    if (polys == NULL)
        return NULL;
    acb_ode_t ODE = flint_malloc(sizeof(acb_ode_struct));
    if (ODE == NULL)
    {
        flint_printf("Initalisation of the differential equation failed. Please try again.\n");
        return NULL;
    }
    order(ODE) = order;
    if (interpret(polys,ODE) != ORDINARY)
    {
        flint_free(ODE);
        return NULL;
    }
    ODE->polys = flint_malloc((order(ODE)+1)*sizeof(acb_ptr));
    for (slong i = 0; i <= order(ODE); i++)
    {
        diff_eq_poly(ODE,i) = _acb_vec_init(degree(ODE)+1);
        if (polys[i] == NULL)
            continue;
        for (slong j = 0; j <= acb_poly_degree(polys[i]); j++)
        {
            acb_poly_get_coeff_acb(diff_eq_coeff(ODE,i,j),polys[i],j);
        }
    }
    acb_poly_init(ODE->solution);
    if (initial != NULL)
        acb_poly_set(ODE->solution,initial);
    return ODE;
}

void acb_ode_clear (acb_ode_t ODE)
{
    /* Free memory allocated for ODE*/
    if (ODE == NULL) return;
    for (int i = 0; i <= order(ODE); i++)
    {
        _acb_vec_clear((ODE->polys)[i],degree(ODE)+1);
    }
    flint_free(ODE->polys);
    acb_poly_clear(ODE->solution);
    flint_free(ODE);
    flint_cleanup();
    return;
}

acb_ode_t acb_ode_set (acb_ode_t ODE_out, acb_ode_t ODE_in)
{
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
        acb_poly_init(ODE_out->solution);
    }
    for (slong i = 0; i <= order(ODE_out); i++)
    {
        _acb_vec_set(diff_eq_poly(ODE_out,i),diff_eq_poly(ODE_in,i),degree(ODE_out)+1);
    }
    acb_poly_set(ODE_out->solution,ODE_in->solution);
    return ODE_out;
}

/* I/O */

acb_ode_t acb_ode_fread(const char *fileName, ulong maxOrder, slong bits)
{
    /* Read an ODE from a file */
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
    ulong numberOfPols;
    if (fscanf(input,"%*c%lu*(",&derivative) == 0)
    {
        flint_printf("The file format is wrong. Please make sure to declare the degree of derivation first.\n");
        fclose(input);
        return NULL;
    }
    numberOfPols = derivative;
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
    for (ulong i = 0; i <= numberOfPols; i++)
        acb_poly_init(polys[i]);
    do {
        if (fscanf(input,"%[^)*]%n",poly,&length) != 0)
            parsePoly(polys[derivative],poly,length,bits);
    } while (fscanf(input,"%*[^a-z]%*c%lu*(",&derivative) != EOF);
    fclose(input);
    acb_ode_t ODE =  acb_ode_init(polys,NULL,numberOfPols);
    for (ulong i = 0; i <= numberOfPols; i++)
        acb_poly_clear(polys[i]);
    flint_free(polys);
    flint_cleanup();
    return ODE;
}

void acb_ode_dump(acb_ode_t ODE, char *file)
{
    /* Dumps the ODE to file. If file is NULL, dump to stdout */
    if (ODE == NULL)
        return;
    FILE *out = stdout;
    if (file != NULL)
        out = fopen(file,"w");
    if (out == NULL)
        return;
    flint_fprintf(out,"Order: %w\nDegree: %w\n",order(ODE),degree(ODE));
    for (slong i = 0; i <= order(ODE); i++)
    {
        flint_fprintf(out,"diff_eq_poly(ODE,%w) = ",i);
        for (slong j = 0; j <= degree(ODE); j++)
        {
            acb_fprintd(out,diff_eq_coeff(ODE,i,j),20);
            flint_fprintf(out,"\t");
        }
        flint_fprintf(out,"\n");
    }
    flint_fprintf(out,"\nSolution:\n");
    acb_poly_fprintd(out,ODE->solution,10);
    flint_fprintf(out,"\n");
    if (out != stdout)
        fclose(out);
    return;
}

/* Transformations */

void acb_ode_shift (acb_ode_t ODE_out, acb_ode_t ODE_in, acb_srcptr a, slong bits)
{
    /* Shifts the origin to a */
    if (!ODE_in || !ODE_out)
        return;
    acb_ode_set(ODE_out,ODE_in);
    if (acb_is_zero(a))
        return;
    if (degree(ODE_out) == 0)
        return;
    for (slong j = 0; j <= order(ODE_out); j++)
    {
        if (diff_eq_poly(ODE_out,j) == NULL)
            continue;
        _acb_poly_taylor_shift(diff_eq_poly(ODE_out,j),a,degree(ODE_out)+1,bits);
    }
    return;
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
    degree(ODE) -= reduced;
    return reduced;
}

/* This function is so terrible, you probably don't want to use it. Ever. */

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
    acb_clear(coeff);
    return;
}
