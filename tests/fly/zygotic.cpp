/*
 * zygotic.cpp
 *
 *  Created on: Feb 5, 2013
 *      Author: zhlou
 */

#include "zygotic.h"

zygotic::zygotic(maternal& in_maternal, FILE *fp, char* parm_section) :
        TheMaternal(in_maternal), defs(in_maternal.getProblem())
{
    parm = ReadParameters(fp, parm_section);
    D = (double *) calloc(defs.ngenes, sizeof(double));
    vinput = (double *) calloc(defs.ngenes * defs.nnucs, sizeof(double));
    bot2 = (double *) calloc(defs.ngenes * defs.nnucs, sizeof(double));
    bot = (double *) calloc(defs.ngenes * defs.nnucs, sizeof(double));

}

/*** A FUNCTION THAT READS PARAMETERS FROM THE DATA FILE INTO STRUCTS ******/

/*** ReadParamters: reads the parameters for a simulation run from the *****
 *                  eqparms or input section of the data file as indicated *
 *                  by the section_title argument and does the conversion  *
 *                  of protein half lives into lambda parameters.          *
 ***************************************************************************/
EqParms zygotic::ReadParameters(FILE* fp, char* section_title)
{
    EqParms l_parm; /* local copy of EqParm struct */
    double *tempparm, *tempparm1; /* temporary array to read parms */

    int i; /* local loop counter */
    int c; /* used to parse text lines */
    int lead_punct; /* flag for leading punctuation */
    int linecount = 0; /* keep track of # of lines read */
    int Tcount = 0; /* keep track of T lines read */
    int Ecount = 0; /* keep track of E lines read */

    char *base; /* pointer to beginning of line string */
    char *record; /* string for reading whole line of params */

    char **fmt, **fmt1; /* array of format strings for reading params */
    char *skip, *skip1; /* string of values to be skipped */

    const char read_fmt[] = "%lg"; /* read a double */
    const char skip_fmt[] = "%*lg "; /* ignore a double */

    base = (char *) calloc(MAX_RECORD, sizeof(char));

    skip = (char *) calloc(MAX_RECORD, sizeof(char));

    skip1 = (char *) calloc(MAX_RECORD, sizeof(char));

    fmt = (char **) calloc(defs.ngenes, sizeof(char *));

    if (defs.egenes > 0)
        fmt1 = (char **) calloc(defs.egenes, sizeof(char *));

    tempparm = (double *) calloc(defs.ngenes, sizeof(double));

    tempparm1 = (double *) calloc(defs.egenes, sizeof(double));

    /* create format strings according to the number of genes */

    for (i = 0; i < defs.ngenes; i++) {
        fmt[i] = (char *) calloc(MAX_RECORD, sizeof(char));
        fmt[i] = strcpy(fmt[i], skip);
        fmt[i] = strcat(fmt[i], read_fmt);
        skip = strcat(skip, skip_fmt);
    }

    /* create format strings according to the number of external inputs */

    for (i = 0; i < defs.egenes; i++) {
        fmt1[i] = (char *) calloc(MAX_RECORD, sizeof(char));
        fmt1[i] = strcpy(fmt1[i], skip1);
        fmt1[i] = strcat(fmt1[i], read_fmt);
        skip1 = strcat(skip1, skip_fmt);
    }

    /* initialize the EqParm struct */

    l_parm.R = (double *) calloc(defs.ngenes, sizeof(double));
    l_parm.T = (double *) calloc(defs.ngenes * defs.ngenes, sizeof(double));

    if (defs.egenes > 0)
        l_parm.E = (double *) calloc(defs.ngenes * defs.egenes, sizeof(double));

    l_parm.m = (double *) calloc(defs.ngenes, sizeof(double));
    l_parm.h = (double *) calloc(defs.ngenes, sizeof(double));
    if ((defs.diff_schedule == 'A') || (defs.diff_schedule == 'C')) {
        l_parm.d = (double *) malloc(sizeof(double));
    } else {
        l_parm.d = (double *) calloc(defs.ngenes, sizeof(double));
    }
    l_parm.lambda = (double *) calloc(defs.ngenes, sizeof(double));
    l_parm.tau = (double *) calloc(defs.ngenes, sizeof(double));

    fp = FindSection(fp, section_title); /* find eqparms section */
    if (!fp)
        error("ReadParameters: cannot locate %s section", section_title);

    while (strncmp((base = fgets(base, MAX_RECORD, fp)), "$$", 2)) {

        record = base;
        lead_punct = 0; /* of string */

        c = (int) *record;
        while (c != '\0') {

            if (isdigit(c)) { /* line contains data */
                record = base;

                /* usually read ngenes parameters, but for diff. schedule A or C only read *
                 * one d parameter                                                         */

                /* If no external inputs, we must be on the next set of parameters, therefore increase linecount by one. */
                if ((linecount == 2) && (defs.egenes == 0)) {
                    linecount++;
                }

                if ((linecount == 5)
                        && ((defs.diff_schedule == 'A')
                                || (defs.diff_schedule == 'C'))) {
                    if (1 != sscanf(record, fmt[0], &tempparm[0]))
                        error("ReadParameters: error reading parms");
                } else if (linecount == 2) {
                    for (i = 0; i < defs.egenes; i++) {
                        if (1 != sscanf(record, fmt1[i], &tempparm1[i]))
                            error("ReadParameters: error reading parms");
                    }
                } else {
                    for (i = 0; i < defs.ngenes; i++) {
                        if (1 != sscanf(record, fmt[i], &tempparm[i]))
                            error("ReadParameters: error reading parms");
                    }
                }

                switch (linecount) { /* copy read parameters into the right array */
                case 0:
                    for (i = 0; i < defs.ngenes; i++) /* R */
                        l_parm.R[i] = tempparm[i];
                    linecount++;
                    break;
                case 1: /* T: keep track of read lines with Tcount */
                    for (i = 0; i < defs.ngenes; i++)
                        l_parm.T[i + Tcount * defs.ngenes] = tempparm[i];
                    Tcount++;
                    if (Tcount == defs.ngenes)
                        linecount++;
                    break;
                case 2: /* E: keep track of read lines with Ecount */
                    if (defs.egenes > 0) {
                        for (i = 0; i < defs.egenes; i++)
                            l_parm.E[i + Ecount * defs.egenes] = tempparm1[i];
                        Ecount++;
                        if (Ecount == defs.ngenes)
                            linecount++;
                    }
                    break;
                case 3: /* m */
                    for (i = 0; i < defs.ngenes; i++)
                        l_parm.m[i] = tempparm[i];
                    linecount++;
                    break;
                case 4:
                    for (i = 0; i < defs.ngenes; i++) /* h */
                        l_parm.h[i] = tempparm[i];
                    linecount++;
                    break;
                case 5: /* d: consider diff. schedule */
                    if ((defs.diff_schedule == 'A')
                            || (defs.diff_schedule == 'C')) {
                        l_parm.d[0] = tempparm[0];
                    } else {
                        for (i = 0; i < defs.ngenes; i++)
                            l_parm.d[i] = tempparm[i];
                    }
                    linecount++;
                    break;
                case 6: /* lambda */
                    for (i = 0; i < defs.ngenes; i++) {
                        l_parm.lambda[i] = tempparm[i];
                        l_parm.lambda[i] = log(2.) / l_parm.lambda[i];
                    } /* conversion done here */
                    linecount++;
                    break;
                case 7: /* tau */
                    for (i = 0; i < defs.ngenes; i++)
                        l_parm.tau[i] = tempparm[i];
                    linecount++;
                    break;
                default:
                    error("ReadParameters: too many lines in parameter section");
                    break;
                }
                break; /* don't do rest of loop anymore! */
            } else if (isalpha(c)) { /* letter means comment */
                break;
            } else if (c == '-') { /* next two elsifs for punct */
                if (((int) *(record + 1)) == '.')
                    record++;
                lead_punct = 1;
                c = (int) *(++record);
            } else if (c == '.') {
                lead_punct = 1;
                c = (int) *(++record);
            } else if (ispunct(c)) { /* other punct means comment */
                break;
            } else if (isspace(c)) { /* ignore leading white space */
                if (lead_punct) /* white space after punct means */
                    break; /* comment */
                else {
                    c = (int) *(++record); /* get next character in record */
                }
            } else {
                error("ReadParameters: illegal character in %s");
            }
        }
    }

    free(tempparm);
    free(tempparm1);
    free(base);
    free(skip);
    free(skip1);

    for (i = 0; i < defs.ngenes; i++)
        free(fmt[i]);
    free(fmt);

    for (i = 0; i < defs.egenes; i++)
        free(fmt1[i]);

    if (defs.egenes > 0)
        free(fmt1);

    return l_parm;
}

zygotic::~zygotic()
{
    // TODO free parm
    // TODO free bot
    // TODO free bot2
    // TODO free vinput
    // TODO free D
}

/*** p_deriv: This is the DvdtOrig, the original derivative function; ******
 *            implements the equations  as published in Reinitz & Sharp    *
 *            (1995), Mech Dev 49, 133-58 plus different g(u) functions    *
 *            as used by Yousong Wang in spring 2002.                      *
 ***************************************************************************/
void zygotic::p_deriv(double* v, double t, double* vdot, int n)
{
    double vinput1 = 0;
    int m; /* number of nuclei */
    int ap; /* nuclear index on AP axis [0,1,...,m-1] */
    int i, j; /* local loop counters */
    int k; /* index of gene k in a specific nucleus */
    int base, base1; /* index of first gene in a specific nucleus */
    int *l_rule; /* for autonomous implementation */
#ifdef ALPHA_DU
    int incx = 1; /* increment step size for vsqrt input array */
    int incy = 1; /* increment step size for vsqrt output array */
#endif

    // The following 3 static variables have been added dvdt_ prefix and
    // become (non-static) member of zygotic
    //static int dvdt_num_nucs = 0; /* store the number of nucs for next step */
    //static int dvdt_bcd_index = 0; /* the *next* array in bicoid struct for bcd */
    //static DArrPtr dvdt_bcd; /* pointer to appropriate bicoid struct */
    double *v_ext; /* array to hold the external input
     concentrations at time t */

    /* get D parameters and bicoid gradient according to cleavage cycle */

    m = n / defs.ngenes; /* m is the number of nuclei */
    if (m != dvdt_num_nucs) { /* time-varying quantities only vary by ccycle */
        TheMaternal.GetD(t, lparm.d, defs.diff_schedule, D);
        /* get diff coefficients, according to diff schedule */
        if (dvdt_num_nucs > m) /* started a new iteration in score */
            dvdt_bcd_index = 0; /* -> start all over again! */
        dvdt_num_nucs = m; /* store # of nucs for next step */
        dvdt_bcd = TheMaternal.GetBicoid(t, genindex); /* get bicoid gradient */
        if (dvdt_bcd.size != dvdt_num_nucs)
            error("DvdtOrig: %d nuclei don't match Bicoid!", dvdt_num_nucs);
        dvdt_bcd_index++; /* store index for next bicoid gradient */
    }

    l_rule = (int *) calloc(defs.ngenes, sizeof(int));
    for (i = 0; i < defs.ngenes; i++)
        l_rule[i] = !(TheMaternal.Theta(t));

    /* Here we retrieve the external input concentrations into v_ext */
    /* 01/13/10 Manu: If there are no external inputs, don't
     * allocate v_ext or populate it with external input
     * concentrations */

    if (defs.egenes > 0) {
        v_ext = (double *) calloc(m * defs.egenes, sizeof(double));
        ExternalInputs(t, t, v_ext, m * defs.egenes);
    }

    /* This is how it works (by JR):

     ap      nucleus position on ap axis
     base    index of first gene (0) in current nucleus
     k       index of gene k in current nucleus

     Protein synthesis terms are calculated according to g(u)

     First we do loop for vinput contributions; vinput contains
     the u that goes into g(u)

     Then we do a separate loop or vector func for sqrt or exp

     Then we do vdot contributions (R and lambda part)

     Then we do the spatial part (Ds); we do a special case for
     each end

     Note, however, that before the real loop, we have to do a
     special case for when rhere is one nuc, hence no diffusion

     These loops look a little funky 'cause we don't want any
     divides'                                                       */

    /* 01/13/10 Manu: If there are no external inputs, the for
     * loops for updating vinput1 with external input terms
     * are never entered */

    /***************************************************************************
     *                                                                         *
     *        g(u) = 1/2 * ( u / sqrt(1 + u^2) + 1)                            *
     *                                                                         *
     ***************************************************************************/
    if (gofu == Sqrt) {

        for (base = 0, base1 = 0, ap = 0; base < n;
                base += defs.ngenes, base1 += defs.egenes, ap++) {
            for (i = base; i < base + defs.ngenes; i++) {

                k = i - base;

                vinput1 = lparm.h[k];
                vinput1 += lparm.m[k] * dvdt_bcd.array[ap]; /* ap is nuclear index */

                for (j = 0; j < defs.egenes; j++)
                    vinput1 += lparm.E[(k * defs.egenes) + j]
                            * v_ext[base1 + j];

                for (j = 0; j < defs.ngenes; j++)
                    vinput1 += lparm.T[(k * defs.ngenes) + j] * v[base + j];

                bot2[i] = 1 + vinput1 * vinput1;
                vinput[i] = vinput1;
            }
        }

        /* now calculate sqrt(1+u2); store it in bot[] */
#ifdef ALPHA_DU
        vsqrt_(bot2,&incx,bot,&incy,&n); /* superfast DEC vector function */
#else
        for (i = 0; i < n; i++) /* slow traditional style sqrt */
            bot[i] = sqrt(bot2[i]);
#endif
        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if (n == defs.ngenes) { /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for (base = 0; base < n; base += defs.ngenes) {
                for (i = base; i < base + defs.ngenes; i++) {

                    k = i - base;
                    vdot1 = -lparm.lambda[k] * v[i];
                    g1 = 1 + vinput[i] / bot[i];
                    vdot1 += l_rule[k] * lparm.R[k] * 0.5 * g1;
                    vdot[i] = vdot1;
                }
            }

        } else { /* then for multiple nuclei -> diffusion */

            register double vdot1, g1;

            for (i = 0; i < defs.ngenes; i++) { /* first anterior-most nucleus */
                k = i;
                vdot1 = -lparm.lambda[k] * v[i];
                g1 = 1 + vinput[i] / bot[i];
                vdot1 += l_rule[k] * lparm.R[k] * 0.5 * g1;
                vdot1 += D[i] * (v[i + defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for (base = defs.ngenes; base < n - defs.ngenes; base +=
                    defs.ngenes) {
                for (i = base; i < base + defs.ngenes; i++) {
                    k = i - base;
                    vdot1 = -lparm.lambda[k] * v[i];
                    g1 = 1 + vinput[i] / bot[i];
                    vdot1 += l_rule[k] * lparm.R[k] * 0.5 * g1;
                    vdot1 += D[k]
                            * ((v[i - defs.ngenes] - v[i])
                                    + (v[i + defs.ngenes] - v[i]));
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for (i = base; i < base + defs.ngenes; i++) {
                k = i - base;
                vdot1 = -lparm.lambda[k] * v[i];
                g1 = 1 + vinput[i] / bot[i];
                vdot1 += l_rule[k] * lparm.R[k] * 0.5 * g1;
                vdot1 += D[k] * (v[i - defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
        }

        /***************************************************************************
         *                                                                         *
         *        g(u) = 1/2 * (tanh(u) + 1) )                                     *
         *                                                                         *
         ***************************************************************************/

    } else if (gofu == Tanh) {

        for (base = 0, base1 = 0, ap = 0; base < n;
                base += defs.ngenes, base1 += defs.egenes, ap++) {
            for (i = base; i < base + defs.ngenes; i++) {

                k = i - base;

                vinput1 = lparm.h[k];
                vinput1 += lparm.m[k] * dvdt_bcd.array[ap]; /* ap is nuclear index */

                for (j = 0; j < defs.egenes; j++)
                    vinput1 += lparm.E[(k * defs.egenes) + j]
                            * v_ext[base1 + j];

                for (j = 0; j < defs.ngenes; j++)
                    vinput1 += lparm.T[(k * defs.ngenes) + j] * v[base + j];

                vinput[i] = vinput1;
            }
        }

        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if (n == defs.ngenes) { /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for (base = 0; base < n; base += defs.ngenes) {
                for (i = base; i < base + defs.ngenes; i++) {

                    k = i - base;
                    vdot1 = -lparm.lambda[k] * v[i];
                    g1 = tanh(vinput[i]) + 1;
                    vdot1 += l_rule[k] * lparm.R[k] * 0.5 * g1;
                    vdot[i] = vdot1;
                }
            }

        } else { /* then for multiple nuclei -> diffusion */

            register double vdot1, g1;

            for (i = 0; i < defs.ngenes; i++) { /* first anterior-most nucleus */
                k = i;
                vdot1 = -lparm.lambda[k] * v[i];
                g1 = tanh(vinput[i]) + 1;
                vdot1 += l_rule[k] * lparm.R[k] * 0.5 * g1;
                vdot1 += D[i] * (v[i + defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for (base = defs.ngenes; base < n - defs.ngenes; base +=
                    defs.ngenes) {
                for (i = base; i < base + defs.ngenes; i++) {
                    k = i - base;
                    vdot1 = -lparm.lambda[k] * v[i];
                    g1 = tanh(vinput[i]) + 1;
                    vdot1 += l_rule[k] * lparm.R[k] * 0.5 * g1;
                    vdot1 += D[k]
                            * ((v[i - defs.ngenes] - v[i])
                                    + (v[i + defs.ngenes] - v[i]));
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for (i = base; i < base + defs.ngenes; i++) {
                k = i - base;
                vdot1 = -lparm.lambda[k] * v[i];
                g1 = tanh(vinput[i]) + 1;
                vdot1 += l_rule[k] * lparm.R[k] * 0.5 * g1;
                vdot1 += D[k] * (v[i - defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
        }

        /***************************************************************************
         *                                                                         *
         *        g(u) = 1 / (1 + exp(-2u))                                        *
         *                                                                         *
         ***************************************************************************/

    } else if (gofu == Exp) {

        for (base = 0, base1 = 0, ap = 0; base < n;
                base += defs.ngenes, base1 += defs.egenes, ap++) {
            for (i = base; i < base + defs.ngenes; i++) {

                k = i - base;

                vinput1 = lparm.h[k];
                vinput1 += lparm.m[k] * dvdt_bcd.array[ap]; /* ap is nuclear index */

                for (j = 0; j < defs.egenes; j++)
                    vinput1 += lparm.E[(k * defs.egenes) + j]
                            * v_ext[base1 + j];

                for (j = 0; j < defs.ngenes; j++)
                    vinput1 += lparm.T[(k * defs.ngenes) + j] * v[base + j];

                vinput[i] = -2.0 * vinput1;
            }
        }

        /* now calculate exp(-u); store it in bot[] */
#ifdef ALPHA_DU
        vexp_(vinput,&incx,bot,&incy,&n); /* superfast DEC vector function */
#else
        for (i = 0; i < n; i++) /* slow traditional style exp */
            bot[i] = exp(vinput[i]);
#endif

        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if (n == defs.ngenes) { /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for (base = 0; base < n; base += defs.ngenes) {
                for (i = base; i < base + defs.ngenes; i++) {

                    k = i - base;
                    vdot1 = -lparm.lambda[k] * v[i];
                    g1 = 1 / (1 + bot[i]);
                    vdot1 += l_rule[k] * lparm.R[k] * g1;
                    vdot[i] = vdot1;
                }
            }

        } else { /* then for multiple nuclei -> diffusion */

            register double vdot1, g1;

            for (i = 0; i < defs.ngenes; i++) { /* first anterior-most nucleus */
                k = i;
                vdot1 = -lparm.lambda[k] * v[i];
                g1 = 1 / (1 + bot[i]);
                vdot1 += l_rule[k] * lparm.R[k] * g1;
                vdot1 += D[i] * (v[i + defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for (base = defs.ngenes; base < n - defs.ngenes; base +=
                    defs.ngenes) {
                for (i = base; i < base + defs.ngenes; i++) {
                    k = i - base;
                    vdot1 = -lparm.lambda[k] * v[i];
                    g1 = 1 / (1 + bot[i]);
                    vdot1 += l_rule[k] * lparm.R[k] * g1;
                    vdot1 += D[k]
                            * ((v[i - defs.ngenes] - v[i])
                                    + (v[i + defs.ngenes] - v[i]));
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for (i = base; i < base + defs.ngenes; i++) {
                k = i - base;
                vdot1 = -lparm.lambda[k] * v[i];
                g1 = 1 / (1 + bot[i]);
                vdot1 += l_rule[k] * lparm.R[k] * g1;
                vdot1 += D[k] * (v[i - defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
        }

        /***************************************************************************
         *                                                                         *
         *        g(u) = 0 if u<0, 1 if u>=0 (Heaviside function)                  *
         *                                                                         *
         *        this makes the model quasi boolean and the equations locally     *
         *        linear for both u>0 and u<0                                      *
         ***************************************************************************/

    } else if (gofu == Hvs) {

        for (base = 0, base1 = 0, ap = 0; base < n;
                base += defs.ngenes, base1 += defs.egenes, ap++) {
            for (i = base; i < base + defs.ngenes; i++) {

                k = i - base;

                vinput1 = lparm.h[k];
                vinput1 += lparm.m[k] * dvdt_bcd.array[ap]; /* ap is nuclear index */

                for (j = 0; j < defs.egenes; j++)
                    vinput1 += lparm.E[(k * defs.egenes) + j]
                            * v_ext[base1 + j];

                for (j = 0; j < defs.ngenes; j++)
                    vinput1 += lparm.T[(k * defs.ngenes) + j] * v[base + j];
                vinput[i] = vinput1;
            }
        }

        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if (n == defs.ngenes) { /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for (base = 0; base < n; base += defs.ngenes) {
                for (i = base; i < base + defs.ngenes; i++) {

                    k = i - base;
                    vdot1 = -lparm.lambda[k] * v[i];
                    if (vinput[i] >= 0.)
                        g1 = 1.;
                    else
                        g1 = 0.;
                    vdot1 += l_rule[k] * lparm.R[k] * g1;
                    vdot[i] = vdot1;
                }
            }

        } else { /* then for multiple nuclei -> diffusion */

            register double vdot1, g1;

            for (i = 0; i < defs.ngenes; i++) { /* first anterior-most nucleus */
                k = i;
                vdot1 = -lparm.lambda[k] * v[i];
                if (vinput[i] >= 0.)
                    g1 = 1.;
                else
                    g1 = 0.;
                vdot1 += l_rule[k] * lparm.R[k] * g1;
                vdot1 += D[i] * (v[i + defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for (base = defs.ngenes; base < n - defs.ngenes; base +=
                    defs.ngenes) {
                for (i = base; i < base + defs.ngenes; i++) {
                    k = i - base;
                    vdot1 = -lparm.lambda[k] * v[i];
                    if (vinput[i] >= 0.)
                        g1 = 1.;
                    else
                        g1 = 0.;
                    vdot1 += l_rule[k] * lparm.R[k] * g1;
                    vdot1 += D[k]
                            * ((v[i - defs.ngenes] - v[i])
                                    + (v[i + defs.ngenes] - v[i]));
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for (i = base; i < base + defs.ngenes; i++) {
                k = i - base;
                vdot1 = -lparm.lambda[k] * v[i];
                if (vinput[i] >= 0.)
                    g1 = 1.;
                else
                    g1 = 0.;
                vdot1 += l_rule[k] * lparm.R[k] * g1;
                vdot1 += D[k] * (v[i - defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
        }

    } else
        error("DvdtOrig: unknown g(u)");

    /* during mitosis only diffusion and decay happen */

    free(l_rule);

    if (defs.egenes > 0)
        free(v_ext);

    return;
}
