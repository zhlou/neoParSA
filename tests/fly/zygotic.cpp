/*
 * zygotic.cpp
 *
 *  Created on: Feb 5, 2013
 *      Author: zhlou
 */

#include "zygotic.h"
#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <typeinfo>
#include <cassert>
#include <cstring>
#include <cctype>
#include <cmath>
#include <cassert>
using namespace std;


const int zygotic::ADD_BIAS = 1;
const int zygotic::NO_OP = 2;
const int zygotic::DIVIDE = 4;
const int zygotic::PROPAGATE = 8;
const int zygotic::MITOTATE = 16;

zygotic::zygotic(maternal& in_maternal, FILE *fp, const char* parm_section,
                 GFunc in_gofu, int in_debug, const char *solver_name) :
        TheMaternal(in_maternal), defs(in_maternal.getProblem()),
        gofu(in_gofu), debug(in_debug)
{
    dvdt_num_nucs = 0;
    dvdt_bcd_index = 0;
    parm = ReadParameters(fp, parm_section);
    D = (double *) calloc(defs.ngenes, sizeof(double));
    vinput = (double *) calloc(defs.ngenes * defs.nnucs, sizeof(double));
    bot2 = (double *) calloc(defs.ngenes * defs.nnucs, sizeof(double));
    bot = (double *) calloc(defs.ngenes * defs.nnucs, sizeof(double));
    solve = SolverFactory(*this, debug, solver_name);

}

/*** A FUNCTION THAT READS PARAMETERS FROM THE DATA FILE INTO STRUCTS ******/

/*** ReadParamters: reads the parameters for a simulation run from the *****
 *                  eqparms or input section of the data file as indicated *
 *                  by the section_title argument and does the conversion  *
 *                  of protein half lives into lambda parameters.          *
 ***************************************************************************/
EqParms zygotic::ReadParameters(FILE* fp, const char* section_title)
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
    free(D);
    free(vinput);
    free(bot2);
    free(bot);
    free(parm.R);
    free(parm.T);
    if (defs.egenes > 0)
        free(parm.E);
    free(parm.m);
    free(parm.h);
    free(parm.d);
    free(parm.lambda);
    free(parm.tau);
    // TODO free parm
    delete solve;
}

/*** THE FOLLOWING FUNCTIONS ARE FOR HANDLING TLIST, a linked list used to *
 *   initialize the structure that tells the solver for which time there's *
 *   data, how many nuclei there are and what to do.                       *
 ***************************************************************************/

/*** InitTList: initializes TList and adds first (t=0) and last ************
 *              (t=gastrulation) element of Tlist.                         *
 ***************************************************************************/
zygotic::TList* zygotic::InitTList(void)
{
    TList *first, *last;                  /* first and last element of TList */

  /* initialize the first element (t=0) */

    first = (TList *)malloc(sizeof(TList));
    first->time = 0.+EPSILON;
    first->n = defs.ngenes * TheMaternal.GetNNucs(0);
    first->op = PROPAGATE;

  /* initialize the last element (t=gast_time) */

    last = (TList *)malloc(sizeof(TList));
    if ( !(last->time = TheMaternal.GetGastTime()) )
      error("InitTList: error getting gastrulation time");
    last->n = defs.ngenes * TheMaternal.GetNNucs(last->time);
    last->op = NO_OP;

  /* link the two elements and return the list */

    first->next = last;
    last->next = NULL;

    return first;
}

/*** InsertTList: takes pointer to first element of TList plus time and ****
 *                desired op for new TList element and inserts a new TList *
 *                element for time at the appropriate place within the     *
 *                linked list. This function returns a pointer to the      *
 *                first element of the TList if insertion worked out fine. *
 ***************************************************************************/

zygotic::TList* zygotic::InsertTList(TList* first, double time, int op)
{
    TList *current;                            /* used to step through TList */
    TList *newlist;                   /* used to allocate memory for new element */

    int n;                                      /* how many nuclei at 'time' */



    n = defs.ngenes * TheMaternal.GetNNucs(time);     /* get number of nuclei for 'time' */

    if (first == NULL)
        error("InsertTList: TList pointer is NULL at time %f!", time);

    /* the following loop steps through the linked list and places the new     *
     * element at the appropriate position depending on its time               */

    current = first;
    do {
        if( fabs(time - current->time) < BIG_EPSILON ) {
            if ( current->n  ==  n ) {  /* if new time really close to an exist- */
                if ( current->op == MITOTATE) {
                    /* if we are skipping epsilon
                                      after MITOSIS start */
                    newlist = (TList *)malloc(sizeof(TList));
                    newlist->time = time;  /* -> allocate new element for right after */
                    newlist->n = n;               /* cell division has occurred */
                    newlist->op = op;
                    newlist->next = current->next;
                    current->next = newlist;
                    return first;
                }

                current->op |= op; /* ing time point and no cell division happened */
                return first;                /* -> just add new op to existing one */
            }
            else if ( current->n < n ) {                /* if, on the other hand */
                newlist = (TList *)malloc(sizeof(TList));     /* cell div HAS happened */
                newlist->time = time;       /* -> allocate new element for right after */
                newlist->n = n;                          /* cell division has occurred */
                newlist->op = op;
                newlist->next = current->next;
                current->next = newlist;
                return first;
            }
            else         /* a sudden reduction of nuclei will be hard to explain */
                error("InsertTList: sudden reduction of nuclei at time %g!",
                        current->time);
        }
        else if ( time > current->time ) {     /* new time > than current time */
            if ( fabs(time - current->next->time) < BIG_EPSILON ) {     /* is it */
                current = current->next;   /* really close to the next time point? */
            }                                                    /* -> go there! */
            else if ( time > current->next->time ) {         /* or if time is >> */
                current = current->next;               /* than the next time point */
            }                                                /* -> go there too! */
            else if ( time < current->next->time ) {         /* but if time is < */
                newlist = (TList *)malloc(sizeof(TList));       /* the next time point */
                newlist->time = time;                       /* -> allocate new element */
                newlist->n = n;
                newlist->op = op;
                newlist->next = current->next;
                current->next = newlist;
                return first;
            }
            else {      /* if all the above don't apply, there's something wrong */
                error("InsertTList: impossible error at time %g!", current->time);
            }
        }
        else {          /* if time < current time, there's something wrong too */
            error("InsterTList: we missed our exit at time %g!", current->time);
        }
    } while (current != NULL);

    return current;
}

/*** CountEntries: counts how many entries we have in a TList **************
 ***************************************************************************/
int zygotic::CountEntries(TList* first)
{
    int n = 0;
    while(first != NULL) {
      n++;
      first = first->next;
    }
    return n;
}

void zygotic::FreeTList(TList* first)
{
    if(first->next != NULL)
      FreeTList(first->next);

    free(first);
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

    SoDe *delay_solver = NULL;
    if (typeid(*solve) == typeid(SoDe)){
        delay_solver = dynamic_cast<SoDe *>(solve);
        if (defs.egenes > 0) {
            v_ext = (double *) calloc(m * defs.egenes, sizeof(double));
            delay_solver->ExternalInputs(t, t, v_ext, m * defs.egenes);
        }
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

/*** JACOBIAN FUNCTION(S) **************************************************/

/*** JacobnOrig: Jacobian function for the DvdtOrig model; calculates the **
 *               Jacobian matrix (matrix of partial derivatives) for the   *
 *               equations at a give time t; input concentrations come in  *
 *               v (of size n), the Jacobian is returned in jac; note that *
 *               all dfdt's are zero in our case since our equations are   *
 *               autonomous (i.e. have no explicit t in them)              *
 ***************************************************************************
 *                                                                         *
 * The Equation: dv/dx = Rg(u) + D() + lambda()                            *
 *                                                                         *
 * The Jacobian: (in this example: nnucs = 3, ngenes = 3                   *
 *                                                                         *
 *                        i = 1       i = 2        i = 3                   *
 *         b =          1   2   3   1   2   3    1   2   3                 *
 *                                                                         *
 *         a = 1      X11 Y12 Y13  D1   0   0    0   0   0                 *
 * i = 1   a = 2      Y21 X22 Y23   0  D2   0    0   0   0                 *
 *         a = 3      Y31 Y32 X33   0   0  D3    0   0   0                 *
 *                                                                         *
 *         a = 1       D1   0   0 X11 Y12 Y13   D1   0   0                 *
 * i = 2   a = 2        0  D2   0 Y21 X22 Y23    0  D2   0                 *
 *         a = 3        0   0  D3 Y31 Y32 X33    0   0  D3                 *
 *                                                                         *
 *         a = 1        0   0   0  D1   0   0  X11 Y12 Y13                 *
 * i = 3   a = 2        0   0   0   0  D2   0  Y21 X22 Y23                 *
 *         a = 3        0   0   0   0   0   0  Y31 Y32 Y33                 *
 *                                                                         *
 * Where:  Xab = Rg'(u) - 2D{a} - lambda{a}                                *
 *         Yab = Rg'(u)                                                    *
 *         Da  = D{a}                                                      *
 *                                                                         *
 ***************************************************************************/

void zygotic::p_jacobn(double t, double *v, double *dfdt, double **jac, int n)
{
  int        m;                                        /* number of nuclei */
  int        ap;                 /* nuclear index on AP axis [0,1,...,m-1] */
  int        i, j;                                  /* local loop counters */
  int        k, kk;               /* index of gene k in a specific nucleus */
  int        base;            /* indices of 1st gene in a specific nucleus */
#ifdef ALPHA_DU
  int        incx = 1;        /* increment step size for vsqrt input array */
  int        incy = 1;       /* increment step size for vsqrt output array */
#endif

  static int num_nucs  = 0;      /* store the number of nucs for next step */
  static int bcd_index = 0;   /* the *next* array in bicoid struct for bcd */
  static DArrPtr bcd;              /* pointer to appropriate bicoid struct */


/* get D parameters and bicoid gradient according to cleavage cycle */

  m = n / defs.ngenes;                        /* m is the number of nuclei */
  if (m != num_nucs) {      /* time-varying quantities only vary by ccycle */
    TheMaternal.GetD(t, lparm.d, defs.diff_schedule, D);
                      /* get diff coefficients, according to diff schedule */
    if (num_nucs > m)                  /* started a new iteration in score */
      bcd_index = 0;                           /* -> start all over again! */
    num_nucs = m;                         /* store # of nucs for next step */
    bcd = TheMaternal.GetBicoid(t, genindex);                   /* get bicoid gradient */
    if( bcd.size != num_nucs)
     error("JacobnOrig: %d nuclei don't match Bicoid!", num_nucs);
    bcd_index++;                   /* store index for next bicoid gradient */
  }

/*** INTERPHASE rule *******************************************************/

  if (rule == INTERPHASE) {

    register double vinput1 = 0;                    /* used to calculate u */
    register double gdot1, vdot1;     /* used to calculate Xab's and Yab's */


/***************************************************************************
 *                                                                         *
 *  g(u)  = 1/2 * ( u / sqrt(1 + u^2) + 1)                                 *
 *  g'(u) = 1/2 * ( 1 / ((1 + u^2)^3/2)) * T{ab}                           *
 *                                                                         *
 ***************************************************************************/

    if ( gofu == Sqrt ) {

/* this looks confusing, but it's actually quite simple: the two loops be- *
 * low are in reality just one that loops over the first dimension of the  *
 * Jacobian (size n) 'nuclear block' by 'nuclear block' (the i's in the    *
 * long introductory comment above); ap keeps track of the nucleus number, *
 * k keeps track of which gene we're dealing with; bot2 saves 1+u^2        */

      for (base=0, ap=0; base<n ; base+=defs.ngenes, ap++) {
    for (i=base; i < base+defs.ngenes; i++) {

      k = i - base;

      vinput1  = lparm.h[k];
      vinput1 += lparm.m[k] * bcd.array[ap];

      for(j=0; j < defs.ngenes; j++)
        vinput1 += lparm.T[(k*defs.ngenes)+j] * v[base + j];

      bot2[i] = 1 + vinput1 * vinput1;
    }
      }

/* now calculate sqrt(1+u^2); store it in bot[] */

#ifdef ALPHA_DU
      vsqrt_(bot2,&incx,bot,&incy,&n);    /* superfast DEC vector function */
#else
      for(i=0; i < n; i++)                  /* slow traditional style sqrt */
    bot[i] = sqrt(bot2[i]);
#endif

/* resume loop after vector sqrt above; we finish calculating g'(u) and    *
 * place D's in the diagonal of the off-diagonal blocks (cf diagram above) */

      for (base=0; base < n; base+=defs.ngenes) {
          for (i=base; i < base+defs.ngenes; i++) {

              k = i - base;

              gdot1  = 1 / (bot[i] * bot2[i]);
              gdot1 *= lparm.R[k] * 0.5;

              for (j=base; j < base+defs.ngenes; j++) {

                  kk = j - base;

                  vdot1 = lparm.T[(k*defs.ngenes)+kk] * gdot1;
                  if ( k == kk ) {
                      if ( n > defs.ngenes ) {
                          if ( base > 0 && base < n-defs.ngenes )
                              vdot1 -= 2. * D[k];
                          else
                              vdot1 -= D[k];
                      }
                      vdot1 -= lparm.lambda[k];
                  }
                  jac[i][j] = vdot1;

              }

              if ( base > 0 )
                  jac[i][i-defs.ngenes] = D[k];

              if ( base < n-defs.ngenes )
                  jac[i][i+defs.ngenes] = D[k];

          }
      }

/***************************************************************************
 *                                                                         *
 * g(u)  = 1/2 * (tanh(u) + 1) )  or  g(u)  = 1 / ( 1 + e^(-2u))           *
 * g'(u) = Sech(u)^2 /2           or  g'(u) = 2e^(-2u) / (1 + e^(-2u))^2   *
 *                                                                         *
 ***************************************************************************
 *                                                                         *
 * These are actually the same function in different guises; we implement  *
 * Exp below since it's faster                                             *
 *                                                                         *
 ***************************************************************************/

    } else if ( (gofu == Tanh) || (gofu == Exp) ) {

/* this looks confusing, but it's actually quite simple: the two loops be- *
 * low are in reality just one that loops over the first dimension of the  *
 * Jacobian (size n) 'nuclear block' by 'nuclear block' (the i's in the    *
 * long introductory comment above); ap keeps track of the nucleus number, *
 * k keeps track of which gene we're dealing with; bot2 saves 1+u^2        */

      for (base=0, ap=0; base<n ; base+=defs.ngenes, ap++) {
    for (i=base; i < base+defs.ngenes; i++) {

      k = i - base;

      vinput1  = lparm.h[k];
      vinput1 += lparm.m[k] * bcd.array[ap];

      for(j=0; j < defs.ngenes; j++)
        vinput1 += lparm.T[(k*defs.ngenes)+j] * v[base + j];

      bot[i] = -2.0 * vinput1;
    }
      }

/* now calculate exp(-u); store it in bot[] */
#ifdef ALPHA_DU
      vexp_(bot,&incx,bot2,&incy,&n);   /* superfast DEC vector function */
#else
      for (i=0; i<n; i++)                    /* slow traditional style exp */
    bot2[i] = exp(bot[i]);
#endif

/* resume loop after vector exp above; we finish calculating g'(u) and     *
 * place D's in the diagonal of the off-diagonal blocks (cf diagram above) */

      for (base=0; base < n; base+=defs.ngenes) {
          for (i=base; i < base+defs.ngenes; i++) {

              k = i - base;

              gdot1  = 2. * bot2[i];
              gdot1 /= (1. + bot2[i]) * (1. + bot2[i]);
              gdot1 *= lparm.R[k];

              for (j=base; j < base+defs.ngenes; j++) {

                  kk = j - base;

                  vdot1 = lparm.T[(k*defs.ngenes)+kk] * gdot1;
                  if ( k == kk ) {
                      if ( n > defs.ngenes ) {
                          if ( base > 0 && base < n-defs.ngenes )
                              vdot1 -= 2. * D[k];
                          else
                              vdot1 -= D[k];
                      }
                      vdot1 -= lparm.lambda[k];
                  }
                  jac[i][j] = vdot1;

              }

              if ( base > 0 )
                  jac[i][i-defs.ngenes] = D[k];

              if ( base < n-defs.ngenes )
                  jac[i][i+defs.ngenes] = D[k];

          }
      }

/*** semi-implicit solvers are NOT allowed with heaviside g(u) *************/

    } else if ( gofu == Hvs ) {
      error("JacobnOrig: can't use semi-implicit solver on heaviside g(u)!");

    } else {
      error("JacobnOrig: unknown g(u) function!\n");
    }

/* during mitosis only diffusion and decay happen */

  } else if (rule == MITOSIS) {

    register double vdot1;

    for (base=0; base < n; base+=defs.ngenes) {
        for (i=base; i < base+defs.ngenes; i++) {
            k = i - base;

            for (j=base; j < base+defs.ngenes; j++) {
                kk = j - base;

                if ( k == kk ) {
                    vdot1  = -lparm.lambda[k];
                    if ( n < defs.ngenes ) {
                        if ( base > 0 && base < n-defs.ngenes )
                            vdot1 -= 2. * D[k];
                        else
                            vdot1 -= D[k];
                    }
                } else
                    vdot1 = 0.;

                jac[i][j] = vdot1;

            }

      if ( base > 0 )
        jac[i][i-defs.ngenes] = D[k];

      if ( base < n-defs.ngenes )
        jac[i][i+defs.ngenes] = D[k];

    }
      }

  } else
    error("JacobnOrig: Bad rule %i sent to JacobnOrig", rule);

  return;

}




/*** MUTATOR FUNCTIONS *****************************************************/

/*** Mutate: calls mutator functions according to genotype string **********
 ***************************************************************************/

void zygotic::Mutate(char *g_type)
{
  int  i;
  int  c;

  char *record;

  lparm = CopyParm(parm);   /* make local copy of parameters to be mutated */

  record = g_type;
  c=(int)*record;

  for ( i=0; c != '\0'; i++, c=(int)*(++record) ) {
    if ( c == 'W' )
      continue;
    else if ( c == 'R' )
      R_Mutate(i);
    else if ( c == 'S' )
      RT_Mutate(i);
    else if ( c == 'T' )
      T_Mutate(i);
    else
      error("Mutate: unrecognized letter in genotype string!");
  }
}



/*** T_Mutate: mutates genes by setting all their T matrix entries *********
 *             to zero. Used to simulate mutants that express a            *
 *             non-functional protein.                                     *
 ***************************************************************************/

void zygotic::T_Mutate(int gene)
{
  int      i;

  for( i=0; i<defs.ngenes; i++)
    lparm.T[(i*defs.ngenes)+gene]=0;
}



/*** R_Mutate: mutates genes by setting their promotor strength R **********
 *             to zero, so there will be no transcription at all           *
 *             anymore. Used to simulate mutants that don't pro-           *
 *             duce any protein anymore.                                   *
 ***************************************************************************/

void zygotic::R_Mutate(int gene)
{
    lparm.R[gene]=0;
}



/*** RT_Mutate: mutates gene by setting both promoter strength R ***********
 *              and T matrix entries to zero. Useful, if you want          *
 *              to be really sure that there is NO protein pro-            *
 *              duction and that maternal contribution also don't          *
 *              contribute anything to gene interactions.                  *
 ***************************************************************************/

void zygotic::RT_Mutate(int gene)
{
  int      i;

  for( i=0; i<defs.ngenes; i++)
    lparm.T[(i*defs.ngenes)+gene]=0;
  lparm.R[gene]=0;
  /* Don't need to zero param.thresh in (trans acting) mutants */
}


/*** CopyParm: copies all the parameters into the lparm struct *************
 ***************************************************************************/

EqParms zygotic::CopyParm(EqParms orig_parm)
{
  int           i,j;                                /* local loop counters */

  EqParms       l_parm;              /* copy of parm struct to be returned */

  l_parm.R = (double *)calloc(defs.ngenes, sizeof(double));
  l_parm.T = (double *)calloc(defs.ngenes * defs.ngenes, sizeof(double));

  if (defs.egenes > 0)
      l_parm.E = (double *)calloc(defs.ngenes * defs.egenes, sizeof(double));

  l_parm.m = (double *)calloc(defs.ngenes, sizeof(double));
  l_parm.h = (double *)calloc(defs.ngenes, sizeof(double));
  if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') ) {
    l_parm.d = (double *)malloc(sizeof(double));
  } else {
    l_parm.d = (double *)calloc(defs.ngenes, sizeof(double));
  }
  l_parm.lambda = (double *)calloc(defs.ngenes, sizeof(double));
  l_parm.tau = (double *)calloc(defs.ngenes, sizeof(double));

  for (i=0; i<defs.ngenes; i++ ) {
    l_parm.R[i] = orig_parm.R[i];
    for (j=0; j<defs.ngenes; j++ )
      l_parm.T[(i*defs.ngenes)+j] = orig_parm.T[(i*defs.ngenes)+j];
    for (j=0; j<defs.egenes; j++ )
      l_parm.E[(i*defs.egenes)+j] = orig_parm.E[(i*defs.egenes)+j];
    l_parm.m[i] = orig_parm.m[i];
    l_parm.h[i] = orig_parm.h[i];
    l_parm.lambda[i] = orig_parm.lambda[i];
    l_parm.tau[i] = orig_parm.tau[i];
  }

  if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') ) {
    l_parm.d[0] = orig_parm.d[0];
  } else {
    for (i=0; i<defs.ngenes; i++)
      l_parm.d[i] = orig_parm.d[i];
  }

  return l_parm;
}

/*** FreeMutant: frees mutated parameter struct ****************************
 ***************************************************************************/
void zygotic::FreeMutant(void)
{
  free(lparm.R);
  free(lparm.T);
  if (defs.egenes > 0)
      free(lparm.E);
  free(lparm.m);
  free(lparm.h);
  free(lparm.d);
  free(lparm.lambda);
  free(lparm.tau);
}

/*** Blastoderm: runs embryo model and returns an array of concentration ***
 *               arrays for each requested time (given by TabTimes) using  *
 *               stephint as a suggested stepsize and accuracy as the re-  *
 *               quired numerical accuracy (in case of adaptive stepsize   *
 *               solvers or global stepsize control) for the solver.       *
 *         NOTE: TabTimes *must* start from 0 and have increasing times.   *
 *               It includes times when bias is added, cell division times *
 *               and times for which we have data and ends with the time   *
 *               of gastrulation.                                          *
 ***************************************************************************/
NArrPtr zygotic::Blastoderm(int in_genindex, char *genotype,
                        InterpObject hist_interrp,
                        InterpObject extinp_interrp, DArrPtr tabtimes,
                        double stephint, double accuracy, FILE *slog)
{
  const  double    epsilon     = EPSILON;      /* epsilons: very small in- */
  const  double    big_epsilon = BIG_EPSILON; /* creases used for division */

  static NArrPtr   solution;               /* solution will be an array of */
                                          /* concs for each requested time */
  static int       *what2do;               /* what to do at each time step */

  static double    *divtable;           /* cell div times in reverse order */
  static double    *transitions;           /* this is when cell divs start */
  static double    *durations;              /* durations of cell divisions */


  static DArrPtr   biastimes;               /* times a which bias is added */
  static DArrPtr   oldtimes;                  /* statically saves tabtimes */

  static int       allocate;             /* flag: need to allocate or not? */

  int              i,ii,j;                                /* loop counters */
  int              k;                    /* index of gene k in current nuc */
  int              ap;                         /* nuc. position on AP axis */
  int              lin;             /* first lineage number at each ccycle */

  int              rule;                         /* MITOSIS or INTERPHASE? */

  DArrPtr          bias;                 /* bias for given time & genotype */

  TList            *entries    = NULL;   /* temp linked list for times and */
  TList            *current;                         /* ops for the solver */
  SoDe *delay_solver = get_delay_solver();


/* allocate transitions array (skip if ndivs == 0) */

  if ( (transitions == NULL) && (defs.ndivs > 0) )
    transitions = (double *)calloc(defs.ndivs, sizeof(double));

/* for each genotype, the 'genotype' variable has to be made static to zy- *
 * gotic.c so that the derivative functions know which genotype they're    *
 * dealing with (i.e. they need to get the appropriate bcd gradient)       */

  genindex = in_genindex;



  if (delay_solver) {
      if (defs.egenes > 0)
            delay_solver->SetExternalInputInterp(extinp_interrp);
      delay_solver->SetHistoryInterp(hist_interrp);
  }

/* Blastoderm() first checks if tabtimes has the same size as in the pre-  *
 * vious run (remembered by the static oldtimes array) and if all tab ti-  *
 * mes are the same; if this is the case, the initialization is skipped    */

  if ( tabtimes.size == oldtimes.size ) {
    for ( i=0; i<tabtimes.size; i++ )
      if ( tabtimes.array[i] != oldtimes.array[i] )
    allocate = 1;
  } else
    allocate = 1;

/* INITIALIZATION OF THE MODEL STRUCTS AND ARRAYS **************************/

/* free structures that are already present */

  if ( allocate ) {
    if ( solution.array != NULL ) {
      FreeSolution(&solution);
      solution.array = NULL;
      solution.size  = 0;
    }
    if (what2do != NULL) {
      free(what2do);
      what2do = NULL;
    }
    if (biastimes.array != NULL) {
      free(biastimes.array);
      biastimes.array = NULL;
      biastimes.size  = 0;
    }
    if (oldtimes.array != NULL) {
      free(oldtimes.array);
      oldtimes.array = NULL;
      oldtimes.size  = 0;
    }

/* get bias times and initialize information about cell divisions */

    biastimes = TheMaternal.GetBTimes(genotype);
    if ( !(biastimes.array) )
        error("Blastoderm: error getting bias times");

    if ( defs.ndivs > 0 ) {
        if ( !(divtable = TheMaternal.GetDivtable()) )
            error("Blastoderm: error getting division table");
        if ( !(durations = TheMaternal.GetDurations()) )
            error("Blastoderm: error getting division durations");
        for(i=0; i<defs.ndivs; i++)
            transitions[i] = divtable[i] - durations[i];
    }

/* entries is a linked list, which we use to set up the solution struct    *
 * and the what2do array; each of these needs an entry for:                *
 * - start and end (gastrulation) time                                     *
 * - times for mitoses: - beginning of mitosis                             *
 *                      - cell division time (still belongs to previous    *
 *                        cleavage cycle)                                  *
 *                      - time right after cell division (+EPSILON), be-   *
 *                        longs to new cell cycle with doubled nnucs       *
 * - times at which we add bias                                            *
 * - tabulated times for which we have data or which we want to display    */

/* add start and end (gastrulation) time */

    entries = InitTList();

/* add all times required for mitoses (skip this for 0 div schedule) */

    if ( defs.ndivs > 0 )
      for (i=0; i<defs.ndivs; i++) {
    entries = InsertTList(entries, divtable[i], DIVIDE);
    if( (TheMaternal.GetNNucs(divtable[i])) == (TheMaternal.GetNNucs(divtable[i] + epsilon)) )
      error("Blastoderm: epsilon of %g too small!", epsilon);
    entries = InsertTList(entries, divtable[i] + epsilon, PROPAGATE);
    entries = InsertTList(entries, transitions[i], MITOTATE);
    if( (TheMaternal.GetNNucs(transitions[i])) !=
                            (TheMaternal.GetNNucs(transitions[i] + epsilon)) )
      error("Blastoderm: division within epsilon of %g!", epsilon);
    entries = InsertTList(entries, transitions[i]+epsilon, PROPAGATE);


      }

/* add bias times */

    for(i=0; i < biastimes.size; i++)
      entries =
    InsertTList(entries, biastimes.array[i], ADD_BIAS | PROPAGATE);

/* tabulated times */

    for(i=0; i < tabtimes.size; i++)
      entries = InsertTList(entries, tabtimes.array[i], PROPAGATE);

/* now we know the number of solutions we have to calculate, so we can     *
 * allocate and initialize the solution struct and the what2do array       */

    solution.size  = CountEntries(entries);
    solution.array = (NucState *)calloc(solution.size, sizeof(NucState));

    what2do = (int *)calloc(solution.size, sizeof(int));

    current = entries;
    for(i=0; i<solution.size; i++) {
      solution.array[i].time = current->time;
      solution.array[i].state.size = current->n;
      solution.array[i].state.array =
    (double *)calloc(current->n, sizeof(double));
      for(j=0; j<current->n; j++)
    solution.array[i].state.array[j] = 0;
      what2do[i] = current->op;
      current = current->next;
    }

/* free the linked list and save new tab times in oldtimes array */

    FreeTList(entries);
    entries = NULL;

    oldtimes.size  = tabtimes.size;
    oldtimes.array = (double *)calloc(tabtimes.size, sizeof(double));
    for(i=0; i< tabtimes.size; i++)
      oldtimes.array[i] = tabtimes.array[i];

    allocate = 0;

/* if we're not allocating: reset initial conditions; we'll *add* them     *
 * again below                                                             */

  } else

    for(i=0; i < solution.array[0].state.size; i++)
      solution.array[0].state.array[i] = 0;


/* RUNNING THE MODEL *******************************************************/

/* Before running the model, mutate zygotic params appropriately */

  Mutate(genotype);

  if ( debug ) {
    fprintf(slog, "\n-----------------------------------------------");
    fprintf(slog, "----------------------------------------------\n");
    fprintf(slog, "Blastoderm: mutated genotype to %s.\n", genotype);
  }

/* Below is the loop that evaluates the solution and saves it in the solu- *
 * tion struct                                                             */

  for(i=0; i<solution.size; i++) {

/* First we have to set the propagation rule in zygotic.c (either MITOSIS  *
 * or INTERPHASE) to make sure that the correct differential equations are *
 * propagated (rules are defined in zygotic.h)                             */

    rule = TheMaternal.GetRule(solution.array[i].time);
    SetRule(rule);

    if ( debug ) {
        if ( rule == 0 )
            fprintf(slog, "Blastoderm: rule is INTERPHASE.\n");
        else if ( rule == 1 && !(what2do[i] & DIVIDE) )
            fprintf(slog, "Blastoderm: rule is MITOSIS.\n");
        else if ( what2do[i] & DIVIDE )
            fprintf(slog, "Blastoderm: rule is DIVISION.\n");
    }

/* ADD_BIAS is a special op in that it can be combined with any other op   *
 * (see also next comment); we can add (or subtract) protein conentrations *
 * at any time; this is for adding initial conditions for genes that are   *
 * both maternal and zygotic (e.g. cad and hb) or to simulated perturba-   *
 * tions like heat shocks or induced overexpression of certain genes;      *
 * note that the bias lives in maternal.c and has to be fetched from there */

    if (what2do[i] & ADD_BIAS) {
      for (j=0; j<biastimes.size; j++) {
    if ( fabs(solution.array[i].time - biastimes.array[j])
         < big_epsilon ) {
      bias = TheMaternal.GetBias(biastimes.array[j], in_genindex);
      for (ii=0; ii < bias.size; ii++)
        solution.array[i].state.array[ii] += bias.array[ii];
    }
      }

      if ( debug )
    fprintf(slog, "Blastoderm: added bias at time %f.\n",
        solution.array[i].time);
    }

/* The ops below can be executed in addition to ADD_BIAS but they cannot   *
 * be combined between themselves; if more than one op is set, the prio-   *
 * rities are as follows: NO_OP > DIVIDE > PROPAGATE; please make sure     *
 * that you do not set double ops, which will only obfuscate the code;     *
 *                                                                         *
 * NO_OP simply does nothing (used for gastrulation time)                  */

    if (what2do[i] & NO_OP)
      ;

/* This is the DIVIDE rule: only symmetrical division is supported yet,    *
 * i.e. both daughter nuclei inherit all concentrations from their mother; *
 * we have to take special care of the fact that we might not need the     *
 * most anterior and most posterior daughter nuclei, depending on the ran- *
 * ge of nuclei in cycle 14 we've chosen to anneal on; e.g. if the most    *
 * anterior nucleus in cycle 14 has an odd lineage number, we need to push *
 * its sibling off the anterior limit when we divide after cycle 13; or in *
 * other words: our most anterior cell in cycle 14 is the posterior daugh- *
 * ter of our most anterior cell in cycle 13; therefore, we need to 'loose'*
 * its anterior sibling; the same applies for the posterior end, where we  *
 * want to loose the most posterior daughter cell if we don't need it any  *
 * more at the later cycle                                                 *
 *                                                                         *
 * This is implemented as follows below: lin is the lineage number of the  *
 * most anterior cell of the next cell cycle; if it's odd numbered, we'll  *
 * push off the most anterior daughter by subtracting the number of genes  *
 * from the daughter indices (i.e. we shift the solution array for the     *
 * next cycle posteriorly by one nucleus); if on the other hand, the most  *
 * posterior nucleus lies outside our new array, we just forget about it   */

    else if ( what2do[i] & DIVIDE ) {

      lin = TheMaternal.GetStartLin(solution.array[i+1].time);
      for (j=0; j < solution.array[i].state.size; j++) {

    k  = j % defs.ngenes;     /* k: index of gene k in current nucleus */
    ap = j / defs.ngenes;      /* ap: rel. nucleus position on AP axis */

/* evaluate ii: index of anterior daughter nucleus */

    if ( lin % 2 )
      ii = 2 * ap * defs.ngenes + k - defs.ngenes;
    else
      ii = 2 * ap * defs.ngenes + k;

/* skip the first most anterior daughter nucleus in case lin is odd */

        if ( ii >= 0 )
      solution.array[i+1].state.array[ii] =
        solution.array[i].state.array[j];

/* the second daughter only exists if it is still within the region */

    if ( ii + defs.ngenes < solution.array[i+1].state.size )
      solution.array[i+1].state.array[ii+defs.ngenes] =
        solution.array[i].state.array[j];
      }

/* Divide the history of the delay solver */

      if (delay_solver)
          delay_solver->DivideHistory(solution.array[i].time,
                                      solution.array[i + 1].time);

      if ( debug )
          fprintf(slog, "Blastoderm: nuclear division %d -> %d nuclei.\n",
                  solution.array[i].state.size / defs.ngenes,
                  solution.array[i+1].state.size / defs.ngenes);

    }

/* This is the MITOTATE rule, which advances the blastoderm from
transitions[i] to transitions[i] + epsilon, from where we propagate
again. This is there to ensure that rounding errors from the solver do
not lead to incorrect rules being used at the end-points */

    else if ( what2do[i] & MITOTATE ) {
/*      printf("aAAHHH! in MITOTATE, t=%.16f, t+epilon=%.16f\n",
                        solution.array[i].time,
                                solution.array[i+1].time ); */
      for (j=0; j < solution.array[i].state.size; j++)
        solution.array[i+1].state.array[j] =
                                    solution.array[i].state.array[j];


    }

/* In case we have to PROPAGATE the differential equeations, we call the   *
 * solver; you have to make sure that the appropriate rule has been set in *
 * zygotic.c (use SetRule(), rules are MITOSIS or INTERPHASE), otherwise   *
 * you will propagate the wrong equations; the solver needs pointers to    *
 * arrays of concentration for start and end time as well as those start   *
 * and end times themselves; all solvers need a suggested stepsize, adap-  *
 * tive stepsize solvers need the accuracy argument; we also need the accu-*
 * racy argument for global stepsize control (if ever implemented); lastly *
 * we need to tell the solver how big input and output arrays are          */

    else if ( what2do[i] & PROPAGATE )
      solve->ps(solution.array[i].state.array,
        solution.array[i+1].state.array,
        solution.array[i].time,
        solution.array[i+1].time,
        stephint, accuracy,
        solution.array[i].state.size,
        slog);

/* unknown op? -> error! */

    else
      error("op was %d!?", what2do[i]);
  }

/* After having calculated the solution, free the mutant parameter structs *
 * and return the result                                                   */

  FreeMutant();
  if (delay_solver)
      delay_solver->resetSolver();
  return solution;

}

/*** PrintBlastoderm: writes the output of the model to a stream specified *
 *                    by the fp file pointer; the table is a solution of   *
 *                    the model as returned by Blastoderm(), the id speci- *
 *                    fies the title of the output and ndigits specifies   *
 *                    the floating point precision of the concentrations   *
 *                    to be printed                                        *
 ***************************************************************************/
void zygotic::PrintBlastoderm(FILE *fp, NArrPtr table, char *id, int ndigits,
                     int columns)
{
    int                i, j, k;                       /* local loop counters */
    int                lineage;                /* lineage number for nucleus */

    /* print title (id) */

    fprintf(fp, "$%s\n", id);

    /* print table with correct lineage numbers (obtained from maternal.c) */

    for (i=0; i < table.size; i++) {
        for(j=0; j < (table.array[i].state.size/columns); j++) {
            lineage = TheMaternal.GetStartLin(table.array[i].time) + j;
            fprintf(fp, "%5d %9.3f", lineage, table.array[i].time);
            for (k=0; k < columns; k++)
                fprintf(fp, " %*.*f", ndigits+5, ndigits,
                        table.array[i].state.array[k+(j*columns)]);
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }
    fprintf(fp,"$$\n");
    fflush(fp);
}

void zygotic::d_deriv(double *v, double **vd, double t, double *vdot, int n)
{
    double vinput1 = 0;
    int        m;                                        /* number of nuclei */
    int        ap;                 /* nuclear index on AP axis [0,1,...,m-1] */
    int        i,j;                                   /* local loop counters */
    int        k;                   /* index of gene k in a specific nucleus */
    int        base, base1;            /* index of first gene in a specific nucleus */
    int        *l_rule;       /* for autonomous implementation */
#ifdef ALPHA_DU
    int        incx = 1;        /* increment step size for vsqrt input array */
    int        incy = 1;       /* increment step size for vsqrt output array */
#endif

    static int num_nucs  = 0;      /* store the number of nucs for next step */
    static int bcd_index = 0;   /* the *next* array in bicoid struct for bcd */
    static DArrPtr bcd;              /* pointer to appropriate bicoid struct */
    double **v_ext;           /* array to hold the external input
                                      concentrations at time t */
    SoDe *delay_solver = get_delay_solver(); // This may need to change to the
                                             // parent of the delay solvers if
                                             // we have multiple delay solvers
    assert(delay_solver); // we need this for external input

    /* get D parameters and bicoid gradient according to cleavage cycle */

    m = n / defs.ngenes;                        /* m is the number of nuclei */
    if (m != num_nucs) {      /* time-varying quantities only vary by ccycle */
        TheMaternal.GetD(t, lparm.d, defs.diff_schedule, D);
        /* get diff coefficients, according to diff schedule */
        if (num_nucs > m)                  /* started a new iteration in score */
            bcd_index = 0;                           /* -> start all over again! */
        num_nucs = m;                         /* store # of nucs for next step */
        bcd = TheMaternal.GetBicoid(t, genindex);                   /* get bicoid gradient */
        if( bcd.size != num_nucs)
            error("DvdtDelay: %d nuclei don't match Bicoid!", num_nucs);
        bcd_index++;                   /* store index for next bicoid gradient */
    }

    l_rule = (int *) calloc(defs.ngenes, sizeof(int));

    for (i = 0; i < defs.ngenes; i++)
        l_rule[i] = !TheMaternal.Theta(t - lparm.tau[i]);/*for autonomous equations */

    /* Here we retrieve the external input concentrations into v_ext */
    /* 01/13/10 Manu: If there are no external inputs, don't
     * allocate v_ext or populate it with external input
     * concentrations */


    if (defs.egenes > 0) {
        v_ext = (double **) calloc(defs.ngenes, sizeof(double *));

        for (i = 0; i < defs.ngenes; i++) {

            v_ext[i] = (double *) calloc(m*defs.egenes, sizeof(double));
            delay_solver->ExternalInputs(t - lparm.tau[i], t, v_ext[i], m*defs.egenes);
        }
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
    if ( gofu == Sqrt ) {

        for (base = 0, base1=0, ap=0; base < n ; base += defs.ngenes, base1 += defs.egenes, ap++) {
            for (i=base; i < base + defs.ngenes; i++) {

                k = i - base;

                vinput1  = lparm.h[k];
                vinput1 += lparm.m[k] * bcd.array[ap];    /* ap is nuclear index */

                for(j=0; j < defs.egenes; j++)
                    vinput1 += lparm.E[(k*defs.egenes)+j] * v_ext[k][base1 + j];

                for(j=0; j < defs.ngenes; j++)
                    vinput1 += lparm.T[(k*defs.ngenes)+j] * vd[k][base + j];

                bot2[i]   = 1 + vinput1 * vinput1;
                vinput[i] = vinput1;
            }
        }

        /* now calculate sqrt(1+u2); store it in bot[] */
#ifdef ALPHA_DU
        vsqrt_(bot2,&incx,bot,&incy,&n);    /* superfast DEC vector function */
#else
        for(i=0; i < n; i++)                  /* slow traditional style sqrt */
            bot[i] = sqrt(bot2[i]);
#endif
        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if ( n == defs.ngenes ) {       /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for( base=0; base<n; base+=defs.ngenes ) {
                for( i=base; i<base+defs.ngenes; i++ ) {

                    k       = i - base;
                    vdot1   = -lparm.lambda[k] * v[i];
                    g1      = 1 + vinput[i]/bot[i];
                    vdot1  += l_rule[k] * lparm.R[k] * 0.5 * g1;
                    vdot[i] = vdot1;
                }
            }

        } else {                    /* then for multiple nuclei -> diffusion */

            register double vdot1,g1;

            for(i=0; i < defs.ngenes; i++){     /* first anterior-most nucleus */
                k       = i;
                vdot1   = -lparm.lambda[k] * v[i];
                g1      = 1 + vinput[i]/bot[i];
                vdot1  += l_rule[k] * lparm.R[k] * 0.5 * g1;
                vdot1  += D[i] * (v[i + defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for(base = defs.ngenes; base < n - defs.ngenes; base += defs.ngenes){
                for(i=base; i < base + defs.ngenes; i++){
                    k       = i - base;
                    vdot1   = -lparm.lambda[k] * v[i];
                    g1      =  1 + vinput[i]/bot[i];
                    vdot1  += l_rule[k] * lparm.R[k] * 0.5 * g1;
                    vdot1  += D[k] * ((v[i - defs.ngenes] - v[i]) +
                            (v[i + defs.ngenes] - v[i]) );
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for(i=base; i < base + defs.ngenes; i++){
                k       = i - base;
                vdot1   = -lparm.lambda[k] * v[i];
                g1      =  1 + vinput[i]/bot[i];
                vdot1  += l_rule[k] * lparm.R[k] * 0.5 * g1;
                vdot1  += D[k] * (v[i - defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
        }

        /***************************************************************************
         *                                                                         *
         *        g(u) = 1/2 * (tanh(u) + 1) )                                     *
         *                                                                         *
         ***************************************************************************/

    } else if ( gofu == Tanh ) {

        for (base = 0, base1=0, ap=0; base < n ; base += defs.ngenes, base1 += defs.egenes, ap++) {
            for(i=base; i < base + defs.ngenes; i++){

                k = i - base;

                vinput1  = lparm.h[k];
                vinput1 += lparm.m[k] * bcd.array[ap];    /* ap is nuclear index */

                for(j=0; j < defs.egenes; j++)
                    vinput1 += lparm.E[(k*defs.egenes)+j] * v_ext[k][base1 + j];

                for(j=0; j < defs.ngenes; j++)
                    vinput1  += lparm.T[(k*defs.ngenes)+j] * vd[k][base + j];

                vinput[i] = vinput1;
            }
        }

        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if ( n == defs.ngenes ) {       /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for( base=0; base<n; base+=defs.ngenes ) {
                for( i=base; i<base+defs.ngenes; i++ ) {

                    k       = i - base;
                    vdot1   = -lparm.lambda[k] * v[i];
                    g1      = tanh(vinput[i]) + 1;
                    vdot1  += l_rule[k]*lparm.R[k] * 0.5 * g1;
                    vdot[i] = vdot1;
                }
            }

        } else {                    /* then for multiple nuclei -> diffusion */

            register double vdot1,g1;

            for(i=0; i < defs.ngenes; i++){     /* first anterior-most nucleus */
                k       = i;
                vdot1   = -lparm.lambda[k] * v[i];
                g1      = tanh(vinput[i]) + 1;
                vdot1  += l_rule[k]*lparm.R[k] * 0.5 * g1;
                vdot1  += D[i] * (v[i + defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for(base = defs.ngenes; base < n - defs.ngenes; base += defs.ngenes){
                for(i=base; i < base + defs.ngenes; i++){
                    k       = i - base;
                    vdot1   = -lparm.lambda[k] * v[i];
                    g1      = tanh(vinput[i]) + 1;
                    vdot1  += l_rule[k]*lparm.R[k] * 0.5 * g1;
                    vdot1  += D[k] * ((v[i - defs.ngenes] - v[i]) +
                            (v[i + defs.ngenes] - v[i]) );
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for(i=base; i < base + defs.ngenes; i++){
                k       = i - base;
                vdot1   = -lparm.lambda[k] * v[i];
                g1      = tanh(vinput[i]) + 1;
                vdot1  += l_rule[k]*lparm.R[k] * 0.5 * g1;
                vdot1  += D[k] * (v[i - defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
        }

        /***************************************************************************
         *                                                                         *
         *        g(u) = 1 / (1 + exp(-2u))                                        *
         *                                                                         *
         ***************************************************************************/

    } else if ( gofu == Exp ) {

        for (base = 0, base1=0, ap=0; base < n ; base += defs.ngenes, base1 += defs.egenes, ap++) {
            for(i=base; i < base + defs.ngenes; i++){

                k = i - base;

                vinput1  = lparm.h[k];
                vinput1 += lparm.m[k] * bcd.array[ap];    /* ap is nuclear index */

                for(j=0; j < defs.egenes; j++)
                    vinput1 += lparm.E[(k*defs.egenes)+j] * v_ext[k][base1 + j];

                for(j=0;  j < defs.ngenes; j++)
                    vinput1  += lparm.T[(k*defs.ngenes)+j] * vd[k][base + j];

                vinput[i] = -2.0 * vinput1;
            }
        }

        /* now calculate exp(-u); store it in bot[] */
#ifdef ALPHA_DU
        vexp_(vinput,&incx,bot,&incy,&n);   /* superfast DEC vector function */
#else
        for (i=0; i<n; i++)                    /* slow traditional style exp */
            bot[i] = exp(vinput[i]);
#endif

        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if ( n == defs.ngenes ) {       /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for( base=0; base<n; base+=defs.ngenes ) {
                for( i=base; i<base+defs.ngenes; i++ ) {

                    k       = i - base;
                    vdot1   = -lparm.lambda[k] * v[i];
                    g1      = 1 / (1 + bot[i]);
                    vdot1  += l_rule[k]*lparm.R[k] * g1;
                    vdot[i] = vdot1;
                }
            }

        } else {                    /* then for multiple nuclei -> diffusion */

            register double vdot1,g1;

            for(i=0; i < defs.ngenes; i++){     /* first anterior-most nucleus */
                k       = i;
                vdot1   = -lparm.lambda[k] * v[i];
                g1      = 1 / (1 + bot[i]);
                vdot1  += l_rule[k]*lparm.R[k] * g1;
                vdot1  += D[i] * (v[i + defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for(base = defs.ngenes; base < n - defs.ngenes; base += defs.ngenes){
                for(i=base; i < base + defs.ngenes; i++){
                    k       = i - base;
                    vdot1   = -lparm.lambda[k] * v[i];
                    g1      = 1 / (1 + bot[i]);
                    vdot1  += l_rule[k]*lparm.R[k] * g1;
                    vdot1  += D[k] * ((v[i - defs.ngenes] - v[i]) +
                            (v[i + defs.ngenes] - v[i]) );
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for(i=base; i < base + defs.ngenes; i++){
                k       = i - base;
                vdot1   = -lparm.lambda[k] * v[i];
                g1      = 1 / (1 + bot[i]);
                vdot1  += l_rule[k]*lparm.R[k]  * g1;
                vdot1  += D[k] * (v[i - defs.ngenes] - v[i]);
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

    } else if ( gofu == Hvs ) {

        for (base = 0, base1=0, ap=0; base < n ; base += defs.ngenes, base1 += defs.egenes, ap++) {
            for(i=base; i < base + defs.ngenes; i++){

                k = i - base;

                vinput1  = lparm.h[k];
                vinput1 += lparm.m[k] * bcd.array[ap];    /* ap is nuclear index */

                for(j=0; j < defs.egenes; j++)
                    vinput1 += lparm.E[(k*defs.egenes)+j] * v_ext[k][base1 + j];

                for(j=0; j < defs.ngenes; j++)
                    vinput1  += lparm.T[(k*defs.ngenes)+j] * vd[k][base + j];
                vinput[i] = vinput1;
            }
        }

        /* next loop does the rest of the equation (R, Ds and lambdas) */
        /* store result in vdot[] */

        if ( n == defs.ngenes ) {       /* first part: one nuc, no diffusion */

            register double vdot1, g1;

            for( base=0; base<n; base+=defs.ngenes ) {
                for( i=base; i<base+defs.ngenes; i++ ) {

                    k       = i - base;
                    vdot1   = -lparm.lambda[k] * v[i];
                    if (vinput[i] >= 0.)
                        g1 = 1.;
                    else
                        g1 = 0.;
                    vdot1  += l_rule[k]*lparm.R[k] * g1;
                    vdot[i] = vdot1;
                }
            }

        } else {                    /* then for multiple nuclei -> diffusion */

            register double vdot1,g1;

            for(i=0; i < defs.ngenes; i++){     /* first anterior-most nucleus */
                k       = i;
                vdot1   = -lparm.lambda[k] * v[i];
                if (vinput[i] >= 0.)
                    g1 = 1.;
                else
                    g1 = 0.;
                vdot1  += l_rule[k]*lparm.R[k] * g1;
                vdot1  += D[i] * (v[i + defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
            /* then middle nuclei */
            for(base = defs.ngenes; base < n - defs.ngenes; base += defs.ngenes){
                for(i=base; i < base + defs.ngenes; i++){
                    k       = i - base;
                    vdot1   = -lparm.lambda[k] * v[i];
                    if (vinput[i] >= 0.)
                        g1 = 1.;
                    else
                        g1 = 0.;
                    vdot1  += l_rule[k]*lparm.R[k] * g1;
                    vdot1  += D[k] * ((v[i - defs.ngenes] - v[i]) +
                            (v[i + defs.ngenes] - v[i]) );
                    vdot[i] = vdot1;
                }
            }
            /* last: posterior-most nucleus */
            for(i=base; i < base + defs.ngenes; i++){
                k       = i - base;
                vdot1   = -lparm.lambda[k] * v[i];
                if (vinput[i] >= 0.)
                    g1 = 1.;
                else
                    g1 = 0.;
                vdot1  += l_rule[k]*lparm.R[k] * g1;
                vdot1  += D[k] * (v[i - defs.ngenes] - v[i]);
                vdot[i] = vdot1;
            }
        }

    } else
        error("DvdtDelay: unknown g(u)");

    /* during mitosis only diffusion and decay happen */

    free(l_rule);

    if (defs.egenes > 0) {
        for (i = 0; i < defs.ngenes; i++)
            free(v_ext[i]);

        free(v_ext);
    }

    return;
}

/*** PrintParameters: prints an eqparms section with 'title' to the stream *
 *                    indicated by fp                                      *
 ***************************************************************************/

void zygotic::PrintParameters(FILE *fp, EqParms *p, const char *title, int ndigits)
{
  int    i, j;                                      /* local loop counters */
  double lambda_tmp;                           /* temporary var for lambda */


  fprintf(fp, "$%s\n", title);
  fprintf(fp, "promoter_strengths:\n");             /* Rs are written here */

  for ( i=0; i<defs.ngenes; i++ )
    fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->R[i]);

  fprintf(fp, "\n");
  fprintf(fp, "genetic_interconnect_matrix:\n");        /* Ts written here */

  for ( i=0; i<defs.ngenes; i++ ) {
    for ( j=0; j<defs.ngenes; j++ )
      fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->T[(i*defs.ngenes)+j]);
    fprintf(fp, "\n");
  }

  fprintf(fp, "external_input_strengths:\n");        /* Es written here */

  for ( i=0; i<defs.ngenes; i++ ) {
    for ( j=0; j<defs.egenes; j++ )
      fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->E[(i*defs.egenes)+j]);
    fprintf(fp, "\n");
  }

  fprintf(fp, "maternal_connection_strengths:\n");      /* ms written here */

  for ( i=0; i<defs.ngenes; i++ )
    fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->m[i]);

  fprintf(fp, "\n");
  fprintf(fp, "promoter_thresholds:\n");            /* hs are written here */

  for ( i=0; i<defs.ngenes; i++ )
    fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->h[i]);

  fprintf(fp, "\n");
  fprintf(fp, "diffusion_parameter(s):\n");         /* ds are written here */

  if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') )
    fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->d[0]);
  else
    for ( i=0; i<defs.ngenes; i++ )
      fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->d[i]);

  fprintf(fp, "\n");
  fprintf(fp, "protein_half_lives:\n");        /* lambdas are written here */

  for ( i=0; i<defs.ngenes; i++ ) {
    lambda_tmp = log(2.) / p->lambda[i];           /* conversion done here */
    fprintf(fp, "%*.*f ", ndigits+4, ndigits, lambda_tmp);
  }

  fprintf(fp, "\n");
  fprintf(fp, "translational_transcriptional_delays:\n");        /* taus are written here */

  for ( i=0; i<defs.ngenes; i++ )
    fprintf(fp, "%*.*f ", ndigits+4, ndigits, p->tau[i]);

  fprintf(fp, "\n$$\n");
  fflush(fp);

}


