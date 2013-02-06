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

