/*
 * fly.cpp
 *
 *  Created on: Feb 18, 2013
 *      Author: zhlou
 */

#include "fly.h"

fly_params::fly_params() : // initialize default parameters
        section_title("eqparms"), solver_name("Rk4")
{
    ndigits = 12;
    gutndigits = 6;
    penaltyflag = 0;
    rmsflag = 1;
    gutflag = 0;
    stepsize = 1.;
    accuracy = 0.001;
    infile = NULL;
    dumpptr = NULL;
    slog = NULL;
    debug = 0;
}

fly::fly(const fly_params &params) :
        TheMaternal(params.infile), defs(TheMaternal.getProblem()),
        zygote(TheMaternal, params.infile, params.section_title.c_str(), params.debug,
               params.solver_name.c_str()),
        score(params.infile, *zygote, params.gutflag, params.gutndigits, params.stepsize,
              params.accuracy,params.slog, params.infile_name.c_str(),
              params.debug)
{
    // read tweak

    int *temptweak, *temptweak1; /* temporary arrays to read tweaks */

    int i; /* local loop counter */
    int c; /* used to parse text lines */
    int linecount = 0; /* keep track of # of lines read */
    int Tcount = 0; /* keep track of T lines read */
    int Ecount = 0; /* keep track of E lines read */

    char *base; /* pointer to beginning of line string */
    char *record; /* string for reading whole line of params */

    char **fmt; /* array of format strings for reading params */
    char **fmt1; /* array of format strings for reading
     E tweaks */
    char *skip, *skip1; /* string of values to be skipped */

    const char read_fmt[] = "%d"; /* read an int */
    const char skip_fmt[] = "%*d "; /* ignore an int */

    base = (char *) calloc(MAX_RECORD, sizeof(char *));

    skip = (char *) calloc(MAX_RECORD, sizeof(char *));

    skip1 = (char *) calloc(MAX_RECORD, sizeof(char *));

    fmt = (char **) calloc(defs.ngenes, sizeof(char *));
    if (defs.egenes > 0)
        fmt1 = (char **) calloc(defs.egenes, sizeof(char *));

    temptweak = (int *) calloc(defs.ngenes, sizeof(int *));
    temptweak1 = (int *) calloc(defs.egenes, sizeof(int *));

    /* create format strings according to the number of genes */

    for (i = 0; i < defs.ngenes; i++) {
        fmt[i] = (char *) calloc(MAX_RECORD, sizeof(char));
        fmt[i] = strcpy(fmt[i], skip);
        fmt[i] = strcat(fmt[i], read_fmt);
        skip = strcat(skip, skip_fmt);
    }

    /* create format strings according to the number of external inputs */

    if (defs.egenes > 0) {

        for (i = 0; i < defs.egenes; i++) {
            fmt1[i] = (char *) calloc(MAX_RECORD, sizeof(char));
            fmt1[i] = strcpy(fmt1[i], skip1);
            fmt1[i] = strcat(fmt1[i], read_fmt);
            skip1 = strcat(skip1, skip_fmt);
        }

    }

    /* initialize the Tweak struct */

    tweak.Rtweak = (int *) calloc(defs.ngenes, sizeof(int));
    tweak.Ttweak = (int *) calloc(defs.ngenes * defs.ngenes, sizeof(int));

    if (defs.egenes > 0)
        tweak.Etweak = (int *) calloc(defs.ngenes * defs.egenes, sizeof(int));

    tweak.mtweak = (int *) calloc(defs.ngenes, sizeof(int));
    tweak.htweak = (int *) calloc(defs.ngenes, sizeof(int));
    if ((defs.diff_schedule == 'A') || (defs.diff_schedule == 'C')) {
        tweak.dtweak = (int *) malloc(sizeof(int));
    } else {
        tweak.dtweak = (int *) calloc(defs.ngenes, sizeof(int));
    }
    tweak.lambdatweak = (int *) calloc(defs.ngenes, sizeof(int));
    tweak.tautweak = (int *) calloc(defs.ngenes, sizeof(int));

    params.infile = FindSection(params.infile, "tweak"); /* find tweak section */
    if (!params.infile)
        error("ReadTweak: could not locate tweak\n");

    while (strncmp((base = fgets(base, MAX_RECORD, params.infile)), "$$", 2)) {

        record = base;

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

                if ((linecount == 5) && ((defs.diff_schedule == 'A')
                        || (defs.diff_schedule == 'C'))) {
                    if (1 != sscanf(record, fmt[0], &temptweak[0]))
                        error("ReadTweak: error reading tweaks");
                } else if (linecount == 2) {
                    for (i = 0; i < defs.egenes; i++) {
                        if (1 != sscanf(record, fmt1[i], &temptweak1[i]))
                            error("ReadTweak: error reading tweak variables");
                    }
                } else {
                    for (i = 0; i < defs.ngenes; i++) {
                        if (1 != sscanf(record, fmt[i], &temptweak[i]))
                            error("ReadTweak: error reading tweak variables");
                    }
                }

                switch (linecount) { /* copy read parameters into the right array */
                case 0:
                    for (i = 0; i < defs.ngenes; i++) /* R tweaks */
                        tweak.Rtweak[i] = temptweak[i];
                    linecount++;
                    break;
                case 1: /* T tweaks: keep track of read lines with Tcount */
                    for (i = 0; i < defs.ngenes; i++)
                        tweak.Ttweak[i + Tcount * defs.ngenes] = temptweak[i];
                    Tcount++;
                    if (Tcount == defs.ngenes)
                        linecount++;
                    break;
                case 2: /* E tweaks: keep track of read lines with Ecount */
                    if (defs.egenes > 0) {
                        for (i = 0; i < defs.egenes; i++)
                            tweak.Etweak[i + Ecount * defs.egenes] =
                                    temptweak1[i];
                        Ecount++;
                        if (Ecount == defs.ngenes)
                            linecount++;
                    }

                    break;
                case 3: /* m tweaks */
                    for (i = 0; i < defs.ngenes; i++)
                        tweak.mtweak[i] = temptweak[i];
                    linecount++;
                    break;
                case 4:
                    for (i = 0; i < defs.ngenes; i++) /* h tweaks */
                        tweak.htweak[i] = temptweak[i];
                    linecount++;
                    break;
                case 5: /* d tweaks: consider diff. schedule */
                    if ((defs.diff_schedule == 'A') || (defs.diff_schedule
                            == 'C')) {
                        tweak.dtweak[0] = temptweak[0];
                    } else {
                        for (i = 0; i < defs.ngenes; i++)
                            tweak.dtweak[i] = temptweak[i];
                    }
                    linecount++;
                    break;
                case 6: /* lambda tweaks */
                    for (i = 0; i < defs.ngenes; i++)
                        tweak.lambdatweak[i] = temptweak[i];
                    linecount++;
                    break;
                case 7: /* lambda tweaks */
                    for (i = 0; i < defs.ngenes; i++)
                        tweak.tautweak[i] = temptweak[i];
                    linecount++;
                    break;
                default:
                    error("ReadTweak: too many data lines in tweak section");
                }
                break; /* don't do rest of loop anymore! */
            }

            else if (isspace(c)) { /* ignore leading white space */
                c = (int) *(++record);
            }

            else { /* anything but space or digit means comment */
                break;
            }
        }
    }

    free(temptweak);
    free(temptweak1);
    free(base);
    free(skip);
    free(skip1);

    for (i = 0; i < defs.ngenes; i++)
        free(fmt[i]);
    free(fmt);

    if (defs.egenes > 0) {
        for (i = 0; i < defs.egenes; i++)
            free(fmt1[i]);
        free(fmt1);
    }
}
