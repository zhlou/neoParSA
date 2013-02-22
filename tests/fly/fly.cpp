/*
 * fly.cpp
 *
 *  Created on: Feb 18, 2013
 *      Author: zhlou
 */
#include <stdexcept>
#include "fly.h"

using namespace std;

fly_params readFlyParams(xmlNode *docroot)
{
    static const double MAX_STEPSIZE = 100.8;
    fly_params params;

    // put defaults on
    params.section_title = string("eqparms");
    params.solver_name = string("Rk4");
    params.ndigits = 12;
    params.gutndigits = 6;
    params.olddivstyle = 0;
    // params.penaltyflag = 0;
    // params.rmsflag = 1;
    params.GofU = Sqrt;
    params.gutflag = 0;
    params.stepsize = 1.;
    params.accuracy = 0.001;
    params.infile = NULL;
    params.dumpptr = NULL;
    params.slog = NULL;
    params.debug = 0;

    // now start to read from xml
    xmlNode *section = docroot->children;
    while (section != NULL) {
        if(!xmlStrcmp(section->name,(const xmlChar *)"fly"))
            break;
        section = section->next;
    }
    if (section == NULL)
        throw runtime_error("No fly section specified");

    // infile must be specified
    xmlChar *prop = NULL;
    prop = xmlGetProp(section, (const xmlChar *)"infile");
    if (prop == NULL)
        throw runtime_error("No input file specified");
    params.infile_name = (char *)prop;
    xmlFree(prop);
    prop = NULL;
    params.infile = fopen(params.infile_name.c_str(), "r");
    if (params.infile == NULL)
        throw runtime_error(string("Unable to open file ")+params.infile_name);


    // rest of the parameters are optional
    prop = xmlGetProp(section, (const xmlChar *)"solver");
    if (prop != NULL) {
        params.solver_name = (char *)prop;
        xmlFree(prop);
        prop = NULL;
    }
    prop = xmlGetProp(section, (const xmlChar *)"ndigits");
    if (prop != NULL) {
        params.ndigits = atoi((char *)prop);
        xmlFree(prop);
        prop = NULL;
        if (params.ndigits < 0)
            throw runtime_error("ndigits: what exactly would a negative precision be???");
        params.gutndigits = params.ndigits;
    }
    prop = xmlGetProp(section, (const xmlChar *)"gutflag");
    if (prop != NULL) {
        params.gutflag = atoi((char *)prop);
        xmlFree(prop);
        prop = NULL;
    }
    prop = xmlGetProp(section, (const xmlChar *)"olddivstyle");
    if (prop != NULL) {
        params.olddivstyle = atoi((char *)prop);
        xmlFree(prop);
        prop = NULL;
    }
    prop = xmlGetProp(section, (const xmlChar *)"debug");
    if (prop != NULL) {
        params.debug = atoi((char *)prop);
        xmlFree(prop);
        prop = NULL;
        if (params. debug) {
            string slogfile(infile_name + ".slog");
            params.slog = fopen(slogfile.c_str(), "w");
        }
    }
    prop = xmlGetProp(section, (const xmlChar *)"gofu");
    if (prop != NULL) {
        // The C++ string *copies* the prop data to the string structure. So
        // prop can be freed immediately. Since the string is declared here,
        // the memory space it allocated will be freed as soon as the code
        // leaves this if branch.
        string in_gofu((char *)prop);
        xmlFree(prop);
        prop = NULL;
        if (in_gofu == "Sqrt")
            params.GofU = Sqrt;
        else if (in_gofu == "Tanh")
            params.GofU = Tanh;
        else if (in_gofu == "Exp")
            params.GofU = Exp;
        else if (in_gofu == "Hvs")
            params.GofU = Hvs;
        else
            throw runtime_error(string("Unrecognized GofU: ") + in_gofu);
    }
    prop = xmlGetProp(section, (const xmlChar *)"stepsize");
    if (prop != NULL) {
        params.stepsize = strtod((char *)prop, NULL);
        xmlFree(prop);
        prop = NULL;
        if (params.stepsize <= 0)
            throw runtime_error("stepsize is too small");
        if (params.stepsize > MAX_STEPSIZE)
            throw runtime_error("stepsize is too big");
    }
    prop = xmlGetProp(section, (const xmlChar *)"accuracy");
    if (prop != NULL) {
        params.accuracy = strtod((char *)prop, NULL);
        xmlFree(prop);
        prop = NULL;
        if (params.stepsize <= 0)
            throw runtime_error("accuracy is too small");
    }
    prop = xmlGetProp(section, (const xmlChar *)"section");
    if (prop != NULL) {
        params.section_title = (char *)prop;
        xmlFree(prop);
        prop = NULL;
    }
    return params;
}

fly::fly(const fly_params &params) :
        TheMaternal(params.infile, params.olddivstyle),
        defs(TheMaternal.getProblem()),
        zygote(TheMaternal, params.infile, params.section_title.c_str(),
               params.GofU, params.debug, params.solver_name.c_str()),
        score(params.infile, *zygote, params.gutflag, params.gutndigits,
              params.stepsize, params.accuracy, params.slog,
              params.infile_name.c_str(), params.debug)
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

double fly::get_score()
{
    return score.Score();
}
