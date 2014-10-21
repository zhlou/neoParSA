/*
 * fly.cpp
 *
 *  Created on: Feb 18, 2013
 *      Author: zhlou
 */
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include "fly.h"

using namespace std;

fly_params readFlyParams(xmlNode *docroot, const char* default_section)
{
    static const double MAX_STEPSIZE = 100.8;
    fly_params params;

    // put defaults on
    params.section_title = default_section;
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
            string slogfile(params.infile_name + ".slog");
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
    prop = xmlGetProp(section, (const xmlChar *)"outfile");
    if (prop != NULL) {
        params.outfile_name = (char *)prop;
        xmlFree(prop);
        prop = NULL;
    }
    return params;
}

fly::Tweak fly::ReadTweak(const fly_params& params)
{
    Tweak l_tweak;
    int *temptweak, *temptweak1; /* temporary arrays to read tweaks */
    int i; /* local loop counter */
    int c; /* used to parse text lines */
    int linecount = 0; /* keep track of # of lines read */
    int Tcount = 0; /* keep track of T lines read */
    int Ecount = 0; /* keep track of E lines read */
    char* base; /* pointer to beginning of line string */
    char* record; /* string for reading whole line of params */
    char** fmt; /* array of format strings for reading params */
    char** fmt1; /* array of format strings for reading
     E tweaks */
    char *skip, *skip1; /* string of values to be skipped */
    const char read_fmt[] = "%d"; /* read an int */
    const char skip_fmt[] = "%*d "; /* ignore an int */
    FILE* fp = params.infile;
    base = (char*) (calloc(MAX_RECORD, sizeof(char*)));
    skip = (char*) (calloc(MAX_RECORD, sizeof(char*)));
    skip1 = (char*) (calloc(MAX_RECORD, sizeof(char*)));
    fmt = (char**) (calloc(defs.ngenes, sizeof(char*)));
    if (defs.egenes > 0)
        fmt1 = (char**) (calloc(defs.egenes, sizeof(char*)));

    temptweak = (int*) (calloc(defs.ngenes, sizeof(int*)));
    temptweak1 = (int*) (calloc(defs.egenes, sizeof(int*)));
    /* create format strings according to the number of genes */
    for (i = 0; i < defs.ngenes; i++) {
        fmt[i] = (char*) (calloc(MAX_RECORD, sizeof(char)));
        fmt[i] = strcpy(fmt[i], skip);
        fmt[i] = strcat(fmt[i], read_fmt);
        skip = strcat(skip, skip_fmt);
    }
    /* create format strings according to the number of external inputs */
    if (defs.egenes > 0) {
        for (i = 0; i < defs.egenes; i++) {
            fmt1[i] = (char*) (calloc(MAX_RECORD, sizeof(char)));
            fmt1[i] = strcpy(fmt1[i], skip1);
            fmt1[i] = strcat(fmt1[i], read_fmt);
            skip1 = strcat(skip1, skip_fmt);
        }
    }
    /* initialize the Tweak struct */
    l_tweak.Rtweak = (int*) (calloc(defs.ngenes, sizeof(int)));
    l_tweak.Ttweak = (int*) (calloc(defs.ngenes * defs.ngenes, sizeof(int)));
    if (defs.egenes > 0)
        l_tweak.Etweak = (int*) (calloc(defs.ngenes * defs.egenes, sizeof(int)));

    l_tweak.mtweak = (int*) (calloc(defs.ngenes, sizeof(int)));
    l_tweak.htweak = (int*) (calloc(defs.ngenes, sizeof(int)));
    if ((defs.diff_schedule == 'A') || (defs.diff_schedule == 'C')) {
        l_tweak.dtweak = (int*) (malloc(sizeof(int)));
    } else {
        l_tweak.dtweak = (int*) (calloc(defs.ngenes, sizeof(int)));
    }
    l_tweak.lambdatweak = (int*) (calloc(defs.ngenes, sizeof(int)));
    l_tweak.tautweak = (int*) (calloc(defs.ngenes, sizeof(int)));
    fp = FindSection(fp, "tweak"); /* find tweak section */
    if (!fp)
        error("ReadTweak: could not locate tweak\n");

    while (strncmp((base = fgets(base, MAX_RECORD, fp)), "$$", 2)) {
        record = base;
        c = (int) (*record);
        while (c != '\0') {
            if (isdigit(c)) {
                /* line contains data */
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

                switch (linecount) {
                /* copy read parameters into the right array */
                case 0:
                    for (i = 0; i < defs.ngenes; i++)
                        /* R tweaks */
                        l_tweak.Rtweak[i] = temptweak[i];
                    linecount++;
                    break;
                case 1: /* T tweaks: keep track of read lines with Tcount */
                    for (i = 0; i < defs.ngenes; i++)
                        l_tweak.Ttweak[i + Tcount * defs.ngenes] = temptweak[i];
                    Tcount++;
                    if (Tcount == defs.ngenes)
                        linecount++;

                    break;
                case 2: /* E tweaks: keep track of read lines with Ecount */
                    if (defs.egenes > 0) {
                        for (i = 0; i < defs.egenes; i++)
                            l_tweak.Etweak[i + Ecount * defs.egenes] =
                                    temptweak1[i];
                        Ecount++;
                        if (Ecount == defs.ngenes)
                            linecount++;
                    }
                    break;
                case 3: /* m tweaks */
                    for (i = 0; i < defs.ngenes; i++)
                        l_tweak.mtweak[i] = temptweak[i];
                    linecount++;
                    break;
                case 4:
                    for (i = 0; i < defs.ngenes; i++)
                        /* h tweaks */
                        l_tweak.htweak[i] = temptweak[i];
                    linecount++;
                    break;
                case 5: /* d tweaks: consider diff. schedule */
                    if ((defs.diff_schedule == 'A') || (defs.diff_schedule
                            == 'C')) {
                        l_tweak.dtweak[0] = temptweak[0];
                    } else {
                        for (i = 0; i < defs.ngenes; i++)
                            l_tweak.dtweak[i] = temptweak[i];
                    }
                    linecount++;
                    break;
                case 6: /* lambda tweaks */
                    for (i = 0; i < defs.ngenes; i++)
                        l_tweak.lambdatweak[i] = temptweak[i];
                    linecount++;
                    break;
                case 7: /* lambda tweaks */
                    for (i = 0; i < defs.ngenes; i++)
                        l_tweak.tautweak[i] = temptweak[i];
                    linecount++;
                    break;
                default:
                    error("ReadTweak: too many data lines in tweak section");
                }
                break; /* don't do rest of loop anymore! */
            } else if (isspace(c)) {
                /* ignore leading white space */
                c = (int) (*(++record));
            } else {
                /* anything but space or digit means comment */
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
    return l_tweak;
}

fly::fly(const fly_params &params) :
        TheMaternal(params.infile, params.olddivstyle),
        defs(TheMaternal.getProblem()),
        zygote(TheMaternal, params.infile, params.section_title.c_str(),
               params.GofU, params.debug, params.solver_name.c_str()),
        score(params.infile, zygote, params.gutflag, params.gutndigits,
              params.stepsize, params.accuracy, params.slog,
              params.infile_name.c_str(), params.debug)
{
    infile = params.infile_name;
    outfile = params.outfile_name;
    ndigits = params.ndigits;
    // read tweak

    tweak = ReadTweak(params);
    Translate(ptab);
    nparams = ptab.size();


    // book keeping on scores
    chisq = score.Score();
    score_valid = true;

}

double fly::get_score()
{
    if (score_valid)
        return chisq;
    else
        return updateChisq();
}

double fly::get_rms()
{
    return sqrt((get_score() - score.GetPenalty()) /
                (double) score.GetNDatapoints());
}

fly::~fly()
{
    free(tweak.Rtweak);
    free(tweak.Ttweak);
    if (defs.egenes > 0)
        free(tweak.Etweak);
    free(tweak.mtweak);
    free(tweak.htweak);
    free(tweak.dtweak);
    free(tweak.lambdatweak);
    free(tweak.tautweak);
}

void fly::Translate(vector<ParamList> &tab)
{
    ParamList p;
    int i, j;

    // new ParamList element should always non-restorable
    p.previous = 0;
    p.canRestore = false;

  /* Get limits and parameters */

    EqParms *parm   = zygote.GetParameters();   /* get pointer to EqParm struct in zygotic.c */
    SearchSpace *limits = score.GetLimits();           /* get pointer to SearchSpace in score.c */

  /* if we are using a penalty, not all of Limits will have been allocated   */
  /* and we must take care not to take indices of unallocated arrays!!!      */


    for ( i=0; i < defs.ngenes; i++ )                 /* pointers to R stuff */
        if ( tweak.Rtweak[i] == 1 ) {
            p.param = &parm->R[i];
            p.param_range = &limits->Rlim[i];
            tab.push_back(p);
        }

    for ( i=0; i < defs.ngenes; i++ )                 /* pointers to T stuff */
        for ( j=0; j < defs.ngenes; j++ )
            if ( tweak.Ttweak[i*defs.ngenes+j] == 1 ) {
                p.param = &parm->T[j+i*defs.ngenes];
                if ( limits->pen_vec == NULL )
                    p.param_range = &limits->Tlim[j+i*defs.ngenes];
                else
                    p.param_range = NULL;
                tab.push_back(p);
            } // end if

    /* 01/13/10 Manu: If there are no external inputs, the for
     * loop below should never be entered */

    for ( i=0; i < defs.ngenes; i++ )                 /* pointers to E stuff */
        for ( j=0; j < defs.egenes; j++ )
            if ( tweak.Etweak[i*defs.egenes+j] == 1 ) {
                p.param = &parm->E[j+i*defs.egenes];
                if ( limits->pen_vec == NULL )
                    p.param_range = &limits->Elim[j+i*defs.egenes];
                else
                    p.param_range = NULL;
                tab.push_back(p);
            }

    for ( i=0; i < defs.ngenes; i++ )                 /* pointers to m stuff */
        if ( tweak.mtweak[i] == 1 ) {
            p.param = &parm->m[i];
            if ( limits->pen_vec == NULL )
                p.param_range = &limits->mlim[i];
            else
                p.param_range = NULL;
            tab.push_back(p);
        }

    for ( i=0; i < defs.ngenes; i++ )                 /* pointers to h stuff */
        if ( tweak.htweak[i] == 1 ) {
            p.param = &parm->h[i];
            if ( limits->pen_vec == NULL )
                p.param_range = &limits->hlim[i];
            else
                p.param_range = NULL;
            tab.push_back(p);

        }

    for ( i=0; i < defs.ngenes; i++ )                 /* pointers to d stuff */
        if ( tweak.dtweak[i] == 1 ) {
            p.param = &parm->d[i];
            p.param_range = &limits->dlim[i];
            tab.push_back(p);
            if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') )
                break;
        }

    for ( i=0; i < defs.ngenes; i++ )            /* pointers to lambda stuff */
        if ( tweak.lambdatweak[i] == 1 ) {
            p.param = &parm->lambda[i];
            p.param_range = &limits->lambdalim[i];
            tab.push_back(p);
        }

    for ( i=0; i < defs.ngenes; i++ )            /* pointers to tau stuff */
        if ( tweak.tautweak[i] == 1 ) {
            p.param = &parm->tau[i];
            p.param_range = &limits->taulim[i];
            tab.push_back(p);
        }
}

void fly::generateMove(int idx, double delta)
{
    ptab[idx].previous = *(ptab[idx].param);
    *(ptab[idx].param) += delta;
    ptab[idx].canRestore = true;
    score_valid = false;
}

void fly::restoreMove(int idx)
{
    if(ptab[idx].canRestore) {
        *(ptab[idx].param) = ptab[idx].previous;
        ptab[idx].canRestore = false;
        score_valid = false;
    } else
        runtime_error("cannot restore move");
}

void fly::WriteParameters(const char *filename, EqParms *p, const char *title,
                          int ndigits)
{
  char   *temp;                                     /* temporary file name */
  char   *record;                         /* record to be read and written */
  char   *record_ptr;        /* pointer used to remember record for 'free' */
  char   *saverec;                 /* used to save following section title */
  char   *shell_cmd;                             /* used by 'system' below */

  FILE   *outfile;                                  /* name of output file */
  FILE   *tmpfile;                               /* name of temporary file */


  temp      = (char *)calloc(MAX_RECORD, sizeof(char));
  record    = (char *)calloc(MAX_RECORD, sizeof(char));
  saverec   = (char *)calloc(MAX_RECORD, sizeof(char));
  shell_cmd = (char *)calloc(MAX_RECORD, sizeof(char));

  record_ptr = record;            /* this is to remember record for 'free' */

/* open output and temporary file */

  outfile = fopen(filename, "r");              /* open outfile for reading */
  if ( !outfile )                              /* sorry for the confusion! */
    error("WriteParameters: error opening output file");

  temp = strcpy(temp,"parmXXXXXX");               /* required by mkstemp() */
  if ( mkstemp(temp) == -1 )              /* get unique name for temp file */
    error("WriteParameters: error creating temporary file");

  tmpfile = fopen(temp, "w");               /* ... and open it for writing */
  if ( !tmpfile )
    error("WriteParameters: error opening temporary file");

  if ( FindSection(outfile, title) ) { /* erase section if already present */
    fclose(outfile);                       /* this is a little kludgey but */
    KillSection(filename, title);      /* since KillSection needs to be in */
    outfile = fopen(filename, "r");           /* total control of the file */
  }
  rewind(outfile);

/* the follwoing two loops look for the appropriate file position to write */
/* the eqparms section (alternatives are input and eqparms)                */

  if ( !strcmp(title, "input") ) {
    while ( strncmp(record=fgets(record, MAX_RECORD, outfile),
            "$genotypes", 10) )
      fputs(record, tmpfile);
  } else if ( !strcmp(title, "eqparms") ) {
    while ( strncmp(record=fgets(record, MAX_RECORD, outfile),
            "$input", 6) )
      fputs(record, tmpfile);
  }
  fputs(record, tmpfile);

  while ( strncmp(record=fgets(record, MAX_RECORD, outfile), "$$", 2) )
    fputs(record, tmpfile);
  fputs(record, tmpfile);

  do {
    record = fgets(record, MAX_RECORD, outfile);
    if ( !record ) break;
  } while ( strncmp(record, "$", 1) );

  fputs("\n", tmpfile);

  if ( record )
    saverec = strcpy(saverec, record);

/* now we write the eqparms section into the tmpfile */

  zygote.PrintParameters(tmpfile, p, title, ndigits);

  fprintf(tmpfile, "\n");

/* ... and then write all the rest */

  if ( record )
    fputs(saverec, tmpfile);

  while ( (record=fgets(record, MAX_RECORD, outfile)) )
    fputs(record, tmpfile);

  fclose(outfile);
  fclose(tmpfile);

/* rename tmpfile into new file */
#ifdef BG
  if ( -1 == rename(temp, filename) )
    error("WriteParameters: error renaming temp file %s");
#else

  sprintf(shell_cmd, "cp -f %s %s", temp, filename);

  if ( -1 == system(shell_cmd) )
    error("WriteParameters: error renaming temp file %s");

  if ( remove(temp) )
    warning("WriteParameters: temp file %s could not be deleted", temp);
#endif

/* clean up */

  free(temp);
  free(record_ptr);
  free(saverec);
  free(shell_cmd);
}

void fly::serialize(void *buf) const
{
    double *dest = static_cast<double *>(buf); // new style cast
    for (int i = 0; i < nparams; ++i) {
        dest[i] = *(ptab[i].param);
    }
}

void fly::deserialize(void const *buf)
{
    double const *from = static_cast<double const *>(buf);
    for (int i = 0; i < nparams; ++i) {
        *(ptab[i].param) = from[i];
    }
    score_valid = false;
}

void fly::writeAnswer(const char* title)
{
    if (outfile.empty()){
        outfile = infile;
    } else {
        char * shell_cmd = (char *)calloc(MAX_RECORD, sizeof(char));
        sprintf(shell_cmd, "cp -f %s %s", infile.c_str(), outfile.c_str());
        if ( -1 == system(shell_cmd) )
            error("FinalMove: error creating output file %s", outfile.c_str());
        free(shell_cmd);
    }
    EqParms *parm = zygote.GetParameters();
    WriteParameters(outfile.c_str(), parm, title, ndigits);

}

void fly::state2theta(const void* buf, double* theta)
{
    double const *from = static_cast<double const *>(buf);
    for (int i = 0; i < nparams; ++i) {
        theta[i] = abs(from[i] - *(ptab[i].param));
    }
}

void fly::scramble()
{
    SearchSpace *limits = score.GetLimits();
    bool penaltyflag = false;
    if (limits->pen_vec != NULL) {
        limits = score.Penalty2Limits();
        penaltyflag = true;
    }
    EqParms *param   = zygote.GetParameters();
    double out, T_pen, m_pen;
    srand48(getpid());
    int i,j;
    for ( i=0; i<defs.ngenes; i++ ) {

        if ( tweak.Rtweak[i] == 1 ) {                      /* scramble Rs here */
            out = drand48();
            param->R[i] = limits->Rlim[i].lower +
                    out * (limits->Rlim[i].upper - limits->Rlim[i].lower);
        }

        for ( j=0; j<defs.ngenes; j++ ) {                  /* scramble Ts here */
            if ( tweak.Ttweak[(i*defs.ngenes)+j] == 1 ) {
                out = drand48();
                if ( !penaltyflag ) {
                    param->T[(i*defs.ngenes)+j] =
                            limits->Tlim[(i*defs.ngenes)+j].lower +
                            out * (limits->Tlim[(i*defs.ngenes)+j].upper -
                                    limits->Tlim[(i*defs.ngenes)+j].lower);
                } else {
                    T_pen = limits->Tlim[(i*defs.ngenes)+j].lower +
                           out * (limits->Tlim[(i*defs.ngenes)+j].upper -
                                    limits->Tlim[(i*defs.ngenes)+j].lower);
                    param->T[(i*defs.ngenes)+j] = T_pen / limits->pen_vec[j+2];
                }                                      /* above is T_pen / vmax[i] */
            }
        }

        for ( j=0; j<defs.egenes; j++ ) {                  /* scramble Es here */
            if ( tweak.Etweak[(i*defs.egenes)+j] == 1 ) {
                out = drand48();
                if ( !penaltyflag ) {
                    param->E[(i*defs.egenes)+j] =
                            limits->Elim[(i*defs.egenes)+j].lower +
                            out * (limits->Elim[(i*defs.egenes)+j].upper -
                                    limits->Elim[(i*defs.egenes)+j].lower);
                } else {
                    T_pen = limits->Elim[(i*defs.egenes)+j].lower +
                            out * (limits->Elim[(i*defs.egenes)+j].upper -
                                    limits->Elim[(i*defs.egenes)+j].lower);
                    param->E[(i*defs.egenes)+j] = T_pen /
                            limits->pen_vec[defs.ngenes+j+2];
                }                                      /* above is T_pen / vmax[i] */
            }
        }

        if ( tweak.mtweak[i] == 1 ) {                      /* scramble ms here */
            out = drand48();
            if ( !penaltyflag ) {
                param->m[i] = limits->mlim[i].lower +
                        out * (limits->mlim[i].upper - limits->mlim[i].lower);
            } else {
                m_pen = limits->mlim[i].lower +
                        out * (limits->mlim[i].upper - limits->mlim[i].lower);
                param->m[i] = m_pen / limits->pen_vec[1];           /* m_pen / mmax */
            }
        }

        if ( tweak.htweak[i] == 1 ) {                      /* scramble hs here */
            out = drand48();
            param->h[i] = limits->hlim[i].lower +
                    out * (limits->hlim[i].upper - limits->hlim[i].lower);
        }

        if ( tweak.lambdatweak[i] == 1 ) {            /* scramble lambdas here */
            out = drand48();
            param->lambda[i] = limits->lambdalim[i].lower +
                    out * (limits->lambdalim[i].upper - limits->lambdalim[i].lower);
        }

        if ( tweak.tautweak[i] == 1 ) {                      /* scramble
      delays here */
            out = drand48();
            param->tau[i] = limits->taulim[i].lower +
                    out * (limits->taulim[i].upper - limits->taulim[i].lower);
        }
    }

    /* ds need to be srambled separately, since for diffusion schedules A & C  */
    /* there's just one single d                                               */

    if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') ) {
        if ( tweak.dtweak[0] == 1 ) {
            out        = drand48();
            param->d[0] = limits->dlim[0].lower +
                    out * (limits->dlim[0].upper - limits->dlim[0].lower);
        }
    } else {
        for ( i=0; i<defs.ngenes; i++ ) {
            if ( tweak.dtweak[i] == 1 ) {
                out        = drand48();
                param->d[i] = limits->dlim[i].lower +
                        out * (limits->dlim[i].upper - limits->dlim[i].lower);
            }
        }
    }
}
