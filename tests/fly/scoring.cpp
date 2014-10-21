/*
 * scoring.cpp
 *
 *  Created on: Feb 8, 2013
 *      Author: zhlou
 */
#include <cstring>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <cassert>
#include "maternal.h"
#include "zygotic.h"
#include "scoring.h"
#include "solvers.h"

using namespace std;

scoring::scoring(FILE *fp, zygotic &zy, int flags, int ndigits, double step, double acc, FILE *slog,
                 const char *infile, int in_debug) :
                 Zygote(zy), TheMaternal(Zygote.get_Maternal()),
                 defs(TheMaternal.getProblem()), stepsize(step), accuracy(acc),
                 slogptr(slog), filename(infile), debug(in_debug)
{

    ndatapoints = 0;
    delay_solver = Zygote.get_delay_solver();
    gutparms.flag = flags;
    gutparms.ndigits = ndigits;

    // InitFacts
    int               i, j;                               /* local loop counter */

    Dlist             *d_inlist;            /* temporary linked list for facts */

    Slist             *s_current;                      /* types from data file */

    nalleles = TheMaternal.GetNAlleles();
    if ( !(facttype=(GenoType *)calloc(nalleles,
                                       sizeof(GenoType))) )
      error("InitFacts: could not allocate facttype struct");

  /*** for loop: read the data for each genotype *****************************/

    for(s_current=TheMaternal.GetGenotypes(), i=0; s_current; s_current=s_current->next, i++) {

        if( !(d_inlist=ReadData(fp, s_current->fact_section, &ndatapoints,
                              defs.ngenes)))
            error("InitFacts: no Dlist to initialize facts");
        else {
            if (!(facttype[i].genotype=(char *)calloc(MAX_RECORD, sizeof(char))))
                error("InitFacts: could not allocate facts genotype string");
            facttype[i].genotype = strcpy(facttype[i].genotype, s_current->genotype);
            facttype[i].ptr.facts = List2Facts(d_inlist);

            free_Dlist(d_inlist);
        }
    }
    // initTTs here

    /* tt.ptr.times is a DArrPtr static to score.c.                            */
    /* Each DArrPtr in the array points to an array of times for which we have */
    /* data for each genotype (if you know what I mean...)                     */

    if (!(tt = (GenoType *) calloc(nalleles, sizeof(GenoType))))
        error("InitTTs: could not allocate tt struct");

    /* the following loop copies the times from the Facts section of GenoTab   */
    /* to the tt array of DArrPtrs                                             */

    for (i = 0; i < nalleles; i++) {
        if (!(tt[i].genotype = (char *) calloc(MAX_RECORD, sizeof(char))))
            error("InitTTs: could not allocate tt genotype string");
        tt[i].genotype = strcpy(tt[i].genotype, facttype[i].genotype);
        tt[i].ptr.times.size = 1 + facttype[i].ptr.facts->size;
        if (!(tt[i].ptr.times.array = (double *) calloc(tt[i].ptr.times.size,
                                                        sizeof(double))))
            error("InitTTs: could not allocate bt array");
        ;
        tt[i].ptr.times.array[0] = 0.;

        for (j = 1; j < tt[i].ptr.times.size; j++)
            tt[i].ptr.times.array[j] =
                    facttype[i].ptr.facts->record[j - 1].time;
    }

    // InitLimits
    limits = ReadLimits(fp);

    if (limits->pen_vec != NULL) { // InitPenalty: installs mmax and vmax into
                                   // penalty vector
        char *g_type; /* genotype string of penalty data to be read */

        char *factpen_section; /* title of penalty data section */
        char *extpen_section; /* title of penalty data section for external inputs*/
        char *bcdpen_section; /* title of maternal penalty data section */

        Dlist *d_current;

        Blist *b_inlist; /* linked list for reading maternal */
        Blist *b_current; /* penalty data */

        double *mmax; /* max. bcd protein conc. in maternal penalty data */
        double *vmax; /* array for max. prot concs. in penalty data */

        int ndp = 0; /* dummy for ReadData, no need to count datapts here */

        factpen_section = (char *) calloc(MAX_RECORD, sizeof(char));
        bcdpen_section = (char *) calloc(MAX_RECORD, sizeof(char));
        extpen_section = (char *) calloc(MAX_RECORD, sizeof(char));
        g_type = (char *) calloc(MAX_RECORD, sizeof(char));

        mmax = (limits->pen_vec) + 1; /* mmax stored in limits->pen_vec[1] */
        vmax = (limits->pen_vec) + 2; /* vmax stored in penalty vector array */

        *mmax = -1.;
        for (i = 0; i < defs.ngenes + defs.egenes; i++)
            vmax[i] = -1.;

        /* loop to read penalty data for all genotypes */

        for (s_current = TheMaternal.GetGenotypes(); s_current; s_current =
                s_current->next) {

            g_type = strcpy(g_type, s_current->genotype);

            /* this part reads the facts data */

            if (!(d_inlist = ReadData(fp, s_current->fact_section, &ndp,
                                      defs.ngenes)))
                error("InitPenalty: no Dlist to initialize penalty");
            else {
                for (d_current = d_inlist; d_current;
                        d_current = d_current->next)
                    for (j = 0; j < defs.ngenes; j++)
                        if (d_current->d[j + 1] > vmax[j])
                            vmax[j] = d_current->d[j + 1];

                free_Dlist(d_inlist);
            }

            /* read additional penalty sections if necessary */

            sprintf(factpen_section, "penalty_data.%s", g_type);
            if (FindSection(fp, factpen_section)) {
                if ((d_inlist = ReadData(fp, factpen_section, &ndp, defs.ngenes))) {
                    for (d_current = d_inlist; d_current;
                            d_current = d_current->next)
                        for (j = 0; j < defs.ngenes; j++)
                            if (d_current->d[j + 1] > vmax[j])
                                vmax[j] = d_current->d[j + 1];
                } else {
                    error("InitPenalty: error reading penalty section for genotype %s",
                          g_type);
                }

                free_Dlist(d_inlist);
            }

            /* just in case if all values are -1 */

            for (i = 0; i < defs.ngenes; i++)
                if (vmax[i] == -1.)
                    vmax[i] = 255.; /* set max value to 255 */

            /* this part reads the external inputs data */

            if (defs.egenes > 0) {

                if (!(d_inlist = ReadInterpData(fp, s_current->ext_section,
                                                defs.egenes, &ndp)))
                    error("InitPenalty: no Dlist to initialize penalty");
                else {
                    for (d_current = d_inlist; d_current;
                            d_current = d_current->next)
                        for (j = 0; j < defs.egenes; j++)
                            if (d_current->d[j + 1] > vmax[j + defs.ngenes])
                                vmax[j + defs.ngenes] = d_current->d[j + 1];

                    free_Dlist(d_inlist);
                }

                /* read additional penalty section if necessary */

                sprintf(extpen_section, "external_penalty_data.%s", g_type);
                if (FindSection(fp, extpen_section)) {
                    if ((d_inlist = ReadInterpData(fp, extpen_section,
                                                   defs.egenes, &ndp))) {
                        for (d_current = d_inlist; d_current; d_current =
                                d_current->next)
                            for (j = 0; j < defs.egenes; j++)
                                if (d_current->d[j + 1] > vmax[j + defs.ngenes])
                                    vmax[j + defs.ngenes] = d_current->d[j + 1];
                    } else {
                        error("InitPenalty: error reading penalty section for genotype %s",
                              g_type);
                    }

                    free_Dlist(d_inlist);
                }

                /* just in case if all values are -1 */

                for (i = 0; i < defs.egenes; i++)
                    if (vmax[i + defs.ngenes] == -1.)
                        vmax[i] = 255.; /* set max value to 255 */

            }

            /* this part reads the bicoid data */

            if (!(b_inlist = TheMaternal.ReadBicoid(fp, s_current->bcd_section)))
                error("InitPenalty: no Blist to initialize penalty");
            else {
                for (b_current = b_inlist; b_current;
                        b_current = b_current->next)
                    if (b_current->conc > *mmax)
                        *mmax = b_current->conc;

                free_Blist(b_inlist);
            }

            /* this part read maternal penalty data */

            sprintf(bcdpen_section, "maternal_penalty_data.%s", g_type);
            if (FindSection(fp, bcdpen_section)) {
                if ((b_inlist = TheMaternal.ReadBicoid(fp, bcdpen_section))) {
                    for (b_current = b_inlist; b_current;
                            b_current = b_current->next)
                        if (b_current->conc > *mmax)
                            *mmax = b_current->conc;
                } else {
                    error("InitPenalty: error reading mat. penalty sect. for genot. %s",
                          g_type);
                }

                free_Blist(b_inlist);
            }

            /* just in case there is no bicoid gradient in the data file */

            if (*mmax == -1.)
                *mmax = 255.; /* set max value to 255 */

        }

        free(factpen_section);
        free(bcdpen_section);
        free(extpen_section);
        free(g_type);

    } // end of initPenalty
    if (delay_solver)
        InitHistory(fp);
    else {
        polations = NULL;
        extinp_polations = NULL;
    }





}



scoring::~scoring()
{
    // TODO free facttype?
    // TODO free tt?
    // TODO free polations?
    // TODO free extinp_polations?
    // TODO files?
    delete limits;
}


/*** List2Facts: takes a Dlist and returns the corresponding DataTable *****
 *               structure we use for facts data.                          *
 *                                                                         *
 * An extensive comment about indices: (by JJ) *****************************
 *                                                                         *
 *               All data points with -1 as a value are NOT read from the  *
 *               data file. The difference between such ignored and zero   *
 *               values is crucial: -1 data points WILL NOT BE COMPARED TO *
 *               simulation data, whereas 0 means NO PROTEIN AT THAT TIME  *
 *               IN THAT NUCLEUS.                                          *
 *               Index numbers help maintain the integrity of the data.    *
 *               An index number is defined as the array index at which a  *
 *               protein concentration would be if the data was complete,  *
 *               i.e. available for all nuclei at all times. In this way   *
 *               a sparse set of data can be compared to a complete set of *
 *               simulation output.                                        *
 *               Thus, indices are defined as starting from 1 for each     *
 *               DataRecord (each time step) and increase by one for each  *
 *               gene in each nucleus in the order predefined by JR.       *
 *                                                                         *
 ***************************************************************************/

DataTable *scoring::List2Facts(Dlist *inlist) {
    int i = 0;
    int j; /* local loop counters */

    double now = -999999999.; /* assigns data to specific time */

    Dlist *current; /* holds current element of Dlist */

    DataTable *D; /* local copy of DataTable */

    D = (DataTable *) malloc(sizeof(DataTable));
    /* Initialize DataTable structure */
    D->size = 0;
    D->record = NULL;

    /*** for loop: steps through linked list and transfers facts into Data-    *
     *             Records, one for each time step                             *
     ***************************************************************************/

    for (current = inlist; current; current = current->next) {

        if (current->d[0] != now) { /* a new time point: allocate */
            now = current->d[0]; /* the time is now! */
            D->size++; /* one DataRecord for each time */
            D->record = /* allocate DataRecord */
            (DataRecord *) realloc(D->record, D->size * sizeof(DataRecord));

            D->record[D->size - 1].time = now; /* next three lines define */
            D->record[D->size - 1].size = 0; /* DataRecord for each */
            D->record[D->size - 1].array = NULL; /* time step */
            i = 0;
        }

        for (j = 1; j <= Zygote.get_ngenes(); j++) { /* always: read concs into array */
            if (current->d[j] != IGNORE) { /* valid conc? if IGNORE -> ignore! */
                D->record[D->size - 1].size++; /* one more in this record */
                D->record[D->size - 1].array = /* reallocate memory for array! */
                (DataPoint *)realloc(D->record[D->size - 1].array,
                        D->record[D->size - 1].size * sizeof(DataPoint));

                /* the following two lines assign concentration value and index ************/

                D->record[D->size - 1].array[D->record[D->size - 1].size - 1].conc =
                        current->d[j];
                D->record[D->size - 1].array[D->record[D->size - 1].size - 1].index =
                        i;
            }

            i++;

        }
    }

    return D;
}

/*** ReadLimits: reads the limits section of a data file and returns the  **
 *               approriate SearchSpace struct to the calling function     *
 ***************************************************************************/

SearchSpace *scoring::ReadLimits(FILE *fp)
{
    SearchSpace       *l_limits;

    char              *record;    /* string for reading whole line of limits */

    int               i, j;                            /* local loop counter */

    char              **fmt1;     /* format strings for reading lower limits */
    char              **fmt2;     /* format strings for reading upper limits */
    char              **efmt1;     /* format strings for reading lower
  limits of external inputs*/
    char              **efmt2;     /* format strings for reading upper
  limits of external inputs*/
    char              *skip;               /* string of values to be skipped */

    const char        read_fmt[]  = "%lg";                  /* read a double */
    const char        skip_fmt1[] = "%*lg, %*lg) (";   /* ignore lower limit */
    const char        skip_fmt2[] = "%*lg) (%*lg, ";   /* ignore upper limit */

    l_limits = new SearchSpace;

    record = (char *)calloc(MAX_RECORD, sizeof(char *));

    skip = (char *)calloc(MAX_RECORD, sizeof(char *));
    fmt1 = (char **)calloc(defs.ngenes, sizeof(char *));
    fmt2 = (char **)calloc(defs.ngenes, sizeof(char *));

    if (defs.egenes > 0) {

        efmt1 = (char **)calloc(defs.egenes, sizeof(char *));
        efmt2 = (char **)calloc(defs.egenes, sizeof(char *));
    }

    /* create format strings according to the number of genes */

    skip = strcpy(skip, "(");                      /* first for lower limits */

    for ( i=0; i<defs.ngenes; i++ ) {
        fmt1[i] = (char *)calloc(MAX_RECORD, sizeof(char));
        fmt1[i] = strcpy(fmt1[i], skip);
        fmt1[i] = strcat(fmt1[i], read_fmt);
        skip    = strcat(skip, skip_fmt1);
    }

    skip = strcpy(skip, "(%*lg, ");                 /* then for upper limits */

    for ( i=0; i<defs.ngenes; i++ ) {
        fmt2[i] = (char *)calloc(MAX_RECORD, sizeof(char));
        fmt2[i] = strcpy(fmt2[i], skip);
        fmt2[i] = strcat(fmt2[i], read_fmt);
        skip    = strcat(skip, skip_fmt2);
    }

    /* create format strings according to the number of external inputs  -
     * this could have been done with just fmt1 and fmt2, picking their
     * size to be whatever was bigger - ngenes or egenes, however this
     * maintains the seperation of reading in normal things and external
     * inputs and avoid confusion */

    if (defs.egenes > 0) {

        skip = strcpy(skip, "(");                      /* first for lower limits */

        for ( i=0; i<defs.egenes; i++ ) {
            efmt1[i] = (char *)calloc(MAX_RECORD, sizeof(char));
            efmt1[i] = strcpy(efmt1[i], skip);
            efmt1[i] = strcat(efmt1[i], read_fmt);
            skip    = strcat(skip, skip_fmt1);
        }

        skip = strcpy(skip, "(%*lg, ");                 /* then for upper limits */

        for ( i=0; i<defs.egenes; i++ ) {
            efmt2[i] = (char *)calloc(MAX_RECORD, sizeof(char));
            efmt2[i] = strcpy(efmt2[i], skip);
            efmt2[i] = strcat(efmt2[i], read_fmt);
            skip    = strcat(skip, skip_fmt2);
        }

    }

    /* find limits section and check if penalty or explicit ranges are used    */

    fp = FindSection(fp, "limits");                   /* find limits section */
    if( !fp )
        error("ReadLimits: cannot locate limits section");

    fscanf(fp, "%*s\n");                    /* advance past first title line */

    record = fgets(record, MAX_RECORD, fp);    /* read Lamda for penalty */

    if ( !strncmp(record, "N/A", 3) ) {
        l_limits->pen_vec = NULL;
    } else {
        l_limits->pen_vec = new double[2 + defs.ngenes + defs.egenes];
        if ( 1 != sscanf(record, "%lg", l_limits->pen_vec) )
            error("ReadLimits: error reading Lambda for penalty");
    }


    /* initialize limits struct according to penalty or not   *
     * first the stuff that's not dependent on penalty        */

    l_limits->Rlim = new Range[defs.ngenes];
    if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') ) {
        l_limits->dlim = new Range[1];
    } else {
        l_limits->dlim = new Range[defs.ngenes];
    }
    l_limits->lambdalim = new Range[defs.ngenes];
    l_limits->taulim = new Range[defs.ngenes];

    /* ... and now the stuff that depends on whether we use penalty or not     */

    if ( l_limits->pen_vec == NULL ) {
        l_limits->Tlim = new Range[defs.ngenes * defs.ngenes];
        if (defs.egenes > 0) {
            l_limits->Elim = new Range[defs.ngenes * defs.egenes];
        } else {
            l_limits->Elim = NULL;
        }
        l_limits->mlim = new Range[defs.ngenes];
        l_limits->hlim = new Range[defs.ngenes];
    } else {
        l_limits->Tlim = NULL;
        l_limits->Elim = NULL;
        l_limits->mlim = NULL;
        l_limits->hlim = NULL;
    }

    /* ... and finally read the actual values from the data file               */

    fscanf(fp, "%*s\n");                          /* advance past title line */

    record = fgets(record, MAX_RECORD, fp);     /* loop to read R limits */
    for ( i=0; i<defs.ngenes; i++ ) {
        if ( 1 != sscanf(record, (const char *)fmt1[i],
                &l_limits->Rlim[i].lower) )
            error("ReadLimits: error reading promoter strength (R) limits");
        if ( 1 != sscanf(record, (const char *)fmt2[i],
                &l_limits->Rlim[i].upper) )
            error("ReadLimits: error reading promoter strength (R) limits");
    }

    fscanf(fp, "%*s\n");

    /* using explicit limits? */

    if ( l_limits->pen_vec == 0 ) {
        for ( i=0; i<defs.ngenes; i++ ) {     /* loops to read T matrix limits */
            record = fgets(record, MAX_RECORD, fp);
            for ( j=0; j<defs.ngenes; j++ ) {
                if ( 1 != sscanf(record, (const char *)fmt1[j],
                        &l_limits->Tlim[(i * defs.ngenes) + j].lower ) )
                    error("ReadLimits:: error reading T matrix limits");
                if ( 1 != sscanf(record, (const char *)fmt2[j],
                        &l_limits->Tlim[(i * defs.ngenes) + j].upper ) )
                    error("ReadLimits:: error reading T matrix limits");
            }
        }

        fscanf(fp, "%*s\n");

        if (defs.egenes > 0) {

            for ( i=0; i<defs.ngenes; i++ ) {     /* loops to read E matrix limits */
                record = fgets(record, MAX_RECORD, fp);
                for ( j=0; j<defs.egenes; j++ ) {
                    if ( 1 != sscanf(record, (const char *)efmt1[j],
                            &l_limits->Elim[(i * defs.egenes) + j].lower ) )
                        error("ReadLimits:: error reading E matrix limits");
                    if ( 1 != sscanf(record, (const char *)efmt2[j],
                            &l_limits->Elim[(i * defs.egenes) + j].upper ) )
                        error("ReadLimits:: error reading E matrix limits");
                }
            }
        } else {

            fscanf(fp, "%*s\n");

        }

        fscanf(fp, "%*s\n");

        record = fgets(record, MAX_RECORD, fp);   /* loop to read m limits */
        for ( i=0; i<defs.ngenes; i++ ) {
            if ( 1 != sscanf(record, (const char *)fmt1[i],
                    &l_limits->mlim[i].lower) )
                error("ReadLimits: error reading maternal interconnect (m) limits");
            if ( 1 != sscanf(record, (const char *)fmt2[i],
                    &l_limits->mlim[i].upper) )
                error("ReadLimits: error reading maternal interconnect (m) limits");
        }

        fscanf(fp, "%*s\n");

        record = fgets(record, MAX_RECORD, fp);   /* loop to read h limits */
        for ( i=0; i<defs.ngenes; i++ ) {
            if ( 1 != sscanf(record, (const char *)fmt1[i],
                    &l_limits->hlim[i].lower) )
                error("ReadLimits: error reading promoter threshold (h) limits");
            if ( 1 != sscanf(record, (const char *)fmt2[i],
                    &l_limits->hlim[i].upper) )
                error("ReadLimits: error reading promoter threshold (h) limits");
        }

        /* using penalty? -> ignore this part of the limit section */

    } else {

        for ( i=0; i<7; i++ )
            fscanf(fp, "%*s\n");

    }

    fscanf(fp, "%*s\n");        /* diffusion paramter limit(s) are read here */

    record = fgets(record, MAX_RECORD, fp);
    if ( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C') ) {
        if ( 1 != sscanf(record, (const char *)fmt1[0],
                &l_limits->dlim[0].lower) )
            error("ReadLimits: error reading diffusion parameter limit (d)");
        if ( 1 != sscanf(record, (const char *)fmt2[0],
                &l_limits->dlim[0].upper) )
            error("ReadLimits: error reading diffusion parameter limit (d)");
    } else {
        for ( i=0; i<defs.ngenes; i++ ) {
            if ( 1 != sscanf(record, (const char *)fmt1[i],
                    &l_limits->dlim[i].lower) )
                error("ReadLimits: error reading diffusion parameter limits (d)");
            if ( 1 != sscanf(record, (const char *)fmt2[i],
                    &l_limits->dlim[i].upper) )
                error("ReadLimits: error reading diffusion parameter limits (d)");
        }
    }

    fscanf(fp, "%*s\n");

    record = fgets(record, MAX_RECORD, fp);  /* loop to read lambda lims */
    for ( i=0; i<defs.ngenes; i++ ) {
        if ( 1 != sscanf(record, (const char *)fmt1[i],
                &l_limits->lambdalim[i].upper) )
            error("ReadLimits: error reading lambda limits");
        l_limits->lambdalim[i].upper = log(2.) / l_limits->lambdalim[i].upper;
        if ( 1 != sscanf(record, (const char *)fmt2[i],
                &l_limits->lambdalim[i].lower) )
            error("ReadLimits: error reading lambda limits");
        l_limits->lambdalim[i].lower = log(2.) / l_limits->lambdalim[i].lower;
    }

    fscanf(fp, "%*s\n");

    record = fgets(record, MAX_RECORD, fp);  /* loop to read tau lims */
    for ( i=0; i<defs.ngenes; i++ ) {
        if ( 1 != sscanf(record, (const char *)fmt1[i],
                &l_limits->taulim[i].lower) )
            error("ReadLimits: error reading tau limits");
        if ( 1 != sscanf(record, (const char *)fmt2[i],
                &l_limits->taulim[i].upper) )
            error("ReadLimits: error reading tau limits");
    }

    free(record);
    free(skip);

    for (i=0; i<defs.ngenes; i++ ) {
        free(fmt1[i]);
        free(fmt2[i]);
    }

    for (i=0; i<defs.egenes; i++ ) {
        free(efmt1[i]);
        free(efmt2[i]);
    }

    free(fmt1);
    free(fmt2);

    if (defs.egenes > 0) {
        free(efmt1);
        free(efmt2);
    }

    return l_limits;
}

void scoring::InitHistory(FILE *fp)
{

    int               i,ii;
    DataTable         **temp_table;
    double            *theta_discons;
    int               theta_discons_size;
    double            *temp_divtable;
    double            *temp_durations;
    Slist             *curr;     /* We will read in the
                               genotypes to set
                               the interp_dat,
                               bias_dat tables */
    int nalleles = TheMaternal.GetNAlleles();



    if ( !(temp_table=(DataTable **)calloc(1, sizeof(DataTable *))) )
        error("InitHistory: could not allocate temp_table struct");


    if ( !(polations=(InterpObject *)calloc(nalleles,
                    sizeof(InterpObject))) )
        error("InitHistory: could not allocate interp_tab struct");

    for(curr=TheMaternal.GetGenotypes(), i=0; curr; curr=curr->next, i++) {

        GetInterp(fp, curr->hist_section, defs.ngenes, temp_table);

        DoInterp(*temp_table, polations+i, defs.ngenes);
        FreeFacts(*temp_table);

        theta_discons = TheMaternal.Get_Theta_Discons(&theta_discons_size);


        (polations+i)->fact_discons = (double *) realloc((polations+i)->fact_discons,
                ((polations+i)->fact_discons_size+theta_discons_size)*sizeof(double));

        for (ii=0; ii < theta_discons_size; ii++)
            (polations+i)->fact_discons[(polations+i)->fact_discons_size + ii] = theta_discons[ii];

        (polations+i)->fact_discons_size += theta_discons_size;
        free(theta_discons);

        if ( defs.ndivs > 0 ) {
            if ( !(temp_divtable = TheMaternal.GetDivtable()) )
                error("InitHistory: error getting temp_divtable");
            if ( !(temp_durations = TheMaternal.GetDurations()) )
                error("Inithistory: error getting division temp_durations");

            for (ii=0; ii<defs.ndivs; ii++) {
                (polations+i)->fact_discons = (double *) realloc((polations+i)->fact_discons,
                        ((polations+i)->fact_discons_size+4)*sizeof(double));

                (polations+i)->fact_discons[(polations+i)->fact_discons_size] =
                    temp_divtable[ii];
                (polations+i)->fact_discons[(polations+i)->fact_discons_size+1] =
                    temp_divtable[ii] + EPSILON;
                (polations+i)->fact_discons[(polations+i)->fact_discons_size+2] =
                    temp_divtable[ii]-temp_durations[ii];
                (polations+i)->fact_discons[(polations+i)->fact_discons_size+3] =
                    temp_divtable[ii]-temp_durations[ii] + EPSILON;

                (polations+i)->fact_discons_size += 4;
            }
        }


        qsort((void *) (polations+i)->fact_discons,
                (polations+i)->fact_discons_size,
                sizeof(double),(int (*) (const void*,const void*)) compare);

    }

    /* Initializing the full set of nuclei based on the lineages of the history,
      please make sure that all of the alleles' lineages are the same */



    free(temp_table);

    return;

}
void scoring::GetInterp(FILE *fp, char *title, int num_genes,
                                        DataTable **interp_tables)
{

    int ndatapoints=0;
    Dlist *list_dat;

    if(!(list_dat=ReadInterpData(fp, title, num_genes,
                                                &ndatapoints)))
        error("SetFacts: no facts section for history for"
                                                    "delay solver");

    *interp_tables = TheMaternal.List2Interp(list_dat, num_genes);

    free_Dlist(list_dat);

    return;
}

void scoring::DoInterp(DataTable *interp_dat, InterpObject *interp_res, int
num_genes)
{

    int i,j,kk;
    NArrPtr Nptrfacts;
    int maxind;
    double *x, *y, t;
    int accepted,currsize;
    double *blug;
    gsl_spline *temp_spline;
    gsl_interp_accel *temp_acc;

    Nptrfacts = Dat2NArrPtr(interp_dat,&maxind);

    interp_res->maxsize = Nptrfacts.array[maxind].state.size;
    interp_res->maxtime = Nptrfacts.array[maxind].time;

/*  printf("maxsize:%d, maxtime:%f\n",maxsize,maxtime);*/

/*  PrintBlastoderm(stdout, Nptrfacts, "Wloo", 2, num_genes);*/

    interp_res->func.size = Nptrfacts.size;
    interp_res->func.array = (NucState *) calloc(interp_res->func.size,
                                                sizeof(NucState));
    interp_res->slope.size = Nptrfacts.size;
    interp_res->slope.array = (NucState *) calloc(interp_res->slope.size,
                                                sizeof(NucState));

    interp_res->fact_discons_size = Nptrfacts.size;
    interp_res->fact_discons = (double *) calloc(interp_res->fact_discons_size,
                                                sizeof(double));
    for (i=0; i < interp_res->func.size; i++)
    {

        interp_res->func.array[i].time=Nptrfacts.array[maxind].time;
        interp_res->slope.array[i].time=Nptrfacts.array[i].time;

        interp_res->fact_discons[i] = Nptrfacts.array[i].time;

        interp_res->func.array[i].state.size=
                                Nptrfacts.array[maxind].state.size;
        interp_res->slope.array[i].state.size=
                                Nptrfacts.array[maxind].state.size;

        interp_res->func.array[i].state.array=(double *)
                    calloc(interp_res->func.array[i].state.size,
                                                    sizeof(double));
        interp_res->slope.array[i].state.array=(double *)
                    calloc(interp_res->slope.array[i].state.size,
                                                    sizeof(double));

/*      printf("Going from: %d to %d\n",Nptrfacts.array[i].state.size,
                                interp_res->func.array[i].state.size);
*/
        if (interp_res->func.array[i].state.size >=
                                        Nptrfacts.array[i].state.size)
            delay_solver->Go_Forward(interp_res->func.array[i].state.array,
                Nptrfacts.array[i].state.array,
                TheMaternal.GetStartLinIndex(interp_res->func.array[i].time),
                TheMaternal.GetStartLinIndex(Nptrfacts.array[i].time), num_genes);
        else
            delay_solver->Go_Backward(interp_res->func.array[i].state.array,
                Nptrfacts.array[i].state.array,
                TheMaternal.GetStartLinIndex(interp_res->func.array[i].time),
                TheMaternal.GetStartLinIndex(Nptrfacts.array[i].time), num_genes);

/*      interp_res->func.array[i].time = Nptrfacts.array[i].time;*/

    }



    for (i=0; i < interp_res->maxsize; i++)
    {
        /* If the last data point is zero, make it equal to the last
         * non-zero value, so that the concentration at a time
         * after the last non -1 is maintained at that value */

        if (interp_res->func.array[Nptrfacts.size-1].state.array[i] == -1.)
        {

            kk = Nptrfacts.size-2;
            while ((interp_res->func.array[kk].state.array[i] == -1.) &&
                                        (kk >= 0)) kk--;
            if (kk == -1)
            {

                printf("All data is -1!\n");
                exit(1);

            }

            interp_res->func.array[Nptrfacts.size-1].state.array[i] =
                interp_res->func.array[kk].state.array[i];

        }

        /* If the first value in a data set is -1, set it to 0. */

        if (interp_res->func.array[0].state.array[i] == -1.)
            interp_res->func.array[0].state.array[i] = 0.;

        /* Now weed-out the -1s so you can interpolate linearly
         * between time points for which you have data */

        x = (double *) calloc(Nptrfacts.size, sizeof(double));
        y = (double *) calloc(Nptrfacts.size, sizeof(double));

        accepted = 0;
        currsize = Nptrfacts.size;


        for (j=0; j < Nptrfacts.size; j++)
        {

            if (interp_res->func.array[j].state.array[i] != -1.)
            {

                x[accepted] = Nptrfacts.array[j].time;
                y[accepted] = interp_res->func.array[j].state.array[i];
                accepted++;

            } else {

                currsize--;
                x = (double *)realloc(x, currsize*sizeof(double));
                y = (double *)realloc(y, currsize*sizeof(double));

            }
        }

        temp_acc = gsl_interp_accel_alloc();
        temp_spline = gsl_spline_alloc(gsl_interp_linear, currsize);
        gsl_spline_init(temp_spline,x,y,currsize);
        free(x);
        free(y);

        for (j=0; j < Nptrfacts.size; j++)
        {

            interp_res->func.array[j].state.array[i] =
                gsl_spline_eval(temp_spline, interp_res->slope.array[j].time,
                                                            temp_acc);

        }

        for (j=0; j < Nptrfacts.size; j++)
        {
            if (j < Nptrfacts.size - 1)
                interp_res->slope.array[j].state.array[i] =
                    (interp_res->func.array[j+1].state.array[i]
                    - interp_res->func.array[j].state.array[i])/
                    (interp_res->slope.array[j+1].time
                    - interp_res->slope.array[j].time);
            else interp_res->slope.array[j].state.array[i] = 0.;

        }

        gsl_interp_accel_free(temp_acc);
        gsl_spline_free(temp_spline);
    }


/*  PrintBlastoderm(stdout, interp_res->func, "loo", 2, num_genes);
    PrintBlastoderm(stdout, interp_res->slope, "sloo", 2, num_genes);*/
    FreeSolution(&Nptrfacts);

    return;
}

/*** REAL SCORING CODE HERE ************************************************/


/*** Score: as the name says, score runs the simulation, gets a solution ***
 *          and then compares it to the data using the Eval least squares  *
 *          function                                                       *
 *   NOTE:  both InitZygote and InitScoring have to be called first!       *
 ***************************************************************************/

double scoring::Score(void)
{
    int        i, j, ii;                              /* local loop counters */
    FILE       *fp;                                          /* file pointer */

    char       *debugfile;                             /* name of debug file */

    NArrPtr    answer;                /* stores the Solution from Blastoderm */

    EqParms    *parm;             /* local copy of parameters for comparison */
    /* to limits */

    double     chisq   = 0;                    /* summed squared differences */
    double     penalty = 0;                          /* variable for penalty */


    /* debugging mode: need debugging file name */

    if ( debug )
        debugfile = (char *)calloc(MAX_RECORD, sizeof(char));

    parm = Zygote.GetParameters();                 /* get parameters from zygotic.c */

    /* The following will be called after parms are tweaked, hence it must     *
     * check signs. If it appears cleaner, sign checking could be done by the  *
     * tweaker                                                                 */

    for ( ii=0; ii < defs.ngenes; ii++) { /* NOTE: parm is a global EqParms  */
        /* struct declared in zygotic.h    */
        if (parm->R[ii] < 0)
            parm->R[ii] = -parm->R[ii];       /* Rs and lambdas should ALWAYS be */
        /* greater than zero. So this just */
        if (parm->lambda[ii] < 0)           /* makes sure of that, I guess     */
            parm->lambda[ii] = -parm->lambda[ii];
    }

    /* The following fors and ifs check the searchspace in such a way as to    *
     * return after as few calculations as possible                            */

    for(i=0; i < defs.ngenes; i++) {
        if( parm->R[i] > limits->Rlim[i].upper)
            return(FORBIDDEN_MOVE);
        if( parm->R[i] < limits->Rlim[i].lower)
            return(FORBIDDEN_MOVE);
        if( parm->lambda[i] > limits->lambdalim[i].upper)
            return(FORBIDDEN_MOVE);
        if( parm->lambda[i] < limits->lambdalim[i].lower)
            return(FORBIDDEN_MOVE);
        if( parm->tau[i] > limits->taulim[i].upper)
            return(FORBIDDEN_MOVE);
        if( parm->tau[i] < limits->taulim[i].lower)
            return(FORBIDDEN_MOVE);
    }

    if( (defs.diff_schedule == 'A') || (defs.diff_schedule == 'C' ) ) {
        if( parm->d[0] > limits->dlim[0].upper)
            return(FORBIDDEN_MOVE);
        if( parm->d[0] < limits->dlim[0].lower)
            return(FORBIDDEN_MOVE);
    }
    else {
        for(i=0; i < defs.ngenes; i++){
            if( parm->d[i] > limits->dlim[i].upper)
                return(FORBIDDEN_MOVE);
            if( parm->d[i] < limits->dlim[i].lower)
                return(FORBIDDEN_MOVE);
        }
    }

    /* Penalty stuff below */
    /* If you're using limits on contributors to u, check'em here */

    if( limits->pen_vec == NULL ) {
        for(i=0; i<defs.ngenes; i++) {
            for(j=0; j<defs.ngenes; j++) {
                if (parm->T[(i * defs.ngenes) + j] >
                limits->Tlim[(i * defs.ngenes) + j].upper)
                    return(FORBIDDEN_MOVE);
                if (parm->T[(i * defs.ngenes) + j] <
                        limits->Tlim[(i * defs.ngenes) + j].lower)
                    return(FORBIDDEN_MOVE);
            }

            /* 01/13/10 Manu: This loop is never entered if there are no
             * external inputs */

            for(j=0; j<defs.egenes; j++) {
                if (parm->E[(i * defs.egenes) + j] >
                limits->Elim[(i * defs.egenes) + j].upper)
                    return(FORBIDDEN_MOVE);
                if (parm->E[(i * defs.egenes) + j] <
                        limits->Elim[(i * defs.egenes) + j].lower)
                    return(FORBIDDEN_MOVE);
            }

            if (parm->m[i] > limits->mlim[i].upper)
                return(FORBIDDEN_MOVE);
            if (parm->m[i] < limits->mlim[i].lower)
                return(FORBIDDEN_MOVE);

            if (parm->h[i] > limits->hlim[i].upper)
                return(FORBIDDEN_MOVE);
            if (parm->h[i] < limits->hlim[i].lower)
                return(FORBIDDEN_MOVE);
        }

        /* if you're going to calculate penalty on u, do it here */
        /* following lines calculate exp of sum of squares penalty function */

    } else {
        penalty = GetPenalty();
        if ( penalty == FORBIDDEN_MOVE )
            return FORBIDDEN_MOVE;
        chisq += penalty;
    }

    /* this loop runs the model and sums squared differences for all genotypes */

    InterpObject dummy_polations; // to make Blastoderm happy
    for ( i=0; i<nalleles; i++) {
        if (delay_solver)
            answer = Zygote.Blastoderm(i, facttype[i].genotype, polations[i],
                                       extinp_polations[i],tt[i].ptr.times,
                                       stepsize, accuracy, slogptr);
        else
            answer = Zygote.Blastoderm(i, facttype[i].genotype, dummy_polations,
                                       dummy_polations,tt[i].ptr.times,
                                       stepsize, accuracy, slogptr);
        if ( debug ) {
            sprintf(debugfile, "%s.%s.pout", filename, facttype[i].genotype);
            fp = fopen(debugfile, "w");
            if( !fp ) {
                perror("printscore");
                exit(1);
            }
            Zygote.PrintBlastoderm(fp, answer, (char *)"debug_output", MAX_PRECISION,
                            defs.ngenes);
            fclose(fp);
        }

        if (gutparms.flag)
            chisq += GutEval(answer, i);
        else {
            chisq += Eval(answer, i);
            /*    printf("Genotype %d: %f\n",i,chisq);*/
        }
    }

    if ( debug )
        free( debugfile );

    return chisq;
}

/*** GetPenalty: calculates penalty from static limits, vmax and mmax ******
 *   CAUTION:    InitPenalty must be called first!                         *
 ***************************************************************************/

double scoring::GetPenalty(void)
{
  int           i, j;                                /* local loop counter */
  double        argument = 0;    /* variable for penalty function argument */
  double        penalty  = 0;                      /* holds penalty result */

  EqParms       *parm;                            /* local copy of eqparms */

  double        Lambda;                                  /* penalty Lambda */
  double        mmax;                            /* mmax in penalty vector */
  double        *vmax;          /* pointer to vmax array in penalty vector */

  static int    donethis = 0;         /* KLUDGE: only print this info once */

  if ( limits->pen_vec == NULL )
    return -1;

  Lambda = *((limits->pen_vec));  /* Lambda: first entry in penalty vector */
  mmax   = *((limits->pen_vec) + 1);               /* locate mmax and vmax */
  vmax   =   (limits->pen_vec) + 2;

/* print debugging info */

  if ( debug && !donethis ) {
    printf("Penalty:\n\n");
    printf("Lambda:    %10.8f\n", Lambda);
    printf("mmax:    %6.2f\n", mmax);
    for ( i=0; i<defs.ngenes; i++ )
      printf("vmax[%d]: %6.2f\n", i, vmax[i]);
    for ( i=0; i<defs.egenes; i++ )
      printf("vmax[%d]: %6.2f\n", i+defs.ngenes, vmax[i+defs.ngenes]);
    printf("\n");
  }

/* calculate penalty */

  parm = Zygote.GetParameters();

  for( i=0; i < defs.ngenes; i++) {
    for(j=0; j < defs.ngenes; j++)
      argument += (parm->T[(i*defs.ngenes)+j]*vmax[j])
    * (parm->T[(i*defs.ngenes)+j]*vmax[j]);

    for(j=0; j < defs.egenes; j++)
      argument += (parm->E[(i*defs.egenes)+j]*vmax[j+defs.ngenes])
    * (parm->E[(i*defs.egenes)+j]*vmax[j+defs.ngenes]);

    argument += (parm->m[i] * mmax) * (parm->m[i] * mmax);
    argument += parm->h[i] * parm->h[i];
  }

/* e was 2.718281828 in old code---doubt it matters!! */
/* 88.7228391 is the maximum save value to use with exp (see man exp) */

  if ( (Lambda * argument) > 88.7228391 )
    return FORBIDDEN_MOVE;
  else
    penalty = (exp( Lambda * argument ) - 2.718281828459045);

  if(penalty <= 0)
    penalty = 0;

  donethis = 1;

  return penalty;

}

/*** Eval: scores the summed squared differences between equation solution *
 *         and data. Because the times for states written to the Solution  *
 *         structure are read out of the data file itself, we do not check *
 *         for consistency of times in this function---all times with data *
 *         will be in the table, but the table may also contain additional *
 *         times.                                                          *
 ***************************************************************************/

double scoring::Eval(NArrPtr Solution, int gindex)
{
    DataTable    fact_tab;      /* stores a copy of the Facts (from GenoTab) */
    DataPoint    point;           /* used to extract an element of DataTable */

    int          tindex;                       /* index for facts timepoints */
    int          sindex;                    /* index for Solution timepoints */
    int          vindex;                        /* index for facts datapoint */

    double       time;                      /* time for each facts timepoint */
    double       difference;      /* diff btw data and model (per datapoint) */
    double       *v;                   /* ptr to solution for each timepoint */
    double       chisq = 0;                /* the final score to be returned */

    FILE *score;


    /* extablish a pointer to the appropriate facts section in GenoTab */

    fact_tab = *(facttype[gindex].ptr.facts);
    sindex   = 0;

    /* the following loop adds all squared differences for the whole Solution */

    for (tindex=0; tindex < fact_tab.size; tindex++) {
        time = fact_tab.record[tindex].time;

        /* new time step to be compared? -> get the Solution for this time */

        while( fabs(time - Solution.array[sindex].time) >= BIG_EPSILON )
            sindex++;

        v = Solution.array[sindex].state.array;

        /* this loop steps through the Solution for a specific timepoint and       */
        /* evaluates the squared diffs                                             */

        for (vindex=0; vindex < fact_tab.record[tindex].size; vindex++) {
            point = fact_tab.record[tindex].array[vindex];
            difference = point.conc - v[point.index];
            chisq += difference * difference;
        }
    }

    return chisq;                                              /* that's it! */
}

/*** GutEval: this is the same as Eval, i.e it calculates the summed squa- *
 *            red differences between equation solution and data, with the *
 *            addition that individual squared differences between data-   *
 *            points are written to STDOUT in the unfold output format     *
 ***************************************************************************/

double scoring::GutEval(NArrPtr Solution, int gindex)
{
    int          i, j;                                      /* loop counters */

    DataTable    fact_tab;      /* stores a copy of the Facts (from GenoTab) */
    DataPoint    point;           /* used to extract an element of DataTable */
    NArrPtr      gut;         /* individual square root diff for a datapoint */
    NArrPtr      outgut;                             /* output gut structure */

    int          tindex;                       /* index for facts timepoints */
    int          sindex;                    /* index for Solution timepoints */
    int          vindex;                        /* index for facts datapoint */

    double       time;                      /* time for each facts timepoint */
    double       difference;      /* diff btw data and model (per datapoint) */
    double       *v;                   /* ptr to solution for each timepoint */
    double       chisq = 0;                /* the final score to be returned */
    char         gen_print[MAX_RECORD];   /* for PrintBlastoderm */

    /* extablish a pointer to the appropriate facts section in GenoTab */

    fact_tab = *(facttype[gindex].ptr.facts);
    sindex   = 0;
    gen_print[0] = (char) 48+gindex;

    /* initialize the gut structure */

    gut.array = (NucState *)calloc(Solution.size, sizeof(NucState));
    gut.size  = Solution.size;
    for (i=0; i<gut.size; i++) {
        gut.array[i].time        = Solution.array[i].time;
        gut.array[i].state.size  = Solution.array[i].state.size;
        gut.array[i].state.array =
                (double *)calloc(gut.array[i].state.size, sizeof(double));
        for (j=0; j < gut.array[i].state.size; j++)
            gut.array[i].state.array[j] = -1.0;
    }

    /* the following loop adds all squared differences for the whole Solution */

    for (tindex=0; tindex<fact_tab.size; tindex++) {
        time = fact_tab.record[tindex].time;

        /* new time step to be compared? -> get the Solution for this time */

        while( fabs(time - Solution.array[sindex].time) >= BIG_EPSILON )
            sindex++;

        v = Solution.array[sindex].state.array;

        /* this loop steps through the Solution for a specific time point and      */
        /* evaluates the squared diffs; additionally it stores individual squared  */
        /* diffs in the guts array                                                 */

        for(vindex=0; vindex < fact_tab.record[tindex].size; vindex++) {
            point = fact_tab.record[tindex].array[vindex];
            difference = point.conc - v[point.index];
            chisq += difference * difference;
            gut.array[sindex].state.array[point.index] = difference * difference;
        }
    }

    /* strip gut struct of cell division times and print it to stdout */

    outgut = ConvertAnswer(gut, tt[gindex].ptr.times);
    Zygote.PrintBlastoderm(stdout, outgut, strcat(gen_print, " genotype\n"),
                           gutparms.ndigits, defs.ngenes);

    /* release the guts struct */

    FreeSolution(&gut);
    free(gut.array);

    return chisq;                                              /* that's it! */
}

SearchSpace *scoring::Penalty2Limits()
{
    assert(limits->pen_vec != NULL);
    assert(limits->Tlim == NULL);
    assert(limits->Elim == NULL);
    assert(limits->mlim == NULL);
    assert(limits->hlim == NULL);

    int i,j;
    double Lambda = limits->pen_vec[0];
    Range u, gu, x, y;

    gu.lower = Lambda;
    gu.upper = 1 - Lambda;

    x.lower = ( 2 * gu.lower - 1 );  /* the following calculates the inverse */
    x.upper = ( 2 * gu.upper - 1 );  /* function of g(u) for gu limits above */
                                          /* (see JJs lab notes for details) */
    y.lower = sqrt( 1 - x.lower * x.lower );
    y.upper = sqrt( 1 - x.upper * x.upper );

    u.lower = x.lower / y.lower;
    u.upper = x.upper / y.upper;

    u.lower = u.lower / sqrt(defs.ngenes);  /* this is to compensate for the */
    u.upper = u.upper / sqrt(defs.ngenes);       /* summing up of parameters */



    limits->Tlim = new Range[defs.ngenes * defs.ngenes];


    if (defs.egenes > 0)
        limits->Elim = new Range[defs.ngenes * defs.egenes];

    limits->mlim = new Range[defs.ngenes];
    limits->hlim = new Range[defs.ngenes];

    for ( i=0; i<defs.ngenes; i++ ) {
        for ( j=0; j<defs.ngenes; j++ ) {
            limits->Tlim[(i*defs.ngenes)+j].lower = u.lower;
            limits->Tlim[(i*defs.ngenes)+j].upper = u.upper;
        }

        for ( j=0; j<defs.egenes; j++ ) {
            limits->Elim[(i*defs.egenes)+j].lower = u.lower;
            limits->Elim[(i*defs.egenes)+j].upper = u.upper;
        }


        limits->mlim[i].lower = u.lower;
        limits->mlim[i].upper = u.upper;


        limits->hlim[i].lower = u.lower;
        limits->hlim[i].upper = u.upper;
    }
    return limits;


}
