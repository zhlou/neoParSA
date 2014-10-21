/*
 * maternal.cpp
 *
 *  Created on: Feb 4, 2013
 *      Author: zhlou
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>

#include "maternal.h"
#include "flyData.h"

using namespace std;

const double maternal::old_divtimes[3] = { 52.0, 28.0, 12.0 };

const double maternal::divtimes1[1] = { 21.1 };
const double maternal::divtimes2[2] = { 33.5, 12.4 };
const double maternal::divtimes3[3] = { 43.0, 21.9, 9.5 };
const double maternal::divtimes4[4] = { 50.8, 29.7, 17.3, 7.8 };

/* division durations */

const double maternal::old_div_duration[3] = { 4.0, 4.0, 4.0 };

const double maternal::div_duration1[1] = { 5.1 };
const double maternal::div_duration2[2] = { 5.1, 3.3 };
const double maternal::div_duration3[3] = { 5.1, 3.3, 3.0 };
const double maternal::div_duration4[4] = { 5.1, 3.3, 3.0, 3.3 };

/* gastrulation times */

const double maternal::old_gast_time = 88.;

const double maternal::gast_time0 = 50.;
const double maternal::gast_time1 = 71.1;
const double maternal::gast_time2 = 83.5;
const double maternal::gast_time3 = 93.;
const double maternal::gast_time4 = 100.8;

/* full division times: including t<0 */

const double maternal::full_divtimes0[TOTAL_DIVS] = { 0.0, -21.1, -33.5, -43.0,
        -51.8, -57.8 };
const double maternal::full_divtimes1[TOTAL_DIVS] = { 21.1, 0.0, -12.4, -21.9,
        -30.7, -36.7 };
const double maternal::full_divtimes2[TOTAL_DIVS] = { 33.5, 12.4, 0.0, -9.5,
        -18.3, -24.3 };
const double maternal::full_divtimes3[TOTAL_DIVS] = { 43.0, 21.9, 9.5, 0.0,
        -8.8, -14.8 };
const double maternal::full_divtimes4[TOTAL_DIVS] = { 50.8, 29.7, 17.3, 7.8,
        -1.0, -7.0 };

/* division durations */

const double maternal::full_div_durations[TOTAL_DIVS] =
        { 5.1, 3.3, 3.0, 3.3, 3.0, 3.0 };


/*** ParseLineage: takes lineage number as input and returns the cleavage **
 *                 cycle the nucleus belongs to.                           *
 ***************************************************************************/

unsigned int ParseLineage(unsigned int lin)
{
    static unsigned int CYCLE1 = 1;
    static unsigned int CYCLE2 = 2;
    static unsigned int CYCLE3 = 4;
    static unsigned int CYCLE4 = 8;
    static unsigned int CYCLE5 = 16;
    static unsigned int CYCLE6 = 32;
    static unsigned int CYCLE7 = 64;
    static unsigned int CYCLE8 = 128;
    static unsigned int CYCLE9 = 256;
    static unsigned int CYCLE10 = 512;
    static unsigned int CYCLE11 = 1024;
    static unsigned int CYCLE12 = 2048;
    static unsigned int CYCLE13 = 4096;
    static unsigned int CYCLE14 = 8192;
    if (lin & CYCLE14)
        return 14;
    else if (lin & CYCLE13)
        return 13;
    else if (lin & CYCLE12)
        return 12;
    else if (lin & CYCLE11)
        return 11;
    else if (lin & CYCLE10)
        return 10;
    else
        error("ParseLineage: illegal lineage number %d", lin);

    return 0; /* just to make the compiler happy! */
}

maternal::maternal(FILE* fp, int divstyle) : olddivstyle(divstyle)
{
    custom_gast = 0.;
    full_nnucs = NULL;
    full_lin_start = NULL;
    bt_init_flag = 0;
    rule_flag = 0;
    ReadTheProblem(fp);
    genotypes = ReadGenotypes(fp);
    nalleles = count_Slist(genotypes);
    InitBicoid(fp);
    InitBias(fp);
    InitNNucs();
    InitGetD();
    InitTheta();
}

maternal::~maternal()
{
    free(lin_start);
    free(full_lin_start);
    free(nnucs);
    free(full_nnucs);
    free(defs.gene_ids);
    if (defs.egenes > 0)
        free(defs.egene_ids);
    int i, j;
    for (i = 0; i < nalleles; ++i) {
        // free bcdtype & biastype
        free(bcdtype[i].genotype);
         for (j = 0; j < bcdtype[i].ptr.bicoid.size; ++j) {
            // next version, this HAS to have a built in destructor
            free(bcdtype[i].ptr.bicoid.array[j].gradient.array);
        }
        free(bcdtype[i].ptr.bicoid.array);

        free(biastype[i].genotype);
        for (j = 0; j < biastype[i].ptr.bias.size; ++j) {
            free(biastype[i].ptr.bias.array[j].state.array);
        }
        free(biastype[i].ptr.bias.array);

        free(bt[i].genotype);
        free(bt[i].ptr.times.array);

    }
    free(bcdtype);
    free(biastype);
    free(bt);
    free_Slist(genotypes);
}

/*** ReadTheProblem: reads the problem section of a data file into the *****
 *                   TheProblem struct.                                    *
 ***************************************************************************/
void maternal::ReadTheProblem(FILE* fp)
{

    fp = FindSection(fp, "problem"); /* find problem section */
    if (!fp)
        error("ReadTheProblem: cannot locate problem section");

    fscanf(fp, "%*s\n"); /* advance pointer past first text line */

    if (1 != (fscanf(fp, "%d\n", &defs.ngenes))) /* read # of genes */
        error("ReadTheProblem: error reading problem section (ngenes)");

    fscanf(fp, "%*s\n"); /* advance the pointer past the second text line */

    if (1 != (fscanf(fp, "%d\n", &defs.egenes))) /* read # of
     external inputs */
        error("ReadTheProblem: error reading problem section (egenes)");

    fscanf(fp, "%*s\n"); /* advance the pointer past the second text line */

    defs.gene_ids = (char *) calloc(defs.ngenes + 1, sizeof(char));
    if (1 != (fscanf(fp, "%s\n", defs.gene_ids))) /* read geneID string */
        error("ReadTheProblem: error reading problem section (gene_ids)");

    fscanf(fp, "%*s\n"); /* advance the pointer past the second text line */

    if (defs.egenes > 0) {

        defs.egene_ids = (char *) calloc(defs.egenes + 1, sizeof(char));
        if (1 != (fscanf(fp, "%s\n", defs.egene_ids))) /* read geneID string */
            error("ReadTheProblem: error reading problem section (egene_ids)");
    } else {
        defs.egene_ids = (char *) calloc(1, sizeof(char));
        defs.egene_ids = (char *)"\0";
        fscanf(fp, "%*s\n"); /* next line (ignore empty external gene ID line) */
    }

    fscanf(fp, "%*s\n"); /* next line (ignore comment) */

    if (1 != (fscanf(fp, "%d\n", &defs.ndivs))) /* read # of cell divs */
        error("ReadTheProblem: error reading problem section (ndivs)");

    fscanf(fp, "%*s\n"); /* advance the pointer past the third text line */
    if (1 != (fscanf(fp, "%d\n", &defs.nnucs)))/* read the max # of nucs */
        error("ReadTheProblem: error reading problem section (nnucs)");

    fscanf(fp, "%*s\n"); /* advance the pointer once more */

    if (1 != (fscanf(fp, "%c\n", &defs.diff_schedule)))/* read diff sched */
        error("ReadTheProblem: error reading problem section (diff. schedule)");
}

/*** ReadGenotypes: This function reads all the genotypes in a datafile & **
 *                  returns an SList with genotype number and pointers to  *
 *                  the corresponding section titles for bias, bcd & facts *
 ***************************************************************************/
Slist* maternal::ReadGenotypes(FILE* fp)
{
    /* Buffer strings for reading:       */
    char biasbuf[MAX_RECORD]; /* section title of bias section     */
    char factsbuf[MAX_RECORD]; /* section title of data section     */
    char matbuf[MAX_RECORD]; /* section title of bcd section      */
    char histbuf[MAX_RECORD]; /* section title of history section      */
    char extbuf[MAX_RECORD]; /* section title of external
     input section      */
    char gtbuf[MAX_RECORD]; /* genotype string                   */

    char *record; /* pointer to current data record    */

    Slist *current; /* holds current element of Slist    */
    Slist *first; /* pointer to first element of Slist */
    Slist *last = NULL; /* pointer to last element of Slist  */

    /*** open the data file and locate genotype section ********************/

    fp = FindSection(fp, "genotypes");
    if (!fp)
        error("ReadGenotypes: cannot locate genotypes");

    if (!(record = (char *) calloc(MAX_RECORD, sizeof(char))))
        error("ReadGenotypes: error allocating record");
    first = init_Slist();
    current = first;

    /*** read titles of bias, data and bcd sections and genotype number for *
     *   each genotype ******************************************************/

    while (strncmp((record = fgets(record, MAX_RECORD, fp)), "$$", 2)) {

        if (6
                != sscanf(record, "%s %s %s %s %s %s", biasbuf, factsbuf,
                        matbuf, histbuf, extbuf, gtbuf))
            error("ReadGenotypes: error reading %s", record);

        /* we eventually want to get rid of the hard wired genotype identifiers    *
         * and replace them by some kind of mechanism by which we can include any  *
         * genes we want in a scoring or annealing run. Think about this!          */

        if (strlen(gtbuf) != defs.ngenes)
            error(
                    "ReadGenotypes: bad genotype string %s (does not match ngenes)",
                    gtbuf);

        if (!current)
            current = init_Slist();

        if (!(current->bias_section = (char *) calloc(MAX_RECORD, sizeof(char))))
            error("ReadGenotypes: error allocating bias_section");
        current->bias_section = strcpy(current->bias_section, biasbuf);

        if (!(current->fact_section = (char *) calloc(MAX_RECORD, sizeof(char))))
            error("ReadGenotypes: error allocating fact_section");
        current->fact_section = strcpy(current->fact_section, factsbuf);

        if (!(current->bcd_section = (char *) calloc(MAX_RECORD, sizeof(char))))
            error("ReadGenotypes: error allocating bcd_section");
        current->bcd_section = strcpy(current->bcd_section, matbuf);

        if (!(current->hist_section = (char *) calloc(MAX_RECORD, sizeof(char))))
            error("ReadGenotypes: error allocating hist_section");
        current->hist_section = strcpy(current->hist_section, histbuf);

        if (!(current->ext_section = (char *) calloc(MAX_RECORD, sizeof(char))))
            error("ReadGenotypes: error allocating ext_section");
        current->ext_section = strcpy(current->ext_section, extbuf);

        if (!(current->genotype = (char *) calloc(MAX_RECORD, sizeof(char))))
            error("ReadGenotypes: error allocating genotype string");
        current->genotype = strcpy(current->genotype, gtbuf);

        addto_Slist(last, current);

        last = current;
        current = NULL;
    }

    free(record);

    return first;
}

/*** InitBicoid: Copies the Blist read by ReadBicoid into the DArrPtr ******
 *               structure; the bicoid DArrPtr contains pointers to bcd    *
 *               arrays for each cell division cycle                       *
 ***************************************************************************/
void maternal::InitBicoid(FILE* fp)
{
    int i; /* local loop counter */

    Blist *inlist; /* temporary linked list for bcd */

    Slist *current; /* types from data file */

    if (!(bcdtype = (GenoType *) calloc(nalleles, sizeof(GenoType))))
        error("InitBicoid: Could not allocate bcdtype struct");

    /*** for loop: read bicoid for each genotype *******************************/

    for (current = genotypes, i = 0; current; current = current->next, i++) {

        if (!(inlist = ReadBicoid(fp, current->bcd_section))) /* read bicoid */
            error("InitBicoid: error reading %s", current->bcd_section);
        else {

            if (!(bcdtype[i].genotype = (char *) calloc(MAX_RECORD,
                    sizeof(char))))
                error("InitBicoid: could not allocate bcd genotype string");
            bcdtype[i].genotype = strcpy(bcdtype[i].genotype,
                    current->genotype);
            bcdtype[i].ptr.bicoid = List2Bicoid(inlist);

            free_Blist(inlist);
        }
    }
}

/*** InitBias:  puts bias records in a form where get_bias can use them; ***
 *              it expects times in increasing order; it expects a non-    *
 *              sparse entry, with no genes or nuclei missing              *
 ***************************************************************************/
void maternal::InitBias(FILE* fp)
{
    int i; /* loop counters */

    Dlist *inlist; /* temporary linked list for bias */

    Slist *current; /* types from data file */

    int ndp = 0; /* dummy for ReadData, no need to count datapts here */

    if (!(biastype = (GenoType *) calloc(nalleles, sizeof(GenoType))))
        error("InitBias: Could not allocate biastype struct");

    /*** for loop: read bicoid for each genotype *******************************/

    for (current = genotypes, i = 0; current; current = current->next, i++) {

        if (!(inlist = ReadData(fp, current->bias_section, &ndp, defs.ngenes))) /* read bias */
            error("InitBias: error reading %s", current->bias_section);
        else {

            if (!(biastype[i].genotype = (char *) calloc(MAX_RECORD,
                    sizeof(char))))
                error("InitBias: could not allocate bias genotype string");
            biastype[i].genotype = strcpy(biastype[i].genotype,
                    current->genotype);
            biastype[i].ptr.bias = List2Bias(inlist);

            free_Dlist(inlist);
        }
    }

    InitBTs(); /* initialize the static bias time struct */

}

/*** InitBTs: initializes the static BT struct that holds all times for ****
 *            which we have bias.                                          *
 ***************************************************************************/
void maternal::InitBTs(void)
{
    int i, j;

    if (!bt_init_flag) {
        if (!(bt = (GenoType *) calloc(nalleles, sizeof(GenoType))))
            error("InitBTs: could not allocate bt struct");

        for (i = 0; i < nalleles; i++) {

            if (!(bt[i].genotype = (char *) calloc(MAX_RECORD, sizeof(char))))
                error("InitBTs: could not allocate BT genotype string");
            bt[i].genotype = strcpy(bt[i].genotype, biastype[i].genotype);
            bt[i].ptr.times.size = biastype[i].ptr.bias.size;
            if (!(bt[i].ptr.times.array = (double *) calloc(
                    bt[i].ptr.times.size, sizeof(double))))
                error("InitBTs: could not allocate bt array");

            for (j = 0; j < bt[i].ptr.times.size; j++)
                bt[i].ptr.times.array[j] = biastype[i].ptr.bias.array[j].time;
        }

        bt_init_flag = 1;
    }
}

/*** InitNNucs: takes the global defs.nnucs and calculates number of nucs **
 *              for each cleavage cycle which are then stored in reverse   *
 *              order in the static nnucs[] array                          *
 *   CAUTION:   defs struct and lin_start need to be initialized before!   *
 ***************************************************************************/
void maternal::InitNNucs(void)
{
    int i; /* loop counter */
    int n; /* used to calculate number of nucs for each cell cycle */

    if (!(nnucs = (int *) calloc(defs.ndivs + 1, sizeof(int))))
        error("InitNNucs: could not allocate nnucs array");

    n = defs.nnucs;

    /* below we have to take into account two cases: a) most anterior lineage  *
     * number is odd-numbered -> always add a nucleus to the earlier cycle or  *
     * b) most anterior lineage number is even-numbered: just add an additio-  *
     * nal nucleus if last nucleus is odd-numbered (see also exhaustive com-   *
     * ments about this at the DIVIDE rule in Blastoderm() in integrate.c)     */

    for (i = 0; i <= defs.ndivs; i++) {
        nnucs[i] = n;
        if (lin_start[i] % 2)
            n = n / 2 + 1;
        else
            n = (n % 2) ? n / 2 + 1 : n / 2;
    }
}


/*** ReadBicoid: reads the bcd section of a data file into a linked list; **
 *               also determines maxconc from the bicoid gradient          *
 ***************************************************************************/
Blist* maternal::ReadBicoid(FILE* fp, char* section)
{
    int c; /* holds char for parser          */
    int lead_punct; /* flag for leading punctuation   */

    double maxv = -1.; /* maximum v (protein conc.)      */

    char *base; /* pointer to beginning of string */
    char *record; /* pointer to string, used as     */
    /* counter                        */
    Blist *current; /* holds current element of Blist */
    Blist *inlist; /* holds whole read Blist         */

    if ((fp = FindSection(fp, section))) { /* position the fp */

        base = (char *) calloc(MAX_RECORD, sizeof(char));
        current = NULL;
        inlist = NULL;

        /* while loop: reads and processes lines from file until sections ends *****/

        while (strncmp((base = fgets(base, MAX_RECORD, fp)), "$$", 2)) {

            record = base; /* base always points to start */
            lead_punct = 0; /* of string */

            /* for loop: parses and processes each line from the data file *************/

            c = (int) *record;
            while (c != '\0') {

                if (isdigit(c)) { /* number means data */

                    record = base; /* reset pointer to start of str */
                    current = init_Blist(); /* allocate memory for lnkd list */

                    if (2
                            != sscanf(record, "%d %lg", &(current->lineage),
                                    &(current->conc)))
                        error("ReadBicoid: error reading %s", base);

                    if (current->conc > maxv)
                        maxv = current->conc;

                    inlist = addto_Blist(inlist, current);
                    break;

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
                    error("ReadBicoid: illegal character in %s", base);
                }
            }
        }

        if (maxv > 12.) /* oldstyle or newstyle data? */
            maxconc = 255.;
        else
            maxconc = 12.;

        free(base);
        return inlist;

    } else {

        return NULL;
    }
}

/*** List2Bicoid: takes a Blist and returns the corresponding BArrPtr ******
 *                structure; also initializes the lin_start array which    *
 *                contains the lineage number of the most anterior nucleus *
 *                for each cell cycle (we need this for cell division and  *
 *                printing model output)                                   *
 ***************************************************************************/
BArrPtr maternal::List2Bicoid(Blist* inlist)
{
    int i = 0; /* local loop counter */
    int n; /* used to evaluate # of nuclei */
    int lin_count = 0; /* counter for initializing lin_start */

    unsigned int ccycle; /* what cleavage cycle are we in? */
    unsigned int samecycle = 0; /* same cleavage cycle as before? */

    Blist *current; /* current element of Blist */
    Blist *start; /* used to evaluate # of nuclei */

    BArrPtr bicoid; /* BArrPtr to be returned */

    bicoid.size = 0;
    bicoid.array = NULL;

    if (!(lin_start = (int *) calloc(defs.ndivs + 1, sizeof(int))))
        error("List2Bicoid: could not allocate lin_start array");




    /*** for loop: step through linked list and copies values into an array    *
     *             of BcdGrads; there's one BcdGrad for each cleavage cycle;   *
     *             each BcdGrad has: - ccycle (the cleavage cycle number)      *
     *                               - gradient.size (# of nuclei for grad.)   *
     *                               - gradient.array (pointer to array)       *
     ***************************************************************************/

    for (current = inlist; current; current = current->next) {

        ccycle = ParseLineage(current->lineage); /* which cycle is it? */

        /* allocate new gradient struct for new cell cycle */

        if (ccycle != samecycle) {

            samecycle = ccycle;

            bicoid.size++;
            bicoid.array = (BcdGrad *) realloc(bicoid.array,
                    bicoid.size * sizeof(BcdGrad));

            /* this loop determines the number of nuclei in each cycle */

            n = 0;
            for (start = current; ParseLineage(current->lineage) == samecycle;
                    current = current->next) {
                n++;
                if (!(current->next)) /* don't count garbage -> core dumps... */
                    break;
            }
            current = start; /* reset pointer to where we were */

            /* initialize lin_start: this array is later used by Blastoderm and such   */

            lin_start[defs.ndivs - lin_count] = current->lineage;
            lin_count++;

            /* allocate array for gradient here and reset the counter */

            bicoid.array[bicoid.size - 1].ccycle = samecycle; /* next three */
            /* lines define BcdGrad for each cleavage cycle */
            bicoid.array[bicoid.size - 1].gradient.array = (double *) calloc(n,
                    sizeof(double));
            bicoid.array[bicoid.size - 1].gradient.size = n;
            i = 0;

        }


        /* in any case: read concentration into gradient array */

        bicoid.array[bicoid.size - 1].gradient.array[i] = current->conc;
        i++;
    }

    // InitFullNNucs & part of List2Interp stuff here
    if ( !(full_lin_start = (int *)calloc(defs.ndivs+1, sizeof(int))) )
        error("List2Interp: could not allocate full_lin_start array");

    full_ccycles =  defs.ndivs + 1;


    if ( !(full_nnucs = (int *)calloc(full_ccycles, sizeof(int))) )
      error("InitFullNNucs: could not allocate full_nnucs array");

    n = defs.nnucs;

  /* below we have to take into account two cases: a) most anterior lineage  *
   * number is odd-numbered -> always add a nucleus to the earlier cycle or  *
   * b) most anterior lineage number is even-numbered: just add an additio-  *
   * nal nucleus if last nucleus is odd-numbered (see also exhaustive com-   *
   * ments about this at the DIVIDE rule in Blastoderm() in integrate.c)     */

    for(i=0; i < full_ccycles; i++) {
        full_lin_start[i] = lin_start[i];
        full_nnucs[i] = n;
        if ( full_lin_start[i] % 2 )
            n = n/2 + 1;
        else
            n = (n % 2) ? n/2 + 1 : n/2;
    }

    return bicoid;
}


/*** List2Bias: takes a Dlist and returns the corresponding DArrPtr ********
 *              structure.                                                 *
 ***************************************************************************/
NArrPtr maternal::List2Bias(Dlist* inlist)
{
    int i = 0;
    int j; /* local loop counters */

    int n; /* used to evaluate # of nuclei */
    double now = -999999999.0; /* variable for time */

    Dlist *current; /* holds current element of Dlist */
    Dlist *start; /* pointer used for evaluating # of */

    NArrPtr bias; /* DArrPtr to be returned */

    bias.size = 0;
    bias.array = NULL;

    /*** for loop: steps through linked list and copies values into an array   *
     *             of NucStates. There's one NucState for each time step       *
     *             Each NucState has: - time (the time)                        *
     *                                - state.size (# of genes * # of nucs)    *
     *                                - state.array (pointer to array)         *
     ***************************************************************************/

    for (current = inlist; current; current = current->next) {

        if (current->d[0] != now) { /* a new time point: allocate */
            now = current->d[0]; /* the time is now! */
            bias.size++; /* one NucState for each time step */
            bias.array = /* allocate for one more NucState */
            (NucState *) realloc(bias.array, bias.size * sizeof(NucState));

            /* determine number of nuclei per time step */

            n = 0;
            for (start = current; current->d[0] == now; current =
                    current->next) {
                n++;
                if (!(current->next)) /* don't count garbage and cause dis- */
                    break; /* may and core dumps... */
            }
            current = start; /* reset list ptr to where we were */

            /* allocate a bias array for each biastime */

            bias.array[bias.size - 1].time = now;
            bias.array[bias.size - 1].state.array = (double *) calloc(
                    n * defs.ngenes, sizeof(double));
            bias.array[bias.size - 1].state.size = n * defs.ngenes;

            i = 0;
        }

        /* always: read concs into array */

        for (j = 1; j <= defs.ngenes; j++) {
            bias.array[bias.size - 1].state.array[i] = current->d[j];
            i++;
        }
    }

    return bias;
}

void maternal::InitGetD() // set pointer to the right division table
{
    if (olddivstyle) {
        if (defs.ndivs != 3)
            error("GetD: only 3 cell divisions allowed for oldstyle (-o)");
        getD_table = (double *) old_divtimes;
    } else if (defs.ndivs == 0) {
    } /* no need for table, if there are no cell divs */
    else if (defs.ndivs == 1)
        getD_table = (double *) divtimes1;
    else if (defs.ndivs == 2)
        getD_table = (double *) divtimes2;
    else if (defs.ndivs == 3)
        getD_table = (double *) divtimes3;
    else if (defs.ndivs == 4)
        getD_table = (double *) divtimes4;
    else
        error("GetD: can't handle %d cell divisions!", defs.ndivs);
}

DArrPtr maternal::GetBicoid(double time, int genindex)
{
    int             i;
    unsigned int    ccycle;

    if ( genindex < 0 || genindex >= nalleles )
      error("GetBicoid: invalid genotype index %d", genindex);

    ccycle = GetCCycle(time);

    for (i=0; i <= bcdtype[genindex].ptr.bicoid.size; i++)
      if ( ccycle == bcdtype[genindex].ptr.bicoid.array[i].ccycle )
        return bcdtype[genindex].ptr.bicoid.array[i].gradient;

    error("GetBicoid: no bicoid gradient for ccycle %d", ccycle);

    return bcdtype[genindex].ptr.bicoid.array[i].gradient;
                                         /* just to make the compiler happy! */
}

/*** GetBias: This function returns bias values for a given time and *******
 *            genotype.                                                    *
 ***************************************************************************/


DArrPtr maternal::GetBias(double time, int genindex)
{
  int      i;

  if ( genindex < 0 || genindex >= nalleles )
    error("GetBias: invalid genotype index %d", genindex);

  for(i=0; i < biastype[genindex].ptr.bias.size; i++) {
    if (biastype[genindex].ptr.bias.array[i].time == time)
      break;
  }

  return biastype[genindex].ptr.bias.array[i].state;
}

/*** GetD: returns diffusion parameters D according to the diff. params. ***
 *         in the data file and the diffusion schedule used                *
 *   NOTE: Caller must allocate D_tab                                      *
 ***************************************************************************/
// table has been moved to maternal class, and its initialization become
// init_GetD
void maternal::GetD(double t, double* d, char diff_schedule, double* D_tab)
{

    int i; /* loop counter */
    double gast; /* gastrulation time */
    double cutoff; /* cutoff time */
    double lscale = 1; /* scaling factor */

    /* first time GetD is called: set pointer to the right division table
     * now in InitGetD */

    /* this loop takes lscale square for each cell division, i.e. the earlier  *
     * we are the bigger lscale (and the smaller the D's that we return        */

    for (i = 0; i < defs.ndivs; i++)
        if (t < getD_table[i])
            lscale *= 2;

    /* diffusion schedule A: all Ds always the same */

    if (diff_schedule == 'A')
        for (i = 0; i < defs.ngenes; i++)
            D_tab[i] = d[0];

    /* diffusion schedule B: all genes have different D's that depend on in-   *
     * verse l-square                                                          */

    else if (diff_schedule == 'B')
        for (i = 0; i < defs.ngenes; i++)
            D_tab[i] = d[i] / (lscale * lscale);

    /* diffusion schedule C: all genes have the same D that depends on inverse *
     * l-square                                                                */

    else if (diff_schedule == 'C')
        for (i = 0; i < defs.ngenes; i++)
            D_tab[i] = d[0] / (lscale * lscale);

    /* diffusion schedule D: all genes have different D's which don't change   *
     * over time                                                               */

    else if (diff_schedule == 'D')
        for (i = 0; i < defs.ngenes; i++)
            D_tab[i] = d[i];

    /* diffusion schedule E: used cutoff at gast-12 otherwise just like B */

    else if (diff_schedule == 'E') {
        if (olddivstyle) {
            if (defs.ndivs != 3)
                error("GetD: only 3 cell divisions allowed for oldstyle (-o)");
            gast = old_gast_time;
        } else if (defs.ndivs == 0)
            gast = gast_time0;
        else if (defs.ndivs == 1)
            gast = gast_time1;
        else if (defs.ndivs == 2)
            gast = gast_time2;
        else if (defs.ndivs == 3)
            gast = gast_time3;
        else if (defs.ndivs == 4)
            gast = gast_time4;
        else
            error("GetD: can't handle %d cell divisions!", defs.ndivs);

        cutoff = gast - 12.0; /* This value probably wrong; see Merrill88 */
        for (i = 0; i < defs.ngenes; i++)
            D_tab[i] = (t < cutoff) ? d[i] / (lscale * lscale) : 0.;

        /* any other diffusion schedule: error! */

    } else
        error("GetD: no code for schedule %c!", diff_schedule);
}
/*** GetCCycle: returns cleavage cycle number for a given time *************
 ***************************************************************************/
unsigned int maternal::GetCCycle(double time)
{
    int i; /* loop counter */
    double *table; /* local copy of divtimes table */

    /* assign 'table' to the appropriate division schedule */

    if (olddivstyle) {
        if (defs.ndivs != 3)
            error("GetCCycle: only 3 cell divisions allowed for oldstyle (-o)");
        table = (double *) old_divtimes;
    } else if (defs.ndivs == 0)
        return 14; /* if defs.ndivs == 0: we're always in cycle 14 */
    else if (defs.ndivs == 1)
        table = (double *) divtimes1;
    else if (defs.ndivs == 2)
        table = (double *) divtimes2;
    else if (defs.ndivs == 3)
        table = (double *) divtimes3;
    else if (defs.ndivs == 4)
        table = (double *) divtimes4;
    else
        error("GetCCycle: can't handle %d cell divisions!", defs.ndivs);

    /* evaluate number of cell cycle for current time; note that for the exact *
     * time of cell division, we'll return the number of the previous cell cy- *
     * cyle                                                                    */

    for (i = 0; i < defs.ndivs; i++)
        if (time > table[i])
            return 14 - i;

    return 14 - (i++);
}

void maternal::InitTheta()
{
    // These 2 static variables have become member variables
    //static double *theta_dt; /* pointer to division time table */
    //static double *theta_dd; /* pointer to division duration table */
    if (olddivstyle) {
        /* get pointers to division time table */
        theta_dt = (double*) (old_divtimes); /* and division duration table */
        theta_dd = (double*) (old_div_duration);
    } else {
        if (defs.ndivs == 0) {
            theta_dt = (double*) (full_divtimes0);
            theta_dd = (double*) (full_div_durations);
        } else if (defs.ndivs == 1) {
            theta_dt = (double*) (full_divtimes1);
            theta_dd = (double*) (full_div_durations);
        } else if (defs.ndivs == 2) {
            theta_dt = (double*) (full_divtimes2);
            theta_dd = (double*) (full_div_durations);
        } else if (defs.ndivs == 3) {
            theta_dt = (double*) (full_divtimes3);
            theta_dd = (double*) (full_div_durations);
        } else if (defs.ndivs == 4) {
            theta_dt = (double*) (full_divtimes4);
            theta_dd = (double*) (full_div_durations);
        } else
            error("Theta: can't handle %d cell divisions!", defs.ndivs);
    }
}

int maternal::Theta(double time)
{
    int i;

    // These 2 static variables have become member variables
    //static double *theta_dt; /* pointer to division time table */
    //static double *theta_dd; /* pointer to division duration table */





    /* checks if we're in a mitosis; we need the 10*DBL_EPSILON kludge for gcc *
     * on Linux which can't handle truncation errors very well                 */

    for (i = 0; i < TOTAL_DIVS; i++) {
        /*      printf("Thetatime=%.16f,[%.16f,%.16f]\n",time,(*(theta_dt+i) - *(theta_dd+i) + HALF_EPSILON), (*(theta_dt+i) + HALF_EPSILON)); */
        if ((time <= (*(theta_dt + i) + HALF_EPSILON))
                && (time >= (*(theta_dt + i) - *(theta_dd + i) + HALF_EPSILON)))
            return MITOSIS;
    }

    return INTERPHASE;
}

int maternal::GetNNucs(double t) const
{
    int      i;                                              /* loop counter */
    double   *table;                         /* local copy of divtimes table */

  /* assign 'table' to the appropriate division schedule */

    if ( olddivstyle ) {
      if ( defs.ndivs != 3 )
        error("GetNNucs: only 3 cell divisions allowed for oldstyle (-o)");
      table = (double *)old_divtimes;
    } else
      if ( defs.ndivs == 0 )
        return nnucs[0];
      else if ( defs.ndivs == 1 )
        table = (double *)divtimes1;
      else if ( defs.ndivs == 2 )
        table = (double *)divtimes2;
      else if ( defs.ndivs == 3 )
        table = (double *)divtimes3;
      else if ( defs.ndivs == 4 )
        table = (double *)divtimes4;
      else
        error("GetNNucs: can't handle %d cell divisions!", defs.ndivs);

  /* evaluate nnucs for current time; note that for the *exact* time of cell *
   * division, we'll return the number of nuclei before the division has ac- *
   * tually occurred                                                         */

    for (i=0; i<defs.ndivs; i++)
      if ( t>table[i] )
        return nnucs[i];

    return nnucs[i];
}

/*** GetGastTime: returns time of gastrulation depending on ndivs and ******
 *                olddivstyle; returns 0 in case of an error; if a custom  *
 *                gastrulation time is chosen with -S, it will be returned *
 *                only if it's bigger than the normal gastrulation time    *
 ***************************************************************************/
double maternal::GetGastTime(void) const
{
    if (olddivstyle) {
        if (defs.ndivs != 3)
            error("GetDurations: only 3 cell divisions allowed for oldstyle (-o)");

        if (custom_gast > old_gast_time)
            return custom_gast;
        else
            return old_gast_time;

    } else

    if (defs.ndivs == 0)

        if (custom_gast > gast_time0)
            return custom_gast;
        else
            return gast_time0;

    else if (defs.ndivs == 1)

        if (custom_gast > gast_time1)
            return custom_gast;
        else
            return gast_time1;

    else if (defs.ndivs == 2)

        if (custom_gast > gast_time2)
            return custom_gast;
        else
            return gast_time2;

    else if (defs.ndivs == 3)

        if (custom_gast > gast_time3)
            return custom_gast;
        else
            return gast_time3;

    else if (defs.ndivs == 4)

        if (custom_gast > gast_time4)
            return custom_gast;
        else
            return gast_time4;

    else
        error("GetGastTime: can't handle %d cell divisions!", defs.ndivs);

    return 0;
}

/*** GetBTimes: returns a sized array of times for which there is bias *****
 ***************************************************************************/

DArrPtr maternal::GetBTimes(char *genotype) const
{
    int index;                                               /* loop counter */

    for(index=0; index<nalleles; index++)
        if ( !(strcmp(biastype[index].genotype, genotype)) )
            break;

    /* if no explicit bias times for this genotype -> use wt bias times */

    if ( index == nalleles )
        index = 0;

    /* check if we actually have biastimes at all or otherwise -> error! */

    if ( bt_init_flag )
        return bt[index].ptr.times;
    else
        error("GetBTimes: called without initialized BTs");

    return bt[index].ptr.times;
    /* just to make the compiler happy */
}

/*** GetDivtable: returns times of cell divisions depending on ndivs and ***
 *                olddivstyle; returns NULL in case of an error            *
 ***************************************************************************/

double *maternal::GetDivtable(void)
{
    if ( olddivstyle ) {
        if ( defs.ndivs != 3 )
            error("GetDivtable: only 3 cell divisions allowed for oldstyle (-o)");
        return (double *)old_divtimes;
    } else
        if ( defs.ndivs == 0 )
            return NULL;                           /* return error if ndivs == 0 */
        else if ( defs.ndivs == 1 )
            return (double *)divtimes1;
        else if ( defs.ndivs == 2 )
            return (double *)divtimes2;
        else if ( defs.ndivs == 3 )
            return (double *)divtimes3;
        else if ( defs.ndivs == 4 )
            return (double *)divtimes4;
        else
            error
            ("GetDivtable: can't handle %d cell divisions!", defs.ndivs);

    return NULL;
}

/*** GetDurations: returns pointer to durations of cell divisions de- ******
 *                 pending on ndivs and olddivstyle; returns NULL in case  *
 *                 of an error                                             *
 ***************************************************************************/

double *maternal::GetDurations(void)
{
  if ( olddivstyle ) {
    if ( defs.ndivs != 3 )
      error("GetDurations: only 3 cell divisions allowed for oldstyle (-o)");
    return (double *)old_div_duration;
  } else
    if ( defs.ndivs == 0 )
      return NULL;                           /* return error if ndivs == 0 */
    else if ( defs.ndivs == 1 )
      return (double *)div_duration1;
    else if ( defs.ndivs == 2 )
      return (double *)div_duration2;
    else if ( defs.ndivs == 3 )
      return (double *)div_duration3;
    else if ( defs.ndivs == 4 )
      return (double *)div_duration4;
    else
      error("GetDurations: can't handle %d cell divisions!", defs.ndivs);

  return NULL;
}

/*** GetRule: returns the appropriate rule for a given time; used by the ***
 *            derivative function                                          *
 ***************************************************************************/

int maternal::GetRule(double time)
{
  int i;

  static double *dt;                     /* pointer to division time table */
  static double *dd;                 /* pointer to division duration table */



/* the following block of code are implemented like this (instead of       */
/* calling GetDivtable or GetDurations) to save function calls and increa- */
/* se performance (GetRule is called from within the inner loop)           */


  if ( !rule_flag ) {                                 /* only do this once */
    if ( olddivstyle ) {            /* get pointers to division time table */
      dt=(double *)old_divtimes;            /* and division duration table */
      dd=(double *)old_div_duration;
    } else {
      if ( defs.ndivs == 0 ) {
    return INTERPHASE;                /* no cell division? no MITOSIS! */
      } else if ( defs.ndivs == 1 ) {
    dt=(double *)divtimes1;
    dd=(double *)div_duration1;
      } else if ( defs.ndivs == 2 ) {
    dt=(double *)divtimes2;
    dd=(double *)div_duration2;
      } else if ( defs.ndivs == 3 ) {
    dt=(double *)divtimes3;
    dd=(double *)div_duration3;
      } else if ( defs.ndivs == 4 ) {
    dt=(double *)divtimes4;
    dd=(double *)div_duration4;
      } else
    error("GetRule: can't handle %d cell divisions!", defs.ndivs);
    }
    rule_flag=1;
  }

/* checks if we're in a mitosis; we need the 10*DBL_EPSILON kludge for gcc *
 * on Linux which can't handle truncation errors very well                 */

  for(i=0; i<defs.ndivs; i++)
  {
/*      printf("GetRuletime=%.16f,[%.16f,%.16f]\n",time,(*(dt+i) - *(dd+i) + HALF_EPSILON), (*(dt+i) + HALF_EPSILON)); */
    if ((time <= (*(dt+i) + HALF_EPSILON))
    && (time >= (*(dt+i) - *(dd+i) + HALF_EPSILON)))
      return MITOSIS;
  }

  return INTERPHASE;
}

/*** GetStartLin: returns the lineage number of the most anterior nucleus **
 *                for a given time                                         *
 ***************************************************************************/

int maternal::GetStartLin(double t)
{
  int      i;                                              /* loop counter */
  double   *table;                         /* local copy of divtimes table */

/* assign 'table' to the appropriate division schedule */

  if ( olddivstyle ) {
    if ( defs.ndivs != 3 )
      error("GetStartLin: only 3 cell divisions allowed for oldstyle (-o)");
    table = (double *)old_divtimes;
  } else
    if ( defs.ndivs == 0 )
      return lin_start[0];
    else if ( defs.ndivs == 1 )
      table = (double *)divtimes1;
    else if ( defs.ndivs == 2 )
      table = (double *)divtimes2;
    else if ( defs.ndivs == 3 )
      table = (double *)divtimes3;
    else if ( defs.ndivs == 4 )
      table = (double *)divtimes4;
    else
      error("GetStartLin: can't handle %d cell divisions!", defs.ndivs);

/* evaluate lineage number of most anterior nucleus for current time; note *
 * that for the *exact* time of cell division, we'll return the lineage    *
 * number of the most anterior nucleus of the previous cell cycle          */

  for (i=0; i<defs.ndivs; i++)
    if ( t>table[i] )
      return lin_start[i];

  return lin_start[i];
}

int maternal::GetStartLinIndex(double t)
{
  int      i;                                              /* loop counter */
  double   *table;                         /* local copy of divtimes table */

/* assign 'table' to the appropriate division schedule */

  if ( olddivstyle ) {
    if ( defs.ndivs != 3 )
      error("GetStartLin: only 3 cell divisions allowed for oldstyle (-o)");
    table = (double *)old_divtimes;
  } else
    if ( defs.ndivs == 0 )
      table = (double *)full_divtimes0;
    else if ( defs.ndivs == 1 )
      table = (double *)full_divtimes1;
    else if ( defs.ndivs == 2 )
      table = (double *)full_divtimes2;
    else if ( defs.ndivs == 3 )
      table = (double *)full_divtimes3;
    else if ( defs.ndivs == 4 )
      table = (double *)full_divtimes4;
    else
      error("GetStartLin: can't handle %d cell divisions!", defs.ndivs);

/* evaluate lineage number of most anterior nucleus for current time; note *
 * that for the *exact* time of cell division, we'll return the lineage    *
 * number of the most anterior nucleus of the previous cell cycle          */

  for (i=0; i<full_ccycles-1; i++)
    if ( t > table[i]+HALF_EPSILON )
      return i;

  return i;
}

double *maternal::Get_Theta_Discons(int *theta_discon_size)
{

    int s,i,j;
    double *dt,*dd,*disc_array;

    s = (TOTAL_DIVS - defs.ndivs);
    disc_array = (double *) calloc(2*s, sizeof(double));

    if ( defs.ndivs == 0 ) {
        dt=(double *)full_divtimes0;
        dd=(double *)full_div_durations;
      } else if ( defs.ndivs == 1 ) {
        dt=(double *)full_divtimes1;
        dd=(double *)full_div_durations;
      } else if ( defs.ndivs == 2 ) {
        dt=(double *)full_divtimes2;
        dd=(double *)full_div_durations;
      } else if ( defs.ndivs == 3 ) {
        dt=(double *)full_divtimes3;
        dd=(double *)full_div_durations;
      } else if ( defs.ndivs == 4 ) {
        dt=(double *)full_divtimes4;
        dd=(double *)full_div_durations;
      } else
        error("Get_Theta_Discons: can't handle %d cell divisions!", defs.ndivs);

    for (i=0,j=0;i<s;i++,j+=2) {
        *(disc_array+j) = dt[defs.ndivs+i];
        *(disc_array+j+1) = dt[defs.ndivs+i] -
        dd[defs.ndivs+i];
    }

    *theta_discon_size = 2*s;
    return disc_array;
}

DataTable *maternal::List2Interp(Dlist *inlist, int num_genes)
{


  int       i = 0;
  int       j;                                      /* local loop counters */

  double    now = -999999999.;            /* assigns data to specific time */

  Dlist     *current;                    /* holds current element of Dlist */

  DataTable *D;                                 /* local copy of DataTable */



  D = (DataTable *)malloc(sizeof(DataTable));
                                         /* Initialize DataTable structure */
  D->size = 0;
  D->record = NULL;

/*** for loop: steps through linked list and transfers facts into Data-    *
 *             Records, one for each time step                             *
 ***************************************************************************/

  for (current=inlist; current; current=current->next) {

    if ( current->d[0] != now ) {             /* a new time point: allocate */
      now = current->d[0];                             /* the time is now! */
      D->size++;                           /* one DataRecord for each time */
      D->record =                                   /* allocate DataRecord */
    (DataRecord *)realloc(D->record,D->size*sizeof(DataRecord));

      D->record[D->size-1].time = now;          /* next three lines define */
      D->record[D->size-1].size = 0;                /* DataRecord for each */
      D->record[D->size-1].array = NULL;                      /* time step */
      i = 0;
    }

    for(j=1; j <= num_genes; j++) {     /* always: read concs into array */
    D->record[D->size-1].size++;            /* one more in this record */
    D->record[D->size-1].array =       /* reallocate memory for array! */
      (DataPoint *)realloc(D->record[D->size-1].array,
          D->record[D->size-1].size * sizeof(DataPoint));

/* the following two lines assign concentration value and index ************/

    D->record[D->size-1].array[D->record[D->size-1].size-1].conc =
      current->d[j];
    D->record[D->size-1].array[D->record[D->size-1].size-1].index = i;

      i++;

    }

/* initialize lin_start: this array is later used by Blastoderm and such   */

    if (ParseLineage(full_lin_start[full_ccycles - 1]) !=
            ParseLineage(current->lineage)  ) {
        full_ccycles++;

        if ( !(full_lin_start = (int *)realloc(full_lin_start,
                                            full_ccycles*sizeof(int))) )
            error("List2Interp: could not allocate full_lin_start array");

        full_lin_start[full_ccycles-1] = current->lineage;

    }

  }

    qsort((void *) full_lin_start, full_ccycles, sizeof(int),(int (*)
    (const void*,const void*)) descend);

/*  for (i=0; i < full_ccycles; i++)
    printf("History lineages before removing dups %d\n", full_lin_start[i]);*/

/* Now lets remove duplicates */

    i = 0;

    while (i < full_ccycles-1 )
    {
        if (full_lin_start[i] == full_lin_start[i+1])
        {
            memmove((full_lin_start+i),(full_lin_start+i+1),
                                (full_ccycles-i-1)*sizeof(int));
/*          printf("Shifted %d elements to %d\n",full_ccycles-i-1,i);*/
            full_ccycles--;
            i--;
        }

    i++;
    full_lin_start = (int *) realloc(full_lin_start,
                                    full_ccycles*sizeof(int));
    }

  return D;
}
