/*
 * maternal.cpp
 *
 *  Created on: Feb 4, 2013
 *      Author: zhlou
 */

#include "flyUtils.h"
#include <cstdio>

using namespace std;

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

maternal::maternal(FILE* fp)
{
    ReadTheProblem(fp);
    genotypes = ReadGenotypes(fp);
    nalleles = count_Slist(genotypes);
    InitBicoid(fp);
    InitBias(fp);
    InitNNucs();
}

maternal::~maternal()
{
    // TODO: free biastype
    // TODO: free lin_start
    // TODO: free bcdtype
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
        defs.egene_ids = "\0";
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
