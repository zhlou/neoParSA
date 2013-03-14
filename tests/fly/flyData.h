/*
 * flyData.h
 *
 *  Created on: Feb 14, 2013
 *      Author: zhlou
 */

#ifndef FLYDATA_H_
#define FLYDATA_H_
#include <cstdio>
#include <cfloat>
using namespace std;
/* DBL_EPSILON is about 2 x 10^-16. so BIG_EPSILON is ~10^-11 min. */

const double EPSILON = DBL_EPSILON * 1000.;
const double BIG_EPSILON = EPSILON * 1000.;
const double HALF_EPSILON = EPSILON * .5;

                                  /* This assumes all times < 100. min. !! */

#ifdef       FLOAT_EVERYTHING
#undef       EPSILON
#define      EPSILON     (FLT_EPSILON * 100.)             /* see above */
#endif


const int INTERPHASE = 0;
const int MITOSIS = 1;

const int MAX_RECORD = 256;
const int MAX_ARGS = 25;



/* general struct used for sized array of doubles */
struct DArrPtr
{
    int size;
    double *array;
};

/* two structs for bicoid gradients */

struct BcdGrad {
  int      ccycle;                                   /* the cleavage cycle */
  DArrPtr  gradient;                  /* the array of concs for a gradient */
};

struct BArrPtr {      /* points to an array of BcdStates of 'size' */
  int      size;
  BcdGrad  *array;
};

/* two structures for bias data and Blastoderm() output (solution) */

struct NucState {
  double   time;                                               /* the time */
  DArrPtr  state;                /* the array of v's, and how many of them */
};

struct NArrPtr {
  int      size;                   /* How many NucState elements in array? */
  NucState *array;    /* points to 1st element of array of NucState elems. */
};

/* The following three structs are used to store facts data. ***************
 * The idea is to encode sparse data against non-sparse v[index] state     *
 * arrays. This will in most of the present cases use more memory, but     *
 * that won't always be true. Also, it is more convenient for those        *
 *'don't care' points.                                                     *
 * NOTE: needs to be in this generic header for union definition below     *
 ***************************************************************************/

struct DataPoint {
  int        index;
  double     conc;
};

struct DataRecord {
  int        size;
  double     time;
  DataPoint  *array;
};

struct DataTable {
  int        size;
  DataRecord *record;
};

/* struct and union for genotype number and pointer to data ****************/
/* used for bicoid, bias, facts and time tables                            */

union DataPtr {
  DArrPtr          times;                   /* used for bias and tab times */
  BArrPtr          bicoid;                  /* for bicoid stuff            */
  NArrPtr          bias;                    /* for the bias                */
  DataTable        *facts;                  /* for facts                   */
};

struct GenoType {
  char             *genotype;
  DataPtr          ptr;
};

struct InterpObject {


    double              *fact_discons;
    int                 fact_discons_size;
    NArrPtr             func;
    NArrPtr             slope;
    int                 maxsize;
    double              maxtime;


};


void FreeSolution(NArrPtr *solution);
void FreeFacts(DataTable *D);
int descend(const int *x, const int *y);
NArrPtr Dat2NArrPtr(DataTable *table, int *maxind);
NArrPtr ConvertAnswer(NArrPtr answer, DArrPtr tabtimes);

void error(const char *format, ...);
void warning(const char *format, ...);
FILE *FindSection(FILE *fp, const char *section);
void KillSection(const char *filename, const char *title);

#endif /* FLYDATA_H_ */
