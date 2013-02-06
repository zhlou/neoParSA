/*
 * maternal.h
 *
 *  Created on: Feb 4, 2013
 *      Author: zhlou
 */

#ifndef MATERNAL_H_
#define MATERNAL_H_

#include "DataLists.h"
#include <cstdio>
using namespace std;

/* this is the problem at hand */
struct TheProblem
{
    int ngenes; /* number of genes to consider */
    int egenes; /* number of external inputs */
    char *gene_ids; /* pointer to char with gene IDs */
    char *egene_ids; /* pointer to char with external gene IDs */
    int ndivs; /* number of cell divisions */
    int nnucs; /* number of nuclei at gastrulation */
    char diff_schedule; /* diffusion schedule */
};

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




class maternal
{
private:
    TheProblem defs;
    Slist *genotypes;
    int nalleles;  /* number of alleles (genotypes) in data */
    double maxconc;     /* max prot conc: 12 (old-), 255 (newstyle) */

    /* following is number of nucs in each cleavage cycle, reverse order */
    int *nnucs;
    //int *full_nnucs = NULL;

    /* following contains lineage numbers at which nuclei at each ccycle start */
    int *lin_start;
    //int *full_lin_start = NULL;
    //int full_ccycles;


    /* The following two store bicoid gradients and bias static to maternal.c  */
    /* bias is found here because it contains the maternal contributions to    */
    /* some zygotic genes (such as hunchback), but it can also be used to add  */
    /* heatshocks during the simulation                                        */
    GenoType *bcdtype;
    GenoType *biastype;

    /* these two static structs are used for storing the times for which       */
    /* there is bias; these times can be retrieved by using GetBTimes          */
    GenoType *bt; /* bias times for each genotype */

    int bt_init_flag = 0; /* flag for BTtable */
    //int d_flag = 0; /* flag for first call to GetD */
    //int rule_flag = 0; /* flag for first call to GetRule */
    //int theta_flag = 0; /* flag for first call to Theta*/


    void ReadTheProblem(FILE *fp);
    Slist *ReadGenotypes(FILE *fp);
    void InitBicoid(FILE *fp);
    void InitBias(FILE *fp);
    void InitBTs(void);
    void InitNNucs(void);
    Blist *ReadBicoid(FILE *fp, char *section);
    BArrPtr List2Bicoid(Blist *inlist);
    NArrPtr List2Bias(Dlist *inlist);

public:
    maternal(FILE *fp);
    ~maternal();
    const TheProblem& getProblem() const {return defs;}
};

#endif /* MATERNAL_H_ */
