/*
 * maternal.h
 *
 *  Created on: Feb 4, 2013
 *      Author: zhlou
 */

#ifndef MATERNAL_H_
#define MATERNAL_H_

#include "DataLists.h"
#include "flyData.h"
#include <cfloat>
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


#define TOTAL_DIVS 6

class maternal
{
private:
    static const double old_divtimes[3];

    static const double divtimes1[1];
    static const double divtimes2[2];
    static const double divtimes3[3];
    static const double divtimes4[4];

    /* division durations */

    static const double old_div_duration[3];

    static const double div_duration1[1];
    static const double div_duration2[2];
    static const double div_duration3[3];
    static const double div_duration4[4];

    /* gastrulation times */

    static const double old_gast_time;

    static const double gast_time0;
    static const double gast_time1;
    static const double gast_time2;
    static const double gast_time3;
    static const double gast_time4;

    /* full division times: including t<0 */

    static const double full_divtimes0[TOTAL_DIVS];
    static const double full_divtimes1[TOTAL_DIVS];
    static const double full_divtimes2[TOTAL_DIVS];
    static const double full_divtimes3[TOTAL_DIVS];
    static const double full_divtimes4[TOTAL_DIVS];

    /* division durations */

    static const double full_div_durations[TOTAL_DIVS];

    TheProblem defs;
    Slist *genotypes;
    int nalleles;  /* number of alleles (genotypes) in data */
    double maxconc;     /* max prot conc: 12 (old-), 255 (newstyle) */
    double custom_gast;                  /* custom gastrulation time set by -S */

    /* following is number of nucs in each cleavage cycle, reverse order */
    int *nnucs;
    int *full_nnucs;

    /* following contains lineage numbers at which nuclei at each ccycle start */
    int *lin_start;
    int *full_lin_start;
    int full_ccycles;


    /* The following two store bicoid gradients and bias static to maternal.c  */
    /* bias is found here because it contains the maternal contributions to    */
    /* some zygotic genes (such as hunchback), but it can also be used to add  */
    /* heatshocks during the simulation                                        */
    GenoType *bcdtype;
    GenoType *biastype;

    /* these two static structs are used for storing the times for which       */
    /* there is bias; these times can be retrieved by using GetBTimes          */
    GenoType *bt; /* bias times for each genotype */

    int bt_init_flag; /* flag for BTtable */
    //int d_flag = 0; // we don't need this anymore /* flag for first call to GetD */
    int rule_flag; /* flag for first call to GetRule */
    //int theta_flag = 0; /* flag for first call to Theta*/

    int    olddivstyle;         /* flag: old or new division times? */

    void ReadTheProblem(FILE *fp);
    Slist *ReadGenotypes(FILE *fp);
    void InitBicoid(FILE *fp);
    void InitBias(FILE *fp);
    void InitBTs(void);
    void InitNNucs(void);

    BArrPtr List2Bicoid(Blist *inlist);
    NArrPtr List2Bias(Dlist *inlist);

    double *getD_table; //double *table in old GetD
    void InitGetD(); // initialize getD_table
    void InitTheta();

    double *theta_dt; /* pointer to division time table */
    double *theta_dd; /* pointer to division duration table */

public:
    maternal(FILE *fp, int divstyle);
    ~maternal();
    const TheProblem& getProblem() const {return defs;}

    Blist *ReadBicoid(FILE *fp, char *section); // scoring needs this

    /*** FUNCTIONS THAT RETURN INFO ABOUT THE EMBRYO **************************/

    /*** GetBicoid: returns a bicoid gradients (in form of a DArrPtr) for a ****
     *              specific time and genotype.                                *
     ***************************************************************************/
    DArrPtr GetBicoid(double time, int genindex);

    /*** GetBias: This function returns bias values for a given time and *******
     *            genotype.                                                    *
     ***************************************************************************/
    DArrPtr      GetBias(double time, int genindex);

    /*** GetCCycle: returns cleavage cycle number for a given time *************
     ***************************************************************************/
    unsigned int GetCCycle(double time);

    /*** GetD: returns diffusion parameters D according to the diff. params. ***
     *         in the data file and the diffusion schedule used                *
     *   NOTE: Caller must allocate D_tab                                      *
     ***************************************************************************/
    void GetD(double t, double *d, char diff_schedule, double *D_tab);

    /*** GetRule: returns the appropriate rule for a given time; used by the ***
     *            derivative function                                          *
     ***************************************************************************/
    int          GetRule(double time);

    /*** Theta: Returns the value of theta(t) in the autonomous ***
     *            version of the equations                        *
     **************************************************************/
    int Theta(double time);

    /*** GetNNucs: reads number of nuclei for a given time *********************
     ***************************************************************************/
    int GetNNucs(double t) const;

    /*** GetGastTime: returns time of gastrulation depending on ndivs and ******
     *                olddivstyle; returns 0 in case of an error; if a custom  *
     *                gastrulation time is chosen with -S, it will be returned *
     *                only if it's bigger than the normal gastrulation time    *
     ***************************************************************************/
    double GetGastTime(void) const;

    /*** GetBTimes: returns a sized array of times for which there is bias *****
     ***************************************************************************/
    DArrPtr GetBTimes(char *genotype) const;

    /*** GetDivtable: returns times of cell divisions depending on ndivs and ***
     *                olddivstyle; returns NULL in case of an error            *
     ***************************************************************************/
    double *GetDivtable(void);

    /*** GetDurations: returns pointer to durations of cell divisions de- ******
     *                 pending on ndivs and olddivstyle; returns NULL in case  *
     *                 of an error                                             *
     ***************************************************************************/
    double *GetDurations(void);

    /*** GetStartLin: returns the lineage number of the most anterior nucleus **
     *                for a given time                                         *
     ***************************************************************************/
    int GetStartLin(double t);

    /*** GetStartLinIndex: returns the index of the lineage array for a given **
     *                     time                                                *
     ***************************************************************************/
    int GetStartLinIndex(double t);

    int Index2StartLin(int index) {return full_lin_start[index];}
    int Index2NNuc(int index){ return full_nnucs[index]; }
    int GetNAlleles() const {return nalleles;}
    Slist *GetGenotypes() const {return genotypes;}
    double *Get_Theta_Discons(int *theta_discon_size);
    DataTable *List2Interp(Dlist *inlist, int num_genes);
};


#endif /* MATERNAL_H_ */
