/*
 * zygotic.h
 *
 *  Created on: Feb 5, 2013
 *      Author: zhlou
 */

#ifndef ZYGOTIC_H_
#define ZYGOTIC_H_
#include <cstdio>
#include <typeinfo>
#include "solvers.h"
#include "maternal.h"
using namespace std;


struct EqParms
{
    double *R; /* strength of each promoter--always >= 0. */
    double *T; /* the genetic interconnect matrix */
    double *E; /* the external input regulatory matrix */
    double *m; /* regulatory coefficient of bcd on gene */
    double *h; /* reg. coeff. for generic TFs on synthesis of gene */
    double *d; /* spatial interaction at gastrulation--always >= 0. */
    double *lambda; /*protein half lives--always >= 0. */
    double *tau; /* delay times for the proteins */
};

/*** AN ENUM ***************************************************************/

/* This is the g(u)-function enum which describes the different types of   *
 * g(u) functions we can use in derivative functions                       */

enum  GFunc{
  Sqrt,
  Tanh,
  Exp,
  Hvs
};



/*** A GLOBAL **************************************************************/


class zygotic
{
private:
    maternal &TheMaternal;
    const TheProblem &defs;
    int debug;

    GFunc    gofu;                            /* the g(u) function we're using */

    /* following two structs are used for equation paramters; parm holds       */
    /* parameters BEFORE mutation and is returned by GetParameters (i.e. to    */
    /* Score(), which would produce an error on mutated parameters); lparm is  */
    /* a copy of the struct that gets mutated (either Rs or Ts get set to zero */
    /* and therefore violate limits if sent to Score()) and is then used by    */
    /* the derivative function DvdtOrig                                        */

    EqParms parm; /* static struct for equation parameters */
    EqParms lparm; /* local copy of parameter that gets mutated */

    /* following arrays are used by DvdtOrig */

    double *D; /* contains info about diffusion sched. */
    double *vinput; /* vinput, bot2 and bot are used for */
    double *bot2, *bot; /* storing intermediate stuff for vector */
    /* functions used by the dvdt function */

    /* two local copies for propagation rule and genotype number */

    int rule; /* propagation rule for DvdtOrig */
    int genindex; /* genotype index, needed by DvdtOrig */ // it is set in blastoderm
    /* for getting bicoid from maternal.c */
    EqParms ReadParameters(FILE *fp, const char *section_title);

    // old static varibles used in p_deriv
    int dvdt_num_nucs; /* store the number of nucs for next step */
    int dvdt_bcd_index; /* the *next* array in bicoid struct for bcd */
    DArrPtr dvdt_bcd; /* pointer to appropriate bicoid struct */

    /* A STRUCT ****************************************************************/

    struct TList {          /* Tlist is a linked list we use to set up */
      double       time;    /* the solution struct in integrate.c; the */
      int          n;       /* TList has an element for each time con- */
      int          op;      /* taining the number of elements for the  */
      struct TList *next;   /* solution struct at that time (n) and a  */
    };                      /* rule (op) to tell the solver what to do */

    /* MORE CONSTANTS: OPS FOR BLASTODERM -- WITH PRIORITIES *******************/

    static const int ADD_BIAS;  /* 1, can happen any step                     */
    static const int NO_OP;     /* 2, ops below: first found set executed     */
    static const int DIVIDE;    /* 4, NO_OP does nothing, DIVIDE divides nucs */
    static const int PROPAGATE; /* 8, and PROPAGATE propagates the equations  */
    static const int MITOTATE;  /* 16, and MITOTATE carries forward the       *
                                 * system by epsilon (to handle rounding      *
                                 * errors from the solver                     */

    /* TList Utility Functions */

    /*** InitTList: initializes TList and adds first (t=0) and last ************
     *              (t=gastrulation) element of Tlist.                         *
     ***************************************************************************/

    TList *InitTList(void);

    /*** InsertTList: takes pointer to first element of TList plus time and ****
     *                desired op for new TList element and inserts a new TList *
     *                element for time at the appropriate place within the     *
     *                linked list. This function returns a pointer to the      *
     *                first element of the TList if insertion worked out fine. *
     ***************************************************************************/

    TList *InsertTList(TList *first, double time, int op);

    /*** CountEntries: counts how many entries we have in a TList **************
     ***************************************************************************/

    int CountEntries(TList *first);

    /*** FreeTList: frees the memory of a TList ********************************
     ***************************************************************************/

    void FreeTList(TList *first);


    solver *solve;

    /* Mutator functions */

    /*** Mutate: calls mutator functions according to genotype string **********
     ***************************************************************************/

    void Mutate(char *g_type);

    /*** T_Mutate: mutates genes by setting all their T matrix entries *********
     *             to zero. Used to simulate mutants that express a            *
     *             non-functional protein.                                     *
     ***************************************************************************/

    void T_Mutate(int g_type);

    /*** R_Mutate: mutates genes by setting their promotor strength R **********
     *             to zero, so there will be no transcription at all           *
     *             anymore. Used to simulate mutants that don't pro-           *
     *             duce any protein anymore.                                   *
     ***************************************************************************/

    void R_Mutate(int g_type);

    /*** RT_Mutate: mutates gene by setting both promoter strength R ***********
     *              and T matrix entries to zero. Useful, if you want          *
     *              to be really sure that there is NO protein pro-            *
     *              duction and that maternal contribution also don't          *
     *              contribute anything to gene interactions.                  *
     ***************************************************************************/

    void RT_Mutate(int g_type);

    /*** CopyParm: copies all the parameters into the lparm struct *************
     ***************************************************************************/

    EqParms CopyParm(EqParms orig_parm);

    /*** SetRule: sets the static variable rule to MITOSIS or INTERPHASE *******
     ***************************************************************************/

    void SetRule(int r) { rule = r; }

    /*** FreeMutant: frees mutated parameter struct ****************************
     ***************************************************************************/

    void FreeMutant(void);



public:

    zygotic(maternal &in_maternal, FILE *fp, const char *parm_section,
            GFunc in_gofu, int debug, const char *solver_name);
    virtual ~zygotic();

    // This is the DvdtOrig
    virtual void p_deriv(double *v, double t, double *vdot, int n);

    virtual void p_jacobn(double t, double *v, double *dfdt, double **jac, int n);

    // This is the DvdtDelay
    virtual void d_deriv(double *v, double **vd, double t, double *vdot, int n);


    /*** GetMutParameters: same as above but returns the mutated copy of the ***
     *                     parameter struct; important for writing guts        *
     ***************************************************************************/
    EqParms *GetMutParameters(void) { return &lparm;}
    /*** Blastoderm: runs embryo model and returns an array of concentration ***
     *               arrays for each requested time (given by TabTimes) using  *
     *               stephint as a suggested stepsize and accuracy as the re-  *
     *               quired numerical accuracy (in case of adaptive stepsize   *
     *               solvers or global stepsize control) for the solver.       *
     *         NOTE: TabTimes *must* start from 0 and have increasing times.   *
     *               It includes times when bias is added, cell division times *
     *               and times for which we have data and ends with the time   *
     *               of gastrulation.                                          *
     ***************************************************************************/

    NArrPtr Blastoderm(int genindex, char *genotype,
                       InterpObject hist_interrp,
                       InterpObject extinp_interrp, DArrPtr tabtimes,
                       double stephint, double accuracy, FILE *slog);

    void PrintBlastoderm(FILE *fp, NArrPtr table, char *id,
                 int ndigits, int columns);

    maternal &get_Maternal() const {return TheMaternal;}
    int get_ngenes() const {return defs.ngenes;}
    int get_egenes() const {return defs.egenes;}
    int get_NNucs(double t) const {return TheMaternal.GetNNucs(t);}
    EqParms *GetParameters(void) {return &parm;}
    SoDe *get_delay_solver () const
    {
        return (typeid(*solve) == typeid(SoDe)) ? dynamic_cast<SoDe *>(solve) :
                NULL;
    }

    void PrintParameters(FILE *fp, EqParms *p, const char *title, int ndigits);

};

#endif /* ZYGOTIC_H_ */
