/*
 * zygotic.h
 *
 *  Created on: Feb 5, 2013
 *      Author: zhlou
 */

#ifndef ZYGOTIC_H_
#define ZYGOTIC_H_

class maternal;
struct TheProblem;

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

typedef enum GFunc {
  Sqrt,
  Tanh,
  Exp,
  Hvs
} GFunc;



/*** A GLOBAL **************************************************************/


class zygotic
{
private:
    maternal &TheMaternal;
    const TheProblem &defs;

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

    //int rule; /* propagation rule for DvdtOrig */
    int genindex; /* genotype index, needed by DvdtOrig */ // TODO: init me
    /* for getting bicoid from maternal.c */
    EqParms ReadParameters(FILE *fp, char *section_title);

    // old static varibles used in p_deriv
    int dvdt_num_nucs = 0; /* store the number of nucs for next step */
    int dvdt_bcd_index = 0; /* the *next* array in bicoid struct for bcd */
    DArrPtr dvdt_bcd; /* pointer to appropriate bicoid struct */
public:
    zygotic(maternal &in_maternal, FILE *fp, char *parm_section);
    virtual ~zygotic();

    // This is the DvdtOrig
    virtual void p_deriv(double *v, double t, double *vdot, int n);

    // TODO implement me This is the JacobnOrig
    virtual void p_jacobn(double t, double *v, double *dfdt, double **jac, int n);

    // TODO implement me This is the DvdtDelay
    virtual void d_deriv(double *v, double **vd, double t, double *vdot, int n);
    // TODO add mutate

    /*** GetMutParameters: same as above but returns the mutated copy of the ***
     *                     parameter struct; important for writing guts        *
     ***************************************************************************/
    EqParms *GetMutParameters(void) { return &lparm;}

};

#endif /* ZYGOTIC_H_ */
