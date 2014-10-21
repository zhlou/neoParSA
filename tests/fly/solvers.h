/*
 * solvers.h
 *
 *  Created on: Feb 6, 2013
 *      Author: zhlou
 */

#ifndef SOLVERS_H_
#define SOLVERS_H_

#include <cstdio>
#include "flyData.h"
using namespace std;

class zygotic;

class solver
{
public:
    solver(zygotic &in_zy, int in_debug) : zygote(in_zy), debug(in_debug) {};
    virtual ~solver() {};
    virtual void ps(double *, double *, double, double, double, double, int,
            FILE *) = 0;
protected:
    void WriteSolvLog(const char *solver, double tin, double tout, double h, int n,
              int nderivs, FILE *slog);
    int debug;
    zygotic &zygote;
};

class Euler : public solver
{
public:
    // for newer compilers support c++11, we can simply use
    // using solver::solver;
    Euler(zygotic &in_zy, int in_debug) : solver(in_zy, in_debug) {};
    void ps(double *vin, double *vout, double tin, double tout,
            double stephint, double accuracy, int n, FILE *slog);
};

/*** Meuler: propagates vin (of size n) from tin to tout by the Modified ***
 *           Euler method (this is NOT the midpoint method, see Rk2());    *
 *           the result is returned by vout                                *
 ***************************************************************************/
class Meuler : public solver
{
public:
    Meuler(zygotic &in_zy, int in_debug) : solver(in_zy, in_debug) {};
    void ps(double *vin, double *vout, double tin, double tout,
            double stephint, double accuracy, int n, FILE *slog);
};

/*** Heun: propagates vin (of size n) from tin to tout by Heun's method ****
 *         the result is returned by vout                                  *
 ***************************************************************************/
class Heun : public solver
{
public:
    Heun(zygotic &in_zy, int in_debug) : solver(in_zy, in_debug) {};
    void ps(double *vin, double *vout, double tin, double tout,
            double stephint, double accuracy, int n, FILE *slog);
};

/*** Rk2: propagates vin (of size n) from tin to tout by the Midpoint or ***
 *        Second-Order Runge-Kutta method; the result is returned by vout  *
 ***************************************************************************/
class Rk2 : public solver
{
public:
    Rk2(zygotic &in_zy, int in_debug) : solver(in_zy, in_debug) {};
    void ps(double *vin, double *vout, double tin, double tout,
            double stephint, double accuracy, int n, FILE *slog);
};

/*** Rk4: propagates vin (of size n) from tin to tout by the Fourth-Order **
 *        Runge-Kutta method; the result is returned by vout               *
 ***************************************************************************
 *                                                                         *
 * written by Joel Linton (somewhere around 1998)                          *
 * fixed and modified by Yoginho (somewhere around 2001)                   *
 *                                                                         *
 ***************************************************************************/
class Rk4 : public solver
{
public:
    Rk4(zygotic &in_zy, int in_debug) : solver(in_zy, in_debug) {};
    void ps(double *vin, double *vout, double tin, double tout,
            double stephint, double accuracy, int n, FILE *slog);
};

/*** Rkck: propagates vin (of size n) from tin to tout by the Runge-Kutta **
 *         Cash-Karp method, which is an adaptive-stepsize Rk method; it   *
 *         uses a fifth-order Rk formula with an embedded forth-oder for-  *
 *         mula for calucalting the error; its result is returned by vout  *
 ***************************************************************************
 *                                                                         *
 * This solver was written by Marcel Wolf, Spring 2002.                    *
 *                                                                         *
 ***************************************************************************/
class Rkck : public solver
{
public:
    Rkck(zygotic &in_zy, int in_debug) : solver(in_zy, in_debug) {};
    void ps(double *vin, double *vout, double tin, double tout,
            double stephint, double accuracy, int n, FILE *slog);
};

/*** Rkf: propagates vin (of size n) from tin to tout by the Runge-Kutta ***
 *        Fehlberg method, which is a the original adaptive-stepsize Rk    *
 *        method (Cash-Karp is an improved version of this); it uses a     *
 *        fifth-order Rk formula with an embedded forth-oder formula for   *
 *        calucalting the error; its result is returned by vout            *
 ***************************************************************************
 *                                                                         *
 * This solver was written by Marcel Wolf, Spring 2002.                    *
 *                                                                         *
 ***************************************************************************/
class Rkf : public solver
{
public:
    Rkf(zygotic &in_zy, int in_debug) : solver(in_zy, in_debug) {};
    void ps(double *vin, double *vout, double tin, double tout,
            double stephint, double accuracy, int n, FILE *slog);
};

/***** Milne: propagates vin (of size n) from tin to tout by Milne-Simpson *
 *            which is a predictor-corrector method; the result is retur-  *
 *            ned by vout                                                  *
 ***************************************************************************
 *                                                                         *
 * This solver was implemented by Konstantin Koslov, Dec 2001/Jan 2002     *
 *                                                                         *
 ***************************************************************************/
class Milne : public solver
{
public:
    Milne(zygotic &in_zy, int in_debug) : solver(in_zy, in_debug) {};
    void ps(double *vin, double *vout, double tin, double tout,
            double stephint, double accuracy, int n, FILE *slog);
};

/***** Adams: propagates vin (of size n) from tin to tout by Adams-Moulton *
 *            which is an implicit predictor-corrector method of second    *
 *            order; the result is returned by vout                        *
 ***************************************************************************
 *                                                                         *
 * This solver was implemented by Konstantin Koslov, Spring 2002           *
 * Slightly modified by Manu, July 2002                                    *
 *                                                                         *
 ***************************************************************************/
class Adams : public solver
{
public:
    Adams(zygotic &in_zy, int in_debug) : solver(in_zy, in_debug) {};
    void ps(double *vin, double *vout, double tin, double tout,
            double stephint, double accuracy, int n, FILE *slog);
};

class StepSolver : public solver
{
public:
    StepSolver(zygotic &in_zy, int in_debug) : solver(in_zy, in_debug) { d = hpoints = NULL;};
protected:
    /*          pzextr: implements the Richardson extrapolation (polynomial)   */
    void pzextr(int iest, double hest, double *yest, double *yz, double *dy,
            int nv, const int KMAXX);
    double *d;
    double *hpoints;
};
/***** BuSt: propagates v(t) from t1 to t2 by Bulirsch-Stoer; this method **
 *           uses Richardson extrapolation to estimate v's at a hypothe-   *
 *           tical stepsize of 0; the extrapolation also yields an error   *
 *           estimate, which is used to adapt stepsize and change the or-  *
 *           der of the method as required; the result is returned by vout *
 ***************************************************************************
 *                                                                         *
 * This solver was implemented by Manu, July 2002                          *
 *                                                                         *
 ***************************************************************************/
class BuSt : public StepSolver
{
public:
    BuSt(zygotic &in_zy, int in_debug) : StepSolver(in_zy, in_debug) {};
    void ps(double *vin, double *vout, double tin, double tout,
            double stephint, double accuracy, int n, FILE *slog);
private:
    /* bsstep is the BuSt stepper function */
    void bsstep(double *v, double *deriv, int n, double *t, double htry,
            double accuracy, double *hdid, double *hnext);
    /* func prototypes: these funcs should not be visible outside solvers.c    *
     *            mmid: implements the modified midpoint method                */


    void mmid(double *vin, double *vout, double *deriv, double tin, double htot,
            int nstep, int n);


};

/***** BaDe: propagates v(t) from t1 to t2 by Bader-Deuflhard; this method *
 *           uses Richardson extrapolation to estimate v's at a hypothe-   *
 *           tical stepsize of 0; the extrapolation also yields an error   *
 *           estimate, which is used to adapt stepsize and change the or-  *
 *           der of the method as required; the result is returned by vout *
 ***************************************************************************
 *                                                                         *
 * This solver was implemented by Yogi, based on BuSt, Aug 2002            *
 *                                                                         *
 ***************************************************************************/
class BaDe : public StepSolver
{
public:
    BaDe(zygotic &in_zy, int in_debug) : StepSolver(in_zy, in_debug) {};
    void ps(double *vin, double *vout, double tin, double tout,
            double stephint, double accuracy, int n, FILE *slog);
private:
    /* stifbs is the BaDe stepper function */

    void stifbs(double *v, double *deriv, int n, double *t, double htry,
            double accuracy, double *hdid, double *hnext);
    /* func prototypes: these funcs should not be visible outside solvers.c    *
     *           simpr: implements the semi-implicit mid-point rule            *
     *          pzextr: implements the Richardson extrapolation (polynomial)   */

    void simpr(double *vin, double *vout, double *deriv, double *dfdt,
            double **jac, double tin, double htot, int nstep, int n);
    /* func prototypes: these funcs should not be visible outside solvers.c    *
     *          lubcmp: does DU decomposition of a matrix                      *
     *          lubksb: solves linear system a * b                             */

      void ludcmp(double **a, int n, int *indx, double *d);
      void lubksb(double **a, int n, int *indx, double *b);

};


int compare(const double *x, const double *y);
// forward declaration
struct EqParms;
class maternal;
class TheProblem;


// The delay solver
class SoDe : public solver
{
public:
    SoDe(zygotic &in_zy, int in_debug);
    ~SoDe();
    void ps(double *vin, double *vout, double tin, double tout,
            double stephint, double accuracy, int n, FILE *slog);
    void resetSolver();
    void SetHistoryInterp(InterpObject interp_info);
    void SetExternalInputInterp(InterpObject interp_info);
    void Go_Forward(double *output, double *input, int output_ind,
                    int input_ind, int num_genes);
    void Go_Backward(double *output, double *input, int output_ind,
                     int input_ind, int num_genes);
    void ExternalInputs(double t, double t_size, double *yd, int n);
    void DivideHistory(double t1, double t2);
private:
    InterpObject hist_interp_object, extinp_interp_object;
    EqParms *lp;
    double maxdel, mindel;
    int    numdel;  /* delay parameters used by DCERk32, y_delayed */
    double *delay;      /* static array set in SoDe, used by DCERk32 */
    int gridstart;  /* for the heuristic */
    int gridpos;                 /* where you are in the grid */
    double *tdone;              /* the grid */
    double **derivv1;              /* intermediate derivatives
                                                            for the Cash
                                                            -Karp formula */
    double **derivv2;
    double **derivv3;
    double **derivv4;
    double **vdonne;

    double *fact_discons, fact_discons_size;
    maternal &TheMaternal;
    const TheProblem &defs;
    /*** DCERk3(2): propagates v[0] (of size n) according to tarray  by  **
    the Runge-Kutta 3(2) pair with continuous extension storing the      **
    result in vatt. Initial conditions are specified in vatt[0],         **
    corresponding to tarray[0]                                           **
    ***********************************************************************/
    void DCERk32(double **vatt, int n, double *tarray, int tpoints, double
                *darray, int dpoints, double stephint, double accuracy);
    int y_delayed(double ***vd, int n, double *rktimes, double *tau,
                  double *grid, double **vdone, double **deriv1,
                  double **deriv2, double **deriv3, double **deriv4,
                  int gridsize, double accu);
    double *Construct_Discont_Array(double range, double *taus, int n,
                                    double *starts, int sn, int *disc_size);
    void CE(double t, double *vans, double tbegin, double *v_at_tbegin,
            double ech, double *d1, double *d2, double *d3, double *d4, int n);
    void History(double t, double t_size, double *yd, int n);



};

solver *SolverFactory(zygotic &zygote, int debug,
                      const char *name="Rk4");
#endif /* SOLVERS_H_ */
