/*
 * scoring.h
 *
 *  Created on: Feb 8, 2013
 *      Author: zhlou
 */

#ifndef SCORING_H_
#define SCORING_H_

#include <cstdio>
#include <limits>
#include "flyData.h"
#include "DataLists.h"
using namespace std;

class zygotic;
class maternal;
class TheProblem;
class SoDe;


const double FORBIDDEN_MOVE = numeric_limits<double>::max();
const int MAX_PRECISION = 16;
/*** STRUCTURES ************************************************************/

/* range struct for limits and such */

class Range {
public:
  double   lower;
  double   upper;
};

/***************************************************************************
 * The following describes the search space to the scoring function.       *
 * There are two ways to specify limits. Lambda, R, and d always get a     *
 * range for each variable---an upper & lower limit. Elements of T, h      *
 * and m that contribute to u in g(u) can be given a range.                *
 * This is probably the way to go as it definitely results in an           *
 * ergodic search space. However, as I write this code (10/95) we have     *
 * only one set of runs using this method. The alternative, which has      *
 * been mostly used, is to treat T, h and m with a penalty function on     *
 * u of the form:                                                          *
 *                                                                         *
 *           | 0 if exp(Lambda*\sum((T^ij)v_max_j)^2 + h^2 +               *
 *           |                                     (m_i mmax)^2 - 1 < 0    *
 * penalty = |                                                             *
 *           | exp(Lambda*\sum((T^ij)v_max_j)^2 + h^2 +(m_i mmax)^2        *
 *           |                                                otherwise    *
 *                                                                         *
 * where vmax_i is the highest level of gene i in the data, similarly      *
 * for mmax and bcd. This method can be non-ergodic if h and T or m are    *
 * being altered and variables change one at a time. In any case, one      *
 * scheme or the other must be used and the programmer must insure that    *
 * one or another set of pointers (see below) are NULL.                    *
 *                                                                         *
 * NOTE: - Lambda is NOT the same stuff as lambda in equation params or    *
 *         the lambda of the Lam schedule (don't get confused!)            *
 ***************************************************************************/
// Note by Zhihao
// After careful checking with the code, there is no need to make limit's
// to be Range **. 1d array is enough.
class SearchSpace {
public:
    double    *pen_vec;    /* pointer to array that defines penalty function */
    Range     *Rlim;                      /* limits fore promoter strengths */
    Range     *Tlim;        /* limits for T matrix, NULL if pen_vec != NULL */
    Range     *Elim;        /* limits for E matrix, NULL if pen_vec != NULL */
    Range     *mlim;               /* limits for m, NULL if pen_vec != NULL */
    Range     *hlim;               /* limits for h, NULL if pen_vec != NULL */
    Range     *dlim;                 /* limit(s) for diffusion parameter(s) */
    Range     *lambdalim;           /* limits for lambda (prot. decay rate) */
    Range     *taulim;           /* limits for tau (delays) */
    SearchSpace() {
        pen_vec = NULL;
        Rlim = NULL;
        Tlim = NULL;
        Elim = NULL;
        mlim = NULL;
        hlim = NULL;
        dlim = NULL;
        lambdalim = NULL;
        taulim = NULL;
    }
    ~SearchSpace() {
        delete[] taulim;
        delete[] lambdalim;
        delete[] dlim;
        delete[] hlim;
        delete[] mlim;
        delete[] Elim;
        delete[] Tlim;
        delete[] Rlim;
        delete[] pen_vec;
    }
};

/* Yousong's GutInfo struct for writing square diff guts */

typedef struct GutInfo{
  int      flag;                    /* for setting the gut flag in score.c */
  int      ndigits;                                /* gut output precision */
} GutInfo;


class scoring
{
private:
    zygotic &Zygote;
    maternal &TheMaternal;
    const TheProblem &defs;
    int nalleles;
    SoDe *delay_solver;
    int debug;
    GenoType *facttype; /* array of structs that hold facts for each genotype */
    GenoType *tt; /* tt holds the times for which there's  */
                  /* data (one array for each genotype)    */

    SearchSpace *limits; /* structure that holds the limits for   */
                         /* the search space for the annealer (?) */

    GutInfo gutparms; /* Yousong's gutinfo structure used to   */
                      /* print guts below                      */

    InterpObject *polations; /* array of interpobjects for the
                                interpolating functions for history
                                for each genotype */

    InterpObject *extinp_polations; /* array of interpobjects for the
                                       interpolating functions for
                                       external inputs
                                       for each genotype */

    /* ... and a couple of minor things */

    double stepsize; /* static variable for solver stepsize */
    double accuracy; /* static variable for solver accuracy */
    int ndatapoints; /* number of data points (for RMS) */

    const char *filename; /* infile name */
    FILE *slogptr; /* solver log file pointer */

    /* ... plus an init flag */

    // int tt_init_flag = 0; /* flag: InitTT called yet or not? */

    /*** List2Facts: takes a Dlist and returns the corresponding DataTable *****
     *               structure we use for facts data.                          *
     ***************************************************************************/

    DataTable   *List2Facts   (Dlist *inlist);

    /*** ReadLimits: reads the limits section of a data file and returns the  **
     *               approriate SearchSpace struct to the calling function     *
     ***************************************************************************/

    SearchSpace *ReadLimits(FILE *fp);
    void InitHistory(FILE *fp);
    void DoInterp(DataTable *interp_dat, InterpObject *interp_res,
                  int num_genes);
    void GetInterp(FILE *fp, char *title, int num_genes,
                   DataTable **interp_tables);



    /*** Eval: scores the summed squared differences between equation solution *
     *         and data. Because the times for states written to the Solution  *
     *         structure are read out of the data file itself, we do not check *
     *         for consistency of times in this function---all times with data *
     *         will be in the table, but the table may also contain additional *
     *         times.                                                          *
     ***************************************************************************/
    double Eval(NArrPtr Solution, int gindex);

    /*** GutEval: this is the same as Eval, i.e it calculates the summed squa- *
     *            red differences between equation solution and data, with the *
     *            addition that individual squared differences between data-   *
     *            points are written to STDOUT in the unfold output format     *
     ***************************************************************************/

    double GutEval(NArrPtr Solution, int gindex);
public:
    scoring(FILE *fp, zygotic &zy, int flags, int ndigits, double step, double acc, FILE *slog,
            const char *infile, int in_debug);
    ~scoring();

    /*** Score: as the name says, score runs the simulation, gets a solution ***
     *          and then compares it to the data using the Eval least squares  *
     *          function                                                       *
     *   NOTE:  both InitZygote and InitScoring have to be called first!       *
     ***************************************************************************/

    double Score(void);
    void SetGuts (int gutflag, int ndigits)
    {
        gutparms.flag = gutflag;
        gutparms.ndigits = ndigits;
    }
    int GetNDatapoints(void) {return ndatapoints;}

    /*** GetPenalty: calculates penalty from static limits, vmax and mmax ******
     *   CAUTION:    InitPenalty must be called first!                         *
     ***************************************************************************/
    double GetPenalty(void);
    SearchSpace *GetLimits() {return limits;}
    SearchSpace *Penalty2Limits();
};

#endif /* SCORING_H_ */
