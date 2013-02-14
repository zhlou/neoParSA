/*
 * scoring.h
 *
 *  Created on: Feb 8, 2013
 *      Author: zhlou
 */

#ifndef SCORING_H_
#define SCORING_H_

#include <cstdio>
using namespace std;

/*** STRUCTURES ************************************************************/

/* range struct for limits and such */

typedef struct Range {
  double   lower;
  double   upper;
} Range;

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

typedef struct {
  double    *pen_vec;    /* pointer to array that defines penalty function */
  Range     **Rlim;                      /* limits fore promoter strengths */
  Range     **Tlim;        /* limits for T matrix, NULL if pen_vec != NULL */
  Range     **Elim;        /* limits for E matrix, NULL if pen_vec != NULL */
  Range     **mlim;               /* limits for m, NULL if pen_vec != NULL */
  Range     **hlim;               /* limits for h, NULL if pen_vec != NULL */
  Range     **dlim;                 /* limit(s) for diffusion parameter(s) */
  Range     **lambdalim;           /* limits for lambda (prot. decay rate) */
  Range     **taulim;           /* limits for tau (delays) */
} SearchSpace;

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
    TheProblem &defs;
    SoDe *delay_solver;
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
    int ndatapoints = 0; /* number of data points (for RMS) */

    char *filename; /* infile name */
    FILE *slogptr = NULL; /* solver log file pointer */

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
                   DataTable **interp_tables)
public:
    scoring(FILE *fp, zygotic &zy, double step, double acc, FILE *slog, char *infile);
    ~scoring();
};

#endif /* SCORING_H_ */
