/*
 * solvers.cpp
 *
 *  Created on: Feb 6, 2013
 *      Author: zhlou
 */

#include "solvers.h"
#include "zygotic.h"



void Euler::ps(double* vin, double* vout, double tin, double tout,
        double stephint, double accuracy, int n, FILE* slog)
{
    int    i;                                          /* local loop counter */

    double **v;                 /* intermediate v's, used to toggle v arrays */
    int    toggle = 0;               /* used to toggle between v[0] and v[1] */

    double *vnow;                                /* ptr to v at current time */
    double *vnext;                    /* ptr to v at current time + stepsize */

    double *deriv;             /* the derivatives at the beginning of a step */

    double deltat;                                             /* tout - tin */
    double t;                                                /* current time */

    double m;         /* tmp var to calculate precise stepsize from stephint */
    double stepsize;                                        /* real stepsize */
    int    step;                                   /* loop counter for steps */
    int    nsteps;                  /* total number of steps we have to take */

    int    nd = 0;                            /* number of deriv evaluations */

  /* steps too small: just return */

    if (tin == tout)
      return;

  /* if steps big enough */

    v = (double **)calloc(2, sizeof(double *));

    v[0]  = (double *)calloc(n, sizeof(double));
    v[1]  = (double *)calloc(n, sizeof(double));
    deriv = (double *)calloc(n, sizeof(double));

    deltat = tout - tin;                 /* how far do we have to propagate? */
    m = floor(deltat/stephint + 0.5);       /* see comment on stephint above */
    if ( m < 1. )
      m = 1.;                          /* we'll have to do at least one step */

    stepsize = deltat / m;                  /* real stepsize calculated here */
    nsteps = (int)m;                                  /* int number of steps */
    t = tin;                                             /* set current time */

    if ( t == t + stepsize )
      error("Euler: stephint of %g too small!", stephint);

    vnow = vin;

    if (nsteps == 1)       /* vnext holds the results after the current step */
      vnext = vout;
    else
      vnext = v[0];

    for (step=0; step < nsteps; step++) {                  /* loop for steps */

      zygote.p_deriv(vnow, t, deriv, n);  /* call derivative func to evaluate deriv */

      if ( debug ) nd++;

      for(i=0; i < n; i++)
        vnext[i] = vnow[i] + stepsize * deriv[i];           /* Euler formula */

      t += stepsize;                                      /* go to next step */

      if (step < nsteps - 2) {                   /* CASE 1: many steps to go */
        vnow = v[toggle];               /* toggle results from v[0] and v[1] */
        toggle++;
        toggle %= 2;
        vnext = v[toggle];

      } else if (step == nsteps - 2) {     /* CASE 2: next iteration = final */
        vnow = v[toggle];                                        /* set vout */
        vnext = vout;

      } else if (step > nsteps - 2) {    /* CASE 3: just did final iteration */
        free(v[0]);                                 /* clean up and go home! */
        free(v[1]);
        free(v);
        free(deriv);
      }
    }

    if ( debug )
      WriteSolvLog("Euler", tin, tout, stepsize, nsteps, nd, slog);

    return;

}

void solver::WriteSolvLog(char* solver, double tin, double tout, double h,
        int n, int nderivs, FILE* slog)
{
    double nds;                 /* Number of Derivative evaluations per Step */
    double ttot;                       /* total time of propagation interval */

    nds = (double)nderivs/(double)n;
    ttot = tout - tin;

    fprintf(slog, "%s:   tin = %7.3f   tout = %7.3f   ttot = %7.3f   ",
        solver, tin, tout, ttot);
    fprintf(slog, "h = %5.3f   nsteps = %4d    nderivs/step = %4.1f\n",
        h, n, nds);
}
