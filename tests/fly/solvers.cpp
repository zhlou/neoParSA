/*
 * solvers.cpp
 *
 *  Created on: Feb 6, 2013
 *      Author: zhlou
 */

#include "solvers.h"
#include "zygotic.h"
#include <cfloat>
#include <algorithm>
#include <cmath>
#include <cstring>
using namespace std;

// always use inline function instead of macro
inline double DSQR(double a) {
    return (a == 0.0) ? 0.0 : (a * a);
}

void solver::WriteSolvLog(char* solver, double tin, double tout, double h,
                          int n, int nderivs, FILE* slog)
{
    double nds; /* Number of Derivative evaluations per Step */
    double ttot; /* total time of propagation interval */

    nds = (double) nderivs / (double) n;
    ttot = tout - tin;

    fprintf(slog, "%s:   tin = %7.3f   tout = %7.3f   ttot = %7.3f   ", solver,
            tin, tout, ttot);
    fprintf(slog, "h = %5.3f   nsteps = %4d    nderivs/step = %4.1f\n", h, n,
            nds);
}

void Euler::ps(double* vin, double* vout, double tin, double tout,
               double stephint, double accuracy, int n, FILE* slog)
{
    int i; /* local loop counter */

    double **v; /* intermediate v's, used to toggle v arrays */
    int toggle = 0; /* used to toggle between v[0] and v[1] */

    double *vnow; /* ptr to v at current time */
    double *vnext; /* ptr to v at current time + stepsize */

    double *deriv; /* the derivatives at the beginning of a step */

    double deltat; /* tout - tin */
    double t; /* current time */

    double m; /* tmp var to calculate precise stepsize from stephint */
    double stepsize; /* real stepsize */
    int step; /* loop counter for steps */
    int nsteps; /* total number of steps we have to take */

    int nd = 0; /* number of deriv evaluations */

    /* steps too small: just return */

    if (tin == tout)
        return;

    /* if steps big enough */

    v = (double **) calloc(2, sizeof(double *));

    v[0] = (double *) calloc(n, sizeof(double));
    v[1] = (double *) calloc(n, sizeof(double));
    deriv = (double *) calloc(n, sizeof(double));

    deltat = tout - tin; /* how far do we have to propagate? */
    m = floor(deltat / stephint + 0.5); /* see comment on stephint above */
    if (m < 1.)
        m = 1.; /* we'll have to do at least one step */

    stepsize = deltat / m; /* real stepsize calculated here */
    nsteps = (int) m; /* int number of steps */
    t = tin; /* set current time */

    if (t == t + stepsize)
        error("Euler: stephint of %g too small!", stephint);

    vnow = vin;

    if (nsteps == 1) /* vnext holds the results after the current step */
        vnext = vout;
    else
        vnext = v[0];

    for (step = 0; step < nsteps; step++) { /* loop for steps */

        zygote.p_deriv(vnow, t, deriv, n); /* call derivative func to evaluate deriv */

        if (debug)
            nd++;

        for (i = 0; i < n; i++)
            vnext[i] = vnow[i] + stepsize * deriv[i]; /* Euler formula */

        t += stepsize; /* go to next step */

        if (step < nsteps - 2) { /* CASE 1: many steps to go */
            vnow = v[toggle]; /* toggle results from v[0] and v[1] */
            toggle++;
            toggle %= 2;
            vnext = v[toggle];

        } else if (step == nsteps - 2) { /* CASE 2: next iteration = final */
            vnow = v[toggle]; /* set vout */
            vnext = vout;

        } else if (step > nsteps - 2) { /* CASE 3: just did final iteration */
            free(v[0]); /* clean up and go home! */
            free(v[1]);
            free(v);
            free(deriv);
        }
    }

    if (debug)
        WriteSolvLog("Euler", tin, tout, stepsize, nsteps, nd, slog);

    return;

}

/*** Meuler: propagates vin (of size n) from tin to tout by the Modified ***
 *           Euler method (this is NOT the midpoint method, see Rk2());    *
 *           the result is returned by vout                                *
 ***************************************************************************/
void Meuler::ps(double* vin, double* vout, double tin, double tout,
                double stephint, double accuracy, int n, FILE* slog)
{
    int i; /* local loop counter */

    double **v; /* intermediate v's, used to toggle v arrays */
    int toggle = 0; /* used to toggle between v[0] and v[1] */

    double *vtemp; /* guessed intermediate v's for midpoint */

    double *vnow; /* ptr to v at current time */
    double *vnext; /* ptr to v at current time + stepsize */

    double *deriv1; /* the derivatives at the beginning of a step */
    double *deriv2; /* the derivatives at midpoint of a step */

    double deltat; /* tout - tin */
    double t; /* current time */
    double th; /* time at end of a step (t+stepsize) */

    double m; /* tmp var to calculate precise stepsize from stephint */
    double stepsize; /* real stepsize */
    int step; /* loop counter for steps */
    int nsteps; /* total number of steps we have to take */

    int nd = 0; /* number of deriv evaluations */

    /* the do-nothing case; too small steps dealt with under usual */

    if (tin == tout)
        return;

    /* the usual case: steps big enough */

    v = (double **) calloc(2, sizeof(double *));

    vtemp = (double *) calloc(n, sizeof(double));
    v[0] = (double *) calloc(n, sizeof(double));
    v[1] = (double *) calloc(n, sizeof(double));
    deriv1 = (double *) calloc(n, sizeof(double));
    deriv2 = (double *) calloc(n, sizeof(double));

    deltat = tout - tin; /* how far do we have to propagate? */
    m = floor(deltat / stephint + 0.5); /* see comment on stephint above */
    if (m < 1.)
        m = 1.; /* we'll have to do at least one step */

    stepsize = deltat / m; /* real stepsize calculated here */
    nsteps = (int) m; /* int number of steps */
    t = tin; /* set current time */

    if (t == t + stepsize)
        error("Meuler: stephint of %g too small!", stephint);

    vnow = vin;

    if (nsteps == 1) /* vnext holds the results after the current step */
        vnext = vout;
    else
        vnext = v[0];

    for (step = 0; step < nsteps; step++) { /* loop for steps */

        th = t + stepsize; /* time of step endpoint */

        zygote.p_deriv(vnow, t, deriv1, n); /* first call to deriv */

        if (debug)
            nd++;

        for (i = 0; i < n; i++) /* evaluate v's at endpoint */
            vtemp[i] = vnow[i] + stepsize * deriv1[i];

        zygote.p_deriv(vtemp, th, deriv2, n); /* second call to deriv */

        if (debug)
            nd++;

        for (i = 0; i < n; i++) /* Modified Euler formula */
            vnext[i] = vnow[i] + (stepsize / 2.) * (deriv1[i] + deriv2[i]);

        t += stepsize; /* go to next step */

        if (step < nsteps - 2) { /* CASE 1: many steps to go */
            vnow = v[toggle]; /* toggle results from v[0] and v[1] */
            toggle++;
            toggle %= 2;
            vnext = v[toggle];

        } else if (step == nsteps - 2) { /* CASE 2: next iteration = final */
            vnow = v[toggle]; /* set vout */
            vnext = vout;

        } else if (step > nsteps - 2) { /* CASE 3: just did final iteration */
            free(v[0]); /* clean up and go home! */
            free(v[1]);
            free(v);
            free(vtemp);
            free(deriv1);
            free(deriv2);
        }
    }

    if (debug)
        WriteSolvLog("Meuler", tin, tout, stepsize, nsteps, nd, slog);

    return;
}

/*** Heun: propagates vin (of size n) from tin to tout by Heun's method ****
 *         the result is returned by vout                                  *
 ***************************************************************************/

void Heun::ps(double *vin, double *vout, double tin, double tout,
              double stephint, double accuracy, int n, FILE *slog)
{
    int i; /* local loop counter */

    double **v; /* intermediate v's, used to toggle v arrays */
    int toggle = 0; /* used to toggle between v[0] and v[1] */

    double *vtemp; /* guessed intermediate v's for midpoint */

    double *vnow; /* ptr to v at current time */
    double *vnext; /* ptr to v at current time + stepsize */

    double *deriv1; /* the derivatives at the beginning of a step */
    double *deriv2; /* the derivatives at midpoint of a step */

    double deltat; /* tout - tin */
    double t; /* current time */
    double thh; /* time at 2/3 of the step */

    double m; /* tmp var to calculate precise stepsize from stephint */
    double stepsize; /* real stepsize */
    double hh; /* 2/3 the stepsize */
    int step; /* loop counter for steps */
    int nsteps; /* total number of steps we have to take */

    int nd = 0; /* number of deriv evaluations */

    /* the do-nothing case; too small steps dealt with under usual */

    if (tin == tout)
        return;

    /* the usual case: steps big enough */

    v = (double **) calloc(2, sizeof(double *));

    vtemp = (double *) calloc(n, sizeof(double));
    v[0] = (double *) calloc(n, sizeof(double));
    v[1] = (double *) calloc(n, sizeof(double));
    deriv1 = (double *) calloc(n, sizeof(double));
    deriv2 = (double *) calloc(n, sizeof(double));

    deltat = tout - tin; /* how far do we have to propagate? */
    m = floor(deltat / stephint + 0.5); /* see comment on stephint above */
    if (m < 1.)
        m = 1.; /* we'll have to do at least one step */

    stepsize = deltat / m; /* real stepsize calculated here */
    nsteps = (int) m; /* int number of steps */
    t = tin; /* set current time */

    if (t == t + stepsize)
        error("Heun: stephint of %g too small!", stephint);

    vnow = vin;

    if (nsteps == 1) /* vnext holds the results after the current step */
        vnext = vout;
    else
        vnext = v[0];

    hh = stepsize * 2. / 3.; /* evaluate 2/3 of stepsize */

    for (step = 0; step < nsteps; step++) { /* loop for steps */

        thh = t + hh; /* time at 2/3 of the step */

        zygote.p_deriv(vnow, t, deriv1, n); /* first call to deriv */

        if (debug)
            nd++;

        for (i = 0; i < n; i++) /* evaluate v's at 2/3 of the step */
            vtemp[i] = vnow[i] + hh * deriv1[i];

        zygote.p_deriv(vtemp, thh, deriv2, n); /* second call to deriv */

        if (debug)
            nd++;

        for (i = 0; i < n; i++) /* Heun's formula */
            vnext[i] = vnow[i] + (stepsize / 4.) * (deriv1[i] + 3. * deriv2[i]);

        t += stepsize; /* go to next step */

        if (step < nsteps - 2) { /* CASE 1: many steps to go */
            vnow = v[toggle]; /* toggle results from v[0] and v[1] */
            toggle++;
            toggle %= 2;
            vnext = v[toggle];

        } else if (step == nsteps - 2) { /* CASE 2: next iteration = final */
            vnow = v[toggle]; /* set vout */
            vnext = vout;

        } else if (step > nsteps - 2) { /* CASE 3: just did final iteration */
            free(v[0]); /* clean up and go home! */
            free(v[1]);
            free(v);
            free(vtemp);
            free(deriv1);
            free(deriv2);
        }
    }

    if (debug)
        WriteSolvLog("Heun", tin, tout, stepsize, nsteps, nd, slog);

    return;
}

/*** Rk2: propagates vin (of size n) from tin to tout by the Midpoint or ***
 *        Second-Order Runge-Kutta method; the result is returned by vout  *
 ***************************************************************************/

void Rk2::ps(double *vin, double *vout, double tin, double tout,
             double stephint, double accuracy, int n, FILE *slog)
{
    int i; /* local loop counter */

    double **v; /* intermediate v's, used to toggle v arrays */
    int toggle = 0; /* used to toggle between v[0] and v[1] */

    double *vtemp; /* guessed intermediate v's for midpoint */

    double *vnow; /* ptr to v at current time */
    double *vnext; /* ptr to v at current time + stepsize */

    double *deriv1; /* the derivatives at the beginning of a step */
    double *deriv2; /* the derivatives at midpoint of a step */

    double deltat; /* tout - tin */
    double t; /* current time */
    double thh; /* time at half of the step */

    double m; /* tmp var to calculate precise stepsize from stephint */
    double stepsize; /* real stepsize */
    double hh; /* half the stepsize */
    int step; /* loop counter for steps */
    int nsteps; /* total number of steps we have to take */

    int nd = 0; /* number of deriv evaluations */

    /* the do-nothing case; too small steps dealt with under usual */

    if (tin == tout)
        return;

    /* the usual case: steps big enough */

    v = (double **) calloc(2, sizeof(double *));

    vtemp = (double *) calloc(n, sizeof(double));
    v[0] = (double *) calloc(n, sizeof(double));
    v[1] = (double *) calloc(n, sizeof(double));
    deriv1 = (double *) calloc(n, sizeof(double));
    deriv2 = (double *) calloc(n, sizeof(double));

    deltat = tout - tin; /* how far do we have to propagate? */
    m = floor(deltat / stephint + 0.5); /* see comment on stephint above */
    if (m < 1.)
        m = 1.; /* we'll have to do at least one step */

    stepsize = deltat / m; /* real stepsize calculated here */
    nsteps = (int) m; /* int number of steps */
    t = tin; /* set current time */

    if (t == t + stepsize)
        error("Rk2: stephint of %g too small!", stephint);

    vnow = vin;

    if (nsteps == 1) /* vnext holds the results of current step */
        vnext = vout;
    else
        vnext = v[0];

    hh = stepsize * 0.5; /* evaluate half of stepsize */

    for (step = 0; step < nsteps; step++) { /* loop for steps */

        thh = t + hh; /* time of interval midpoints */

        zygote.p_deriv(vnow, t, deriv1, n); /* first call to deriv */

        if (debug)
            nd++;

        for (i = 0; i < n; i++) /* evaluate v's at midpoint */
            vtemp[i] = vnow[i] + hh * deriv1[i];

        zygote.p_deriv(vtemp, thh, deriv2, n); /* second call to deriv */

        if (debug)
            nd++;

        for (i = 0; i < n; i++)
            vnext[i] = vnow[i] + stepsize * deriv2[i]; /* Midpoint formula */

        t += stepsize; /* go to next step */

        if (step < nsteps - 2) { /* CASE 1: many steps to go */
            vnow = v[toggle]; /* toggle results from v[0] and v[1] */
            toggle++;
            toggle %= 2;
            vnext = v[toggle];

        } else if (step == nsteps - 2) { /* CASE 2: next iteration = final */
            vnow = v[toggle]; /* set vout */
            vnext = vout;

        } else if (step > nsteps - 2) { /* CASE 3: just did final iteration */
            free(v[0]); /* clean up and go home! */
            free(v[1]);
            free(v);
            free(vtemp);
            free(deriv1);
            free(deriv2);
        }
    }

    if (debug)
        WriteSolvLog("Rk2", tin, tout, stepsize, nsteps, nd, slog);

    return;
}

/*** Rk4: propagates vin (of size n) from tin to tout by the Fourth-Order **
 *        Runge-Kutta method; the result is returned by vout               *
 ***************************************************************************
 *                                                                         *
 * written by Joel Linton (somewhere around 1998)                          *
 * fixed and modified by Yoginho (somewhere around 2001)                   *
 *                                                                         *
 ***************************************************************************/
void Rk4::ps(double *vin, double *vout, double tin, double tout,
             double stephint, double accuracy, int n, FILE *slog)
{
    int i; /* local loop counter */

    double **v; /* intermediate v's, used to toggle v arrays */
    int toggle = 0; /* used to toggle between v[0] and v[1] */

    double *vtemp; /* guessed intermediate v's */

    double *vnow; /* ptr to v at current time */
    double *vnext; /* ptr to v at current time + stepsize */

    double *deriv1; /* the derivatives at the beginning of a step */
    double *deriv2; /* the derivatives at test-midpoint 1 */
    double *deriv3; /* the derivatives at test-midpoint 2 */
    double *deriv4; /* the derivatives at test-endpoint */

    double deltat; /* tout - tin */
    double t; /* current time */
    double th; /* time at the end of the step */
    double thh; /* time for interval midpoints */

    double m; /* used to calculate number of steps */
    double stepsize; /* real stepsize */
    double hh; /* half the stepsize */
    double h6; /* one sixth of the stepsize (for Rk4 formula) */
    int step; /* loop counter for steps */
    int nsteps; /* number of steps we have to take */

    int nd = 0; /* number of deriv evaluations */

    /* the do-nothing case; too small steps dealt with under usual */

    if (tin == tout)
        return;

    /* the usual case: steps big enough */

    v = (double **) calloc(2, sizeof(double *));

    vtemp = (double *) calloc(n, sizeof(double));
    v[0] = (double *) calloc(n, sizeof(double));
    v[1] = (double *) calloc(n, sizeof(double));
    deriv1 = (double *) calloc(n, sizeof(double));
    deriv2 = (double *) calloc(n, sizeof(double));
    deriv3 = (double *) calloc(n, sizeof(double));
    deriv4 = (double *) calloc(n, sizeof(double));

    deltat = tout - tin; /* how far do we have to propagate? */
    m = floor(deltat / stephint + 0.5); /* see comment on stephint above */
    if (m < 1.)
        m = 1.; /* we'll have to do at least one step */

    stepsize = deltat / m; /* real stepsize calculated here */
    nsteps = (int) m; /* int number of steps */
    t = tin; /* set current time */

    if (t == t + stepsize)
        error("Rk4: stephint of %g too small!", stephint);

    vnow = vin;

    if (nsteps == 1) /* vnext holds the results of current step */
        vnext = vout;
    else
        vnext = v[0];

    hh = stepsize * 0.5;
    h6 = stepsize / 6.0;

    for (step = 0; step < nsteps; step++) { /* loop for steps */

        thh = t + hh; /* time of interval midpoints */
        th = t + stepsize; /* time at end of the interval */

        /* do the Runge-Kutta thing here: first calculate intermediate derivatives */

        zygote.p_deriv(vnow, t, deriv1, n);

        if (debug)
            nd++;

        for (i = 0; i < n; i++)
            vtemp[i] = vnow[i] + hh * deriv1[i];
        zygote.p_deriv(vtemp, thh, deriv2, n);

        if (debug)
            nd++;

        for (i = 0; i < n; i++)
            vtemp[i] = vnow[i] + hh * deriv2[i];
        zygote.p_deriv(vtemp, thh, deriv3, n);

        if (debug)
            nd++;

        for (i = 0; i < n; i++)
            vtemp[i] = vnow[i] + stepsize * deriv3[i];
        zygote.p_deriv(vtemp, th, deriv4, n);

        if (debug)
            nd++;

        /* ... then feed them to the Fourth-Order Runge-Kutta formula */

        for (i = 0; i < n; i++)
            vnext[i] = vnow[i]
                    + h6 * (deriv1[i] + 2.0 * deriv2[i] + 2.0 * deriv3[i]
                            + deriv4[i]);

        /* next step */

        t += stepsize;

        if (step < nsteps - 2) { /* CASE 1: many steps to go */
            vnow = v[toggle]; /* toggle results from v[0] and v[1] */
            toggle++;
            toggle %= 2;
            vnext = v[toggle];

        } else if (step == nsteps - 2) { /* CASE 2: next iteration = final */
            vnow = v[toggle]; /* set vout */
            vnext = vout;

        } else if (step > nsteps - 2) { /* CASE 3: just did final iteration */
            free(v[0]); /* clean up and go home! */
            free(v[1]);
            free(v);
            free(vtemp);
            free(deriv1);
            free(deriv2);
            free(deriv3);
            free(deriv4);
        }
    }

    if (debug)
        WriteSolvLog("Rk4", tin, tout, stepsize, nsteps, nd, slog);

    return;
}

/*** Rkck: propagates vin (of size n) from tin to tout by the Runge-Kutta **
 *         Cash-Karp method, which is an adaptive-stepsize Rk method; it   *
 *         uses a fifth-order Rk formula with an embedded forth-oder for-  *
 *         mula for calucalting the error; its result is returned by vout  *
 ***************************************************************************
 *                                                                         *
 * This solver was written by Marcel Wolf, Spring 2002.                    *
 *                                                                         *
 ***************************************************************************/

void Rkck::ps(double *vin, double *vout, double tin, double tout,
              double stephint, double accuracy, int n, FILE *slog)
{
    int i; /* local loop counter */

    double **v; /** used for storing intermediate steps */
    int toggle = 0; /* used to toggle between v[0] and v[1] */

    double *vtemp; /* guessed intermediate v's */

    double *vnow; /* ptr to v at current time */
    double *vnext; /* ptr to v at current time + stepsize */

    double *verror; /* error estimate */
    double verror_max; /* the maximum error in verror */

    double *deriv1; /* intermediate derivatives for the Cash-Karp formula */
    double *deriv2;
    double *deriv3;
    double *deriv4;
    double *deriv5;
    double *deriv6;

    double t; /* the current time */

    double h = stephint; /* initial stepsize */
    double hnext; /* used to calculate next stepsize */
    const double SAFETY = 0.9; /* safety margin for decreasing stepsize */

    /* declare and initialize Cash-Karp parameters */
    /* parameters that are zero have been added as comments for clarity */

    static double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875;

    static double b21 = 0.2;
    static double b31 = 3.0 / 40.0;
    static double b32 = 9.0 / 40.0;
    static double b41 = 0.3;
    static double b42 = -0.9;
    static double b43 = 1.2;
    static double b51 = -11.0 / 54.0;
    static double b52 = 2.5;
    static double b53 = -70.0 / 27.0;
    static double b54 = 35.0 / 27.0;
    static double b61 = 1631.0 / 55296.0;
    static double b62 = 175.0 / 512.0;
    static double b63 = 575.0 / 13824;
    static double b64 = 44275.0 / 110592.0;
    static double b65 = 253.0 / 4096.0;

    static double c1 = 37.0 / 378.0;
    /*  static double c2 = 0.0; */
    static double c3 = 250.0 / 621.0;
    static double c4 = 125.0 / 594.0;
    /*  static double c5 = 0.0; */
    static double c6 = 512.0 / 1771.0;

    double dc1 = c1 - 2825.0 / 27648.0;
    /*  double dc2 = 0.0; */
    double dc3 = c3 - 18575.0 / 48384.0;
    double dc4 = c4 - 13525.0 / 55296.0;
    double dc5 = -277.0 / 14336.0;
    double dc6 = c6 - 0.25;

    /* the do-nothing case; too small steps dealt with under usual */

    if (tin == tout)
        return;

    /* the usual case: steps big enough */

    v = (double **) calloc(2, sizeof(double *));

    vtemp = (double *) calloc(n, sizeof(double));
    verror = (double *) calloc(n, sizeof(double));
    v[0] = (double *) calloc(n, sizeof(double));
    v[1] = (double *) calloc(n, sizeof(double));
    deriv1 = (double *) calloc(n, sizeof(double));
    deriv2 = (double *) calloc(n, sizeof(double));
    deriv3 = (double *) calloc(n, sizeof(double));
    deriv4 = (double *) calloc(n, sizeof(double));
    deriv5 = (double *) calloc(n, sizeof(double));
    deriv6 = (double *) calloc(n, sizeof(double));

    t = tin;
    vnow = vin;
    vnext = v[0];

    /* initial stepsize cannot be bigger than total time */

    if (tin + h >= tout)
        h = tout - tin;

    while (t < tout) {

        /* Take one step and evaluate the error. Repeat until the resulting error  *
         * is less than the desired accuracy                                       */

        while (1) {

            /* do the Runge-Kutta thing here: calulate intermediate derivatives */

            zygote.p_deriv(vnow, t, deriv1, n);

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i] + h * (b21 * deriv1[i]);
            zygote.p_deriv(vtemp, t + a2 * h, deriv2, n);

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i] + h * (b31 * deriv1[i] + b32 * deriv2[i]);
            zygote.p_deriv(vtemp, t + a3 * h, deriv3, n);

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i]
                        + h * (b41 * deriv1[i] + b42 * deriv2[i]
                               + b43 * deriv3[i]);
            zygote.p_deriv(vtemp, t + a4 * h, deriv4, n);

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i]
                        + h * (b51 * deriv1[i] + b52 * deriv2[i]
                               + b53 * deriv3[i] + b54 * deriv4[i]);
            zygote.p_deriv(vtemp, t + a5 * h, deriv5, n);

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i]
                        + h * (b61 * deriv1[i] + b62 * deriv2[i]
                               + b63 * deriv3[i] + b64 * deriv4[i]
                               + b65 * deriv5[i]);
            zygote.p_deriv(vtemp, t + a6 * h, deriv6, n);

            /* ... then feed them to the Cash-Karp formula */

            for (i = 0; i < n; i++)
                vnext[i] = vnow[i]
                        + h * (c1 * deriv1[i] + c3 * deriv3[i] + c4 * deriv4[i]
                               + c6 * deriv6[i]);

            /* calculate the error estimate using the embedded formula */

            for (i = 0; i < n; i++)
                verror[i] = h
                        * (dc1 * deriv1[i] + dc3 * deriv3[i] + dc4 * deriv4[i]
                           + dc5 * deriv5[i] + dc6 * deriv6[i]);

            /* find the maximum error */

            verror_max = 0.;
            for (i = 0; i < n; i++)
                if (vnext[i] != 0.)
                    verror_max = max(fabs(verror[i] / vnext[i]), verror_max);
                else
                    verror_max = max(verror[i] / DBL_EPSILON, verror_max);

            /* scale error according to desired accuracy */

            verror_max /= accuracy;

            /* compare maximum error to the desired accuracy; if error < accuracy, we  *
             * are done with this step; otherwise, the stepsize has to be reduced and  *
             * the step repeated; for detailed comments on the approximation involving *
             * SAFETY and -0.25 see 'Numerical Recipes in C', 2nd Edition, p.718       */

            if (verror_max <= 1.0)
                break;

            hnext = SAFETY * h * pow(verror_max, -0.25);

            /* decrease stepsize by no more than a factor of 10; check for underflows */

            h = (hnext > 0.1 * h) ? hnext : 0.1 * h;
            if (h == 0)
                error("Rkck: stepsize underflow");

        }

        /* advance the current time by last stepsize */

        t += h;

        if (t >= tout)
            break; /* that was the last iteration */

        /* increase stepsize according to error (5th order) for next iteration */

        h = h * pow(verror_max, -0.20);

        /* make sure t does not overstep tout */

        if (t + h >= tout)
            h = tout - t;

        /* toggle vnow and vnext between v[0] and v[1] */

        vnow = v[toggle];
        toggle++;
        toggle %= 2;
        vnext = v[toggle];

    }

    /* copy the last result to vout after the final iteration */

    memcpy(vout, vnext, sizeof(*vnext) * n);

    free(v[0]);
    free(v[1]);
    free(v);
    free(vtemp);
    free(verror);
    free(deriv1);
    free(deriv2);
    free(deriv3);
    free(deriv4);
    free(deriv5);
    free(deriv6);

}

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

void Rkf::ps(double *vin, double *vout, double tin, double tout,
             double stephint, double accuracy, int n, FILE *slog)
{
    int i; /* local loop counter */

    double **v; /* used for storing intermediate steps */
    int toggle = 0; /* used to toggle between v[0] and v[1] */

    double *vtemp; /* guessed intermediate v's */

    double *vnow; /* ptr to v at current time */
    double *vnext; /* ptr to v at current time + stepsize */

    double *verror; /* error estimate */
    double verror_max; /* the maximum error in verror */

    double *deriv1; /* intermediate derivatives for the Fehlberg formula */
    double *deriv2;
    double *deriv3;
    double *deriv4;
    double *deriv5;
    double *deriv6;

    double t; /* the current time */

    double h = stephint; /* initial stepsize */
    double hnext; /* used to calculate next stepsize */
    const double SAFETY = 0.9; /* safety margin for decreasing stepsize */

    /* declare and initialize Fehlberg parameters */
    /* parameters that are zero have been added as comments for clarity */

    static double a2 = 0.25, a3 = 0.375, a4 = 12.0 / 13.0, a5 = 1.0, a6 = 0.5;

    static double b21 = 0.25, b31 = 3.0 / 32.0, b32 = 9.0 / 32.0, b41 = 1932.0
            / 2197.0, b42 = -7200.0 / 2197.0, b43 = 7296.0 / 2197.0, b51 = 439.0
            / 216.0, b52 = -8.0, b53 = 3680.0 / 513.0, b54 = -845.0 / 4104.0,
            b61 = -8.0 / 27.0, b62 = 2.0, b63 = -3544.0 / 2565.0, b64 = 1859.0
                    / 4104.0, b65 = -11.0 / 40.0;

    static double c1 = 25.0 / 216.0,
    /*  c2  =     0.0 */
    c3 = 1408.0 / 2565.0, c4 = 2197.0 / 4104.0, c5 = -0.2;
    /*  c6  =     0.0 */

    static double dc1 = 1.0 / 360.0,
    /*  dc2 =     0.0 */
    dc3 = -128.0 / 4275.0, dc4 = -2197.0 / 75240.0, dc5 = 1.0 / 50.0, dc6 = 2.0
            / 55.0;

    /* the do-nothing case; too small steps dealt with under usual */

    if (tin == tout)
        return;

    /* the usual case: steps big enough */

    v = (double **) calloc(2, sizeof(double *));

    vtemp = (double *) calloc(n, sizeof(double));
    verror = (double *) calloc(n, sizeof(double));
    v[0] = (double *) calloc(n, sizeof(double));
    v[1] = (double *) calloc(n, sizeof(double));
    deriv1 = (double *) calloc(n, sizeof(double));
    deriv2 = (double *) calloc(n, sizeof(double));
    deriv3 = (double *) calloc(n, sizeof(double));
    deriv4 = (double *) calloc(n, sizeof(double));
    deriv5 = (double *) calloc(n, sizeof(double));
    deriv6 = (double *) calloc(n, sizeof(double));

    t = tin;
    vnow = vin;
    vnext = v[0];

    /* initial stepsize cannot be bigger than total time */

    if (tin + h >= tout)
        h = tout - tin;

    while (t < tout) {

        /* Take one step and evaluate the error. Repeat until the resulting error  *
         * is less than the desired accuracy                                       */

        while (1) {

            /* do the Runge-Kutta thing here: calulate intermediate derivatives */

            zygote.p_deriv(vnow, t, deriv1, n);

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i] + h * (b21 * deriv1[i]);
            zygote.p_deriv(vtemp, t + a2 * h, deriv2, n);

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i] + h * (b31 * deriv1[i] + b32 * deriv2[i]);
            zygote.p_deriv(vtemp, t + a3 * h, deriv3, n);

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i]
                        + h * (b41 * deriv1[i] + b42 * deriv2[i]
                               + b43 * deriv3[i]);
            zygote.p_deriv(vtemp, t + a4 * h, deriv4, n);

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i]
                        + h * (b51 * deriv1[i] + b52 * deriv2[i]
                               + b53 * deriv3[i] + b54 * deriv4[i]);
            zygote.p_deriv(vtemp, t + a5 * h, deriv5, n);

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i]
                        + h * (b61 * deriv1[i] + b62 * deriv2[i]
                               + b63 * deriv3[i] + b64 * deriv4[i]
                               + b65 * deriv5[i]);
            zygote.p_deriv(vtemp, t + a6 * h, deriv6, n);

            /* ... then feed them to the Fehlberg formula */

            for (i = 0; i < n; i++)
                vnext[i] = vnow[i]
                        + h * (c1 * deriv1[i] + c3 * deriv3[i] + c4 * deriv4[i]
                               + c5 * deriv5[i]);

            /* calculate the error estimate using the embedded formula */

            for (i = 0; i < n; i++)
                verror[i] = fabs(
                        h * (dc1 * deriv1[i] + dc3 * deriv3[i] + dc4 * deriv4[i]
                             + dc5 * deriv5[i] + dc6 * deriv6[i]));

            /* find the maximum error */

            verror_max = 0.;
            for (i = 0; i < n; i++)
                if (vnext[i] != 0.)
                    verror_max = max(verror[i] / vnext[i], verror_max);
                else
                    verror_max = max(verror[i] / DBL_EPSILON, verror_max);

            /* scale error according to desired accuracy */

            verror_max /= accuracy;

            /* compare maximum error to the desired accuracy; if error < accuracy, we  *
             * are done with this step; otherwise, the stepsize has to be reduced and  *
             * the step repeated; for detailed comments on the approximation involving *
             * SAFETY and -0.25 see 'Numerical Recipes in C', 2nd Edition, p.718       */

            if (verror_max <= 1.0)
                break;

            hnext = SAFETY * h * pow(verror_max, -0.25);

            /* decrease stepsize by no more than a factor of 10; check for underflows */

            h = (hnext > 0.1 * h) ? hnext : 0.1 * h;
            if (h == 0)
                error("Rkf: stepsize underflow");

        }

        /* advance the current time by last stepsize */

        t += h;

        if (t >= tout)
            break; /* that was the last iteration */

        /* increase stepsize according to error for next iteration */

        h = h * pow(verror_max, -0.20);

        /* make sure t does not overstep tout */

        if (t + h >= tout)
            h = tout - t;

        /* toggle vnow and vnext between v[0] and v[1] */

        vnow = v[toggle];
        toggle++;
        toggle %= 2;
        vnext = v[toggle];

    }

    /* copy the last result to vout after the final iteration */

    memcpy(vout, vnext, sizeof(*vnext) * n);

    free(v[0]);
    free(v[1]);
    free(v);
    free(vtemp);
    free(verror);
    free(deriv1);
    free(deriv2);
    free(deriv3);
    free(deriv4);
    free(deriv5);
    free(deriv6);

}

/***** Milne: propagates vin (of size n) from tin to tout by Milne-Simpson *
 *            which is a predictor-corrector method; the result is retur-  *
 *            ned by vout                                                  *
 ***************************************************************************
 *                                                                         *
 * This solver was implemented by Konstantin Koslov, Dec 2001/Jan 2002     *
 *                                                                         *
 ***************************************************************************/

void Milne::ps(double *vin, double *vout, double tin, double tout,
               double stephint, double accuracy, int n, FILE *slog)
{
    int i; /* local loop counter */

    double **v; /* array to store intermediate results */
    int toggle = 0; /* used to toggle between v[0] and v[1] */

    double *vtemp; /* guessed intermediate v's */

    double *vnow; /* current protein concs */
    double *vnext; /* next protein concs */

    double *deriv1; /* the derivatives at the beginning of a step */
    double *deriv2; /* derivatives at test-midpoint 1 */
    double *deriv3; /* derivatives at test-midpoint 2 */
    double *deriv4; /* derivatives at test-endpoint */

    double **history_dv; /* history of the derivatives */
    double **history_v; /* history of the v: v_{i-3} v_{i-1} */
    double *dblank; /* pointer for history_dv rotation */

    double deltat; /* tout - tin */
    double t; /* current time */
    double th; /* time at end of interval */
    double thh; /* time for interval midpoints */

    double m; /* used to calculate number of steps */
    double stepsize; /* real stepsize */
    double hh; /* half the stepsize */
    double h6; /* one sixth of the stepsize (for Rk4 formula) */
    int step; /* loop counter for steps */
    int nsteps; /* number of steps we have to take */

    double hp; /* multiplier for the predictor derivs */
    double hc; /* multiplier for the corrector derivs */

    double mistake; /* error estimator */
    double corr; /* temp variable for the corrector */
    double ee; /* temp variable for the error estimator */
    double cc29; /* 1/29: used for calculating ee */

    int nd = 0; /* number of deriv evaluations */

    /* the do-nothing case; too small steps dealt with under usual */

    if (tin == tout) {
        return;
    }

    /* the usual case: steps big enough */

    v = (double **) calloc(2, sizeof(double *));

    vtemp = (double *) calloc(n, sizeof(double));
    v[0] = (double *) calloc(n, sizeof(double));
    v[1] = (double *) calloc(n, sizeof(double));
    deriv1 = (double *) calloc(n, sizeof(double));
    deriv2 = (double *) calloc(n, sizeof(double));
    deriv3 = (double *) calloc(n, sizeof(double));
    deriv4 = (double *) calloc(n, sizeof(double));

    deltat = tout - tin; /* how far do we have to propagate? */
    m = floor(deltat / stephint + 0.5); /* see comment on stephint above */

    if (m < 1.)
        m = 1.; /* we'll have to do at least one step */

    stepsize = deltat / m; /* real stepsize calculated here */
    nsteps = (int) m; /* int number of steps */
    t = tin; /* set current time */

    if (t == t + stepsize)
        error("Milne: stephint of %g too small!", stephint);

    vnow = vin;

    if (nsteps == 1) /* vnext holds the results of current step */
        vnext = vout;
    else
        vnext = v[0];

    hh = stepsize * 0.5;
    h6 = stepsize / 6.0;

    /* do we need Milne? no, not for less than 4 steps! */

    if (nsteps < 4) {

        for (step = 0; step < nsteps; step++) { /* loop for steps */

            thh = t + hh; /* time of interval midpoints */
            th = t + stepsize; /* time at end of the interval */

            /* do the Runge-Kutta thing here: first calculate intermediate derivatives */

            zygote.p_deriv(vnow, t, deriv1, n);

            if (debug)
                nd++;

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i] + hh * deriv1[i];
            zygote.p_deriv(vtemp, thh, deriv2, n);

            if (debug)
                nd++;

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i] + hh * deriv2[i];
            zygote.p_deriv(vtemp, thh, deriv3, n);

            if (debug)
                nd++;

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i] + stepsize * deriv3[i];
            zygote.p_deriv(vtemp, th, deriv4, n);

            if (debug)
                nd++;

            /* ... then feed them to the Fourth-Order Runge-Kutta formula */

            for (i = 0; i < n; i++)
                vnext[i] = vnow[i]
                        + h6 * (deriv1[i] + 2.0 * deriv2[i] + 2.0 * deriv3[i]
                                + deriv4[i]);

            /* next step */

            t += stepsize;

            if (step < nsteps - 2) { /* CASE 1: many steps to go */
                vnow = v[toggle]; /* toggle results from v[0] and v[1] */
                toggle++;
                toggle %= 2;
                vnext = v[toggle];

            } else if (step == nsteps - 2) { /* CASE 2: next iteration = final */
                vnow = v[toggle]; /* set vout */
                vnext = vout;

            } else if (step > nsteps - 2) { /* CASE 3: just did final iteration */
                free(v[0]); /* clean up and go home! */
                free(v[1]);
                free(v);
                free(deriv1);
                free(deriv2);
                free(deriv3);
                free(deriv4);
                free(vtemp);
            }
        }

        /* more than 3 steps: yes, we need Milne */

    } else {

        history_dv = (double **) calloc(3, sizeof(double *));
        history_v = (double **) calloc(4, sizeof(double *));

        for (step = 0; step < 3; step++) {
            history_dv[step] = (double *) calloc(n, sizeof(double));
        }

        for (step = 0; step < 4; step++) {
            history_v[step] = (double *) calloc(n, sizeof(double));
        }

        mistake = 0.;
        cc29 = 1. / 29.;

        /* loop for initial steps using Rk4 */

        for (step = 0; step < 3; step++) {

            thh = t + hh; /* time of interval midpoints */
            th = t + stepsize; /* time at end of the interval */

            /* do the Runge-Kutta thing here: first calculate intermediate derivatives */

            zygote.p_deriv(vnow, t, deriv1, n);

            if (debug)
                nd++;

            if (step > 0)
                for (i = 0; i < n; i++)
                    history_dv[step - 1][i] = deriv1[i];

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i] + hh * deriv1[i];
            zygote.p_deriv(vtemp, thh, deriv2, n);

            if (debug)
                nd++;

            for (i = 0; i < n; i++) /* evaluate deriv at second midpoint */
                vtemp[i] = vnow[i] + hh * deriv2[i];
            zygote.p_deriv(vtemp, thh, deriv3, n);

            if (debug)
                nd++;

            for (i = 0; i < n; i++) /* evaluate deriv at end point */
                vtemp[i] = vnow[i] + stepsize * deriv3[i];
            zygote.p_deriv(vtemp, th, deriv4, n);

            if (debug)
                nd++;

            /* ... then feed them to the Fourth-Order Runge-Kutta formula */

            for (i = 0; i < n; i++)
                vnext[i] = vnow[i]
                        + h6 * (deriv1[i] + 2.0 * deriv2[i] + 2.0 * deriv3[i]
                                + deriv4[i]);

            /* next step */

            t += stepsize;

            for (i = 0; i < n; i++) /* save old v values in history_v */
                history_v[step][i] = vnow[i];

            if (step == 1) {
                vnow = vnext;
                vnext = vout;
            } else {
                vnow = v[toggle]; /* toggle results from v[0] and v[1] */
                toggle++;
                toggle %= 2;
                vnext = v[toggle];
            }
        }

        /* we have calculated 4 initial points with rk4: start multistepping */

        hc = stepsize / 3.0; /* stepsize multipliers: for the corrector */
        hp = 4. * hc; /* ... and for the predictor */

        for (i = 0; i < n; i++) /* save current values of v in history_v */
            history_v[3][i] = vout[i];

        for (step = 3; step < nsteps; step++) { /* loop for Milne steps */

            th = t + stepsize; /* time at end of the interval */

            /* do the Milne thing here */

            zygote.p_deriv(vout, t, deriv1, n); /* evaluate deriv at start of interval */

            if (debug)
                nd++;

            for (i = 0; i < n; i++) /* save current derivs in history_dv */
                history_dv[2][i] = deriv1[i];

            /* Predictor step (P): extrapolate derivatives */

            for (i = 0; i < n; i++)
                v[0][i] = history_v[0][i]
                        + hp * (2 * history_dv[0][i] - history_dv[1][i] + 2
                                * history_dv[2][i]);

            /* Evaluate the extrapolated derivatives (E) */

            zygote.p_deriv(v[0], th, deriv1, n);

            if (debug)
                nd++;

            /* Corrector step (C): interpolate current v (+ calculate error estimator) */

            for (i = 0; i < n; i++) {

                corr = history_v[2][i] + hc
                        * (history_dv[1][i] + 4 * history_dv[2][i] + deriv1[i]);
                ee = fabs(corr - v[0][i]) * cc29; /* error estimator */

                if (ee > mistake)
                    mistake = ee;
                vout[i] = corr;

            }

            t += stepsize; /* go to next step */

            /* update the arrays */

            if (step <= nsteps - 2) { /* CASE 1: many steps to go */

                dblank = history_dv[0];
                history_dv[0] = history_dv[1];
                history_dv[1] = history_dv[2];
                history_dv[2] = dblank;

                dblank = history_v[0];
                history_v[0] = history_v[1];
                history_v[1] = history_v[2];
                history_v[2] = history_v[3];
                history_v[3] = dblank;

                for (i = 0; i < n; i++)
                    history_v[3][i] = vout[i];

            } else { /* CASE 2: just did final iteration */

                free(v[0]); /* clean up and go home! */
                free(v[1]);
                free(v);
                free(history_v[0]);
                free(history_v[1]);
                free(history_v[2]);
                free(history_v[3]);
                free(history_v);
                free(history_dv[0]);
                free(history_dv[1]);
                free(history_dv[2]);
                free(history_dv);
                free(deriv1);
                free(deriv2);
                free(deriv3);
                free(deriv4);
                free(vtemp);
            }
        }

        /* use below for sanity check if Milne behaves unstably */

        /*
         fprintf(stdout,"mistake=%.15f\n", mistake);
         fflush(stdout);
         getchar();
         */

    }

    if (debug)
        WriteSolvLog("Milne", tin, tout, stepsize, nsteps, nd, slog);

    return;

}

/***** Adams: propagates vin (of size n) from tin to tout by Adams-Moulton *
 *            which is an implicit predictor-corrector method of second    *
 *            order; the result is returned by vout                        *
 ***************************************************************************
 *                                                                         *
 * This solver was implemented by Konstantin Koslov, Spring 2002           *
 * Slightly modified by Manu, July 2002                                    *
 *                                                                         *
 ***************************************************************************/

void Adams::ps(double *vin, double *vout, double tin, double tout,
               double stephint, double accuracy, int n, FILE *slog)
{
    int i; /* local loop counter */

    double **v; /* array to store intermediate results */
    int toggle = 0; /* used to toggle between v[0] and v[1] */

    double *vtemp; /* guessed intermediate v's */

    double *vnow; /* current protein concs */
    double *vnext; /* next protein concs */

    double *deriv1; /* the derivatives at the beginning of a step */
    double *deriv2; /* derivatives at test-midpoint 1 */
    double *deriv3; /* derivatives at test-midpoint 2 */
    double *deriv4; /* derivatives at test-endpoint */

    double **history_dv; /* history of the derivatives */
    double *dblank; /* pointer for history_dv rotation */

    double deltat; /* tout - tin */
    double t; /* current time */
    double th; /* time at end of interval */
    double thh; /* time for interval midpoints */

    double m; /* used to calculate number of steps */
    double stepsize; /* real stepsize */
    double hh; /* half the stepsize */
    double h6; /* one sixth of the stepsize (for Rk4 formula) */
    int step; /* loop counter for steps */
    int nsteps; /* number of steps we have to take */

    double mistake; /* error estimators */
    double eps1;
    double eps2;

    int nd = 0; /* number of deriv evaluations */

    /* the do-nothing case; too small steps dealt with under usual */

    if (tin == tout)
        return;

    /* the usual case: steps big enough */

    v = (double **) calloc(2, sizeof(double *));

    vtemp = (double *) calloc(n, sizeof(double));
    v[0] = (double *) calloc(n, sizeof(double));
    v[1] = (double *) calloc(n, sizeof(double));
    deriv1 = (double *) calloc(n, sizeof(double));
    deriv2 = (double *) calloc(n, sizeof(double));
    deriv3 = (double *) calloc(n, sizeof(double));
    deriv4 = (double *) calloc(n, sizeof(double));

    deltat = tout - tin; /* how far do we have to propagate? */
    m = floor(deltat / stephint + 0.5); /* see comment on stephint above */

    if (m < 1.)
        m = 1.; /* we'll have to do at least one step */

    stepsize = deltat / m; /* real stepsize calculated here */
    nsteps = (int) m; /* int number of steps */
    t = tin; /* set current time */

    if (t == t + stepsize)
        error("Adams: stephint of %g too small!", stephint);

    vnow = vin;

    if (nsteps == 1) /* vnext holds the results of current step */
        vnext = vout;
    else
        vnext = v[0];

    hh = stepsize * 0.5;
    h6 = stepsize / 6.0;

    /* do we need Adams? no, not for less than 4 steps! */

    if (nsteps < 4) {

        for (step = 0; step < nsteps; step++) { /* loop for steps */

            thh = t + hh; /* time of interval midpoints */
            th = t + stepsize; /* time at end of the interval */

            /* do the Runge-Kutta thing here: first calculate intermediate derivatives */

            zygote.p_deriv(vnow, t, deriv1, n);

            if (debug)
                nd++;

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i] + hh * deriv1[i];
            zygote.p_deriv(vtemp, thh, deriv2, n);

            if (debug)
                nd++;

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i] + hh * deriv2[i];
            zygote.p_deriv(vtemp, thh, deriv3, n);

            if (debug)
                nd++;

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i] + stepsize * deriv3[i];
            zygote.p_deriv(vtemp, th, deriv4, n);

            if (debug)
                nd++;

            /* ... then feed them to the Fourth-Order Runge-Kutta formula */

            for (i = 0; i < n; i++)
                vnext[i] = vnow[i]
                        + h6 * (deriv1[i] + 2.0 * deriv2[i] + 2.0 * deriv3[i]
                                + deriv4[i]);

            /* next step */

            t += stepsize;

            if (step < nsteps - 2) { /* CASE 1: many steps to go */
                vnow = v[toggle]; /* toggle results from v[0] and v[1] */
                toggle++;
                toggle %= 2;
                vnext = v[toggle];

            } else if (step == nsteps - 2) { /* CASE 2: next iteration = final */
                vnow = v[toggle]; /* set vout */
                vnext = vout;

            } else if (step > nsteps - 2) { /* CASE 3: just did final iteration */
                free(v[0]); /* clean up and go home! */
                free(v[1]);
                free(v);
                free(vtemp);
                free(deriv1);
                free(deriv2);
                free(deriv3);
                free(deriv4);
            }
        }

    } else {

        /* more than 3 steps: yes, we need Adams*/

        history_dv = (double **) calloc(4, sizeof(double *));
        for (step = 0; step < 4; step++) {
            history_dv[step] = (double *) calloc(n, sizeof(double));
        }

        mistake = 0.;

        /* loop for initial steps using Rk4 */

        for (step = 0; step < 3; step++) {

            thh = t + hh; /* time of interval midpoints */
            th = t + stepsize; /* time at end of the interval */

            /* do the Runge-Kutta thing here: first calculate intermediate derivatives */

            zygote.p_deriv(vnow, t, deriv1, n);

            if (debug)
                nd++;

            for (i = 0; i < n; i++)
                history_dv[step][i] = deriv1[i];

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i] + hh * deriv1[i];
            zygote.p_deriv(vtemp, thh, deriv2, n);

            if (debug)
                nd++;

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i] + hh * deriv2[i];
            zygote.p_deriv(vtemp, thh, deriv3, n);

            if (debug)
                nd++;

            for (i = 0; i < n; i++)
                vtemp[i] = vnow[i] + stepsize * deriv3[i];
            zygote.p_deriv(vtemp, th, deriv4, n);

            if (debug)
                nd++;

            /* ... then feed them to the Fourth-Order Runge-Kutta formula */

            for (i = 0; i < n; i++)
                vnext[i] = vnow[i]
                        + h6 * (deriv1[i] + 2.0 * deriv2[i] + 2.0 * deriv3[i]
                                + deriv4[i]);

            /* next step */

            t += stepsize;

            vnow = v[toggle]; /* toggle results from v[0] and v[1] */
            toggle++;
            toggle %= 2;
            vnext = v[toggle];

        }

        /* we have calculated 4 initial points with rk4: start multistepping */

        hh = stepsize / 24.0; /* reset hh to 1/24 of stepsize */

        if (nsteps == 4)
            vnext = vout;

        for (step = 3; step < nsteps; step++) { /* loop for steps */

            th = t + stepsize; /* time at end of the interval */

            /* do the Adams thing here */

            zygote.p_deriv(vnow, t, deriv1, n); /* evaluate deriv at start of interval */

            if (debug)
                nd++;

            for (i = 0; i < n; i++)
                history_dv[3][i] = deriv1[i];

            for (i = 0; i < n; i++)
                vnext[i] =
                        vnow[i] + hh
                                * (55 * history_dv[3][i] - 59 * history_dv[2][i] + 37
                                        * history_dv[1][i]
                                   - 9 * history_dv[0][i]);

            zygote.p_deriv(vnext, th, deriv1, n);

            if (debug)
                nd++;

            /* the following evaluates the error, but then nothing is done with it */

            /*    for(i=0; i<n; i++) {
             eps1 = hh * (           deriv1[i] - 3 * history_dv[3][i]
             + 3 * history_dv[2][i] -     history_dv[1][i]);
             eps2 = fabs(eps1);
             if (eps2 > mistake)
             mistake=eps2;
             }
             */
            for (i = 0; i < n; i++)
                vnext[i] = vnow[i]
                        + hh * (9 * deriv1[i] + 19 * history_dv[3][i] - 5
                                * history_dv[2][i]
                                + history_dv[1][i]);

            t += stepsize; /* go to next step */

            if (step <= nsteps - 2) { /* CASE 1: many steps to go */
                dblank = history_dv[0];
                history_dv[0] = history_dv[1];
                history_dv[1] = history_dv[2];
                history_dv[2] = history_dv[3];
                history_dv[3] = dblank;

                if (step < nsteps - 2) {
                    vnow = v[toggle]; /* toggle results from v[0] and v[1] */
                    toggle++;
                    toggle %= 2;
                    vnext = v[toggle];

                } else if (step == nsteps - 2) { /* CASE 2: next iteration = final */
                    vnow = v[toggle]; /* set vout */
                    vnext = vout;
                }

            } else { /* CASE 2: just did final iteration */

                free(v[0]); /* clean up and go home! */
                free(v[1]);
                free(v);
                free(history_dv[0]);
                free(history_dv[1]);
                free(history_dv[2]);
                free(history_dv[3]);
                free(history_dv);
                free(deriv1);
                free(deriv2);
                free(deriv3);
                free(deriv4);
                free(vtemp);
            }
        }
    }

    if (debug)
        WriteSolvLog("Adams", tin, tout, stepsize, nsteps, nd, slog);

    return;
}

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

void BuSt::ps(double *vin, double *vout, double tin, double tout,
              double stephint, double accuracy, int n, FILE *slog)
{

    int i; /* local loop counter */

    double t; /* current time */

    double ht; /* stepsize of next step */
    double hd; /* stepsize of step we just did */

    double *initderiv; /* array for initial derivatives */

    /* the do-nothing case; too small steps dealt with under usual */

    if (tin == tout)
        return;

    /* the usual case: steps big enough -> initialize variables */

    initderiv = (double *) calloc(n, sizeof(double));

    t = tin;

    /* initial stepsize is either stephint or tout-tin */

    ht = min(stephint, (tout - tin));

    /* we're saving old v's in vin and propagate vout below */

    for (i = 0; i < n; i++)
        vout[i] = vin[i];

    /* this loop steps through the whole interval tout-tin */

    do {

        zygote.p_deriv(vout, t, initderiv, n);
        bsstep(vout, initderiv, n, &t, ht, accuracy, &hd, &ht);

        if (t + ht > tout)
            ht = tout - t;

    } while (t < tout);

    /* clean up and go home */

    free(initderiv);

    return;
}

/*** bsstep: this is the Bulirsch-Stoer stepper function; it propagates ****
 *           the equations at a given overall stepsize and order, then     *
 *           evaluates the error and repeats the step with increased order *
 *           if necessary, until the required accuracy is achieved; its    *
 *           arguments are:                                                *
 *           - v:        input and output v's over a give total stepsize   *
 *           - deriv:    array of initial derivatives (at tin)             *
 *           - n:        size of v and deriv                               *
 *           - t:        pointer to the current time                       *
 *           - htry:     initial stepsize to be tried                      *
 *           - accuracy: relative accuracy to be achieved                  *
 *           - hdid:     returns the stepsize of the step we just did      *
 *           - hnext:    returns the stepsize of the next step to be done  *
 ***************************************************************************/
void BuSt::bsstep(double *v, double *deriv, int n, double *t, double htry,
                  double accuracy, double *hdid, double *hnext)
{

    /*** constants *************************************************************/

    const int KMAXX = 8; /* this is how many successive stepsize */
    const int IMAXX = KMAXX + 1; /* refinements we're going to try */

    const double SAFE1 = 0.25; /* safety margins for errors and h reduction */
    const double SAFE2 = 0.7;

    const double REDMAX = 1.0e-5; /* boundaries for stepsize reduction */
    const double REDMIN = 0.7;

    const double SCALMX = 0.1; /* maximum scaling value for stepsize */

    /*** variables *************************************************************/

    int i, k, km; /* loop counters */

    double *vsav; /* v's at beginning of interval */
    double *vseq; /* v's at end of interval for a given stepsize */

    double *verror; /* error estimates */
    double verror_max; /* the maximum error in verror */
    double *err; /* normalized error estimates */

    double h; /* stepsize for solver method used below */
    double hest; /* square of h passed to extrapolation function */
    double red; /* reduction factor for reducing stepsize */
    double scale; /* scaling factor for next stepsize */

    double work; /* temp vars for calculating row for convergence */
    double wrkmin;
    double fact;

    static double old_accuracy = -1.0; /* used to save old accuracy */
    static double tnew; /* used to save old start time */

    /* the following two arrays are used for Deuflhard's error estimation; a   *
     * contains the work coefficients and alf (alpha) the correction factors   */

    static double *a = NULL;
    static double **alf = NULL;

    double accuracy1; /* error (< accuracy) used to calculate alphas */

    static int kmax; /* used for finding kopt */
    static int kopt; /* optimal row number for convergence */

    /* sequence of separate attempts to cross interval htot with increasing    *
     * values of nsteps as suggested by Deuflhard (Num Rec, p. 726)            */

    static int nseq[] = { 0, 2, 4, 6, 8, 10, 12, 14, 16, 18 };

    /* miscellaneous flags */

    static int first = 1; /* is this the first try for a given step? */
    int reduct; /* flag indicating if we have reduced stepsize yet */
    int exitflag = 0; /* exitflag: when set, we exit (!) */

    /* static global arrays */

//    extern double *d; /* D's used for extrapolation in pzextr */
//    extern double *hpoints; /* stepsizes h (=H/n) which we have tried */

    /* allocate arrays */

    if (!(err = (double *) calloc(KMAXX, sizeof(double))))
        error("BuSt: error allocating err.\n");

    if (!(verror = (double *) calloc(n, sizeof(double))))
        error("BuSt: error allocating verror.\n");

    if (!(vsav = (double *) calloc(n, sizeof(double))))
        error("BuSt: error allocating vsav.\n");

    if (!(vseq = (double *) calloc(n, sizeof(double))))
        error("BuSt: error allocating vseq.\n");

    if (!(d = (double *) calloc(KMAXX * n, sizeof(double))))
        error("BuSt: error allocating d.\n");

    if (!(hpoints = (double *) calloc(KMAXX, sizeof(double))))
        error("BuSt: error allocating hpoints.\n");

    /* set initial stepsize to try to htry */

    h = htry;

    /* new accuracy? -> initialize (here, this only applies at start) */

    if (accuracy != old_accuracy) {

        *hnext = tnew = -1.0e29; /* init these to impossible value */

        /* allocate memory if necessary */

        if (!a) {
            a = (double *) calloc(IMAXX + 1, sizeof(double));
            alf = (double **) calloc(KMAXX + 1, sizeof(double *));
            for (i = 0; i < KMAXX + 1; i++)
                alf[i] = (double *) calloc(KMAXX + 1, sizeof(double));
        }

        /* initialize the work coefficients */

        a[1] = nseq[1] + 1;
        for (k = 1; k <= KMAXX; k++)
            a[k + 1] = a[k] + nseq[k + 1];

        /* initialize the correction factors (alpha) */

        accuracy1 = SAFE1 * accuracy; /* accuracy1 used to calculate alphas */
        for (i = 2; i <= KMAXX; i++)
            for (k = 1; k < i; k++)
                alf[k][i] = pow(
                        accuracy1,
                        ((a[k + 1] - a[i + 1]) / ((a[i + 1] - a[1] + 1.0)
                                * (2 * k + 1))));
        old_accuracy = accuracy;

        /* determine optimal row number for convergence */

        for (kopt = 2; kopt < KMAXX; kopt++)
            if (a[kopt + 1] > a[kopt] * alf[kopt - 1][kopt])
                break;
        kmax = kopt;

    }

    /* actual stepping starts here: first save starting values of v[] in vsav  */

    for (i = 0; i < n; i++)
        vsav[i] = v[i];

    /* new integration or stepsize? -> re-establish the order window */

    if (*t != tnew || h != (*hnext)) {
        first = 1;
        kopt = kmax;
    }

    reduct = 0; /* we have not reduced stepsize yet */

    for (;;) {

        /* try increasing numbers of steps over the interval h */

        for (k = 1; k <= kmax; k++) {

            tnew = (*t) + h;
            if (tnew == (*t))
                error("BuSt: stepsize underflow in bsstep\n");

            /* propagate the equations from t to t+h with nseq[k] steps using vsav and *
             * deriv as input; mmid returns vseq                                       */

            mmid(vsav, vseq, deriv, *t, h, nseq[k], n);

            /* extrapolate v's for h->0 based on the different stepsizes (hest's) we   *
             * have tried already; the result is returned in v (errors in verror);     *
             * hest is squared since the error series is even                          */

            hest = DSQR(h / nseq[k]);
            pzextr(k - 1, hest, vseq, v, verror, n, KMAXX);

            if (k > 1) {

                /* find the maximum error */

                verror_max = 0.;
                for (i = 0; i < n; i++)
                    if (v[i] != 0.)
                        verror_max = max(fabs(verror[i] / v[i]), verror_max);
                    else
                        verror_max = max(fabs(verror[i] / DBL_EPSILON),
                                         verror_max);

                /* scale error according to desired accuracy */

                verror_max /= accuracy;

                /* compute normalized error estimates (epsilons) */

                km = k - 1;
                err[km - 1] = pow(verror_max / SAFE1, 1.0 / (2 * km + 1));

            }

            /* are we in the order window? -> converged */

            if (k > 1 && (k >= kopt - 1 || first)) {

                if (verror_max < 1.0) { /* exit if accuracy good enough */
                    exitflag = 1;
                    break;
                }

                if (k == kmax || k == kopt + 1) { /* stepsize reduction possible? */
                    red = SAFE2 / err[km - 1];
                    break;
                } else if (k == kopt && alf[kopt - 1][kopt] < err[km - 1]) {
                    red = 1.0 / err[km - 1];
                    break;
                } else if (kopt == kmax && alf[km][kmax - 1] < err[km - 1]) {
                    red = alf[km][kmax - 1] * SAFE2 / err[km - 1];
                    break;
                } else if (alf[km][kopt] < err[km - 1]) {
                    red = alf[km][kopt - 1] / err[km - 1];
                    break;
                }
            }
        }

        if (exitflag)
            break;

        /* reduce stepsize by at least REDMIN and at most REDMAX, then try again */

        red = min(red, REDMIN);
        red = max(red, REDMAX);

        h *= red;
        reduct = 1;
    }

    /* we've taken a successful step */

    *t = tnew;
    *hdid = h;
    first = 0;
    wrkmin = 1.0e35;

    /* compute optimal row for convergence and corresponding stepsize */

    for (i = 1; i <= km; i++) {

        fact = max(err[i - 1], SCALMX);
        work = fact * a[i + 1];

        if (work < wrkmin) {

            scale = fact;
            wrkmin = work;
            kopt = i + 1;

        }
    }

    *hnext = h / scale;

    /* check for possible order increase, but not if stepsize was just reduced */

    if (kopt >= k && kopt != kmax && !reduct) {

        fact = max(scale / alf[kopt - 1][kopt], SCALMX);

        if (a[kopt + 1] * fact <= wrkmin) {
            *hnext = h / fact;
            kopt++;
        }
    }

    /* clean up */

    free(d);
    free(hpoints);
    free(vseq);
    free(verror);
    free(err);
    free(vsav);
}

/*** pzextr: uses polynomial extrapolation (Neville's algorithm) to evalu- *
 *           ate v's at a the hypothetical stepsize 0; this is called Ri-  *
 *           chardson extrapolation; the arguments are:                    *
 *           - iest: how many steps have we done already before?           *
 *           - hest: current stepsize h                                    *
 *           - vest: v's obtained using current stepsize h                 *
 *           - vout: extrapolated v's that will be returned                *
 *           - dv:   array of error estimates to be returned               *
 *           - n:    size of verst, vout and dv arrays                     *
 *           Neville's algorithm uses a recursive approach to determine a  *
 *           suitable Lagrange polynomial for extrapolation by traversing  *
 *           a tableau of differences between Lagrange polynomials of in-  *
 *           creasing order; these differences are called C and D below    *
 *           (see Numerical Recipes in C, Chapter 3.1 for details)         *
 ***************************************************************************/

void StepSolver::pzextr(int iest, double hest, double *vest, double *vout,
                        double *dv, int n, const int KMAXX)
{
    int i, j; /* loop counters */

    double q; /* temp vars for calculating C's and D's */
    double f1, f2;
    double delta;
    double *c; /* C's used for extrapolation */

//    extern double *d; /* D's used for extrapolation */
//    extern double *hpoints; /* stepsizes h (=H/n) which we have tried */

    if (!(c = (double *) calloc(n, sizeof(double))))
        error("pzextr: error allocating c.\n");

    hpoints[iest] = hest; /* stepsizes h (=H/n) which we have tried */

    for (i = 0; i < n; i++)
        dv[i] = vout[i] = vest[i];

    /* first time around? -> store first estimate in first column of d */

    if (iest == 0)
        for (i = 0; i < n; i++)
            d[i * KMAXX] = vest[i];

    /* more than one point? -> do polynomial extrapolation to h->0 */

    else {

        for (i = 0; i < n; i++)
            c[i] = vest[i];

        for (i = 1; i <= iest; i++) {

            /* calculate new values of C's and D's to traverse Neville's tableau */

            delta = 1.0 / (hpoints[iest - i] - hest);
            f1 = hest * delta;
            f2 = hpoints[iest - i] * delta;

            /* propagate tableau one diagonal more; the error is calculated using the  *
             * values of the D's */

            for (j = 0; j < n; j++) {
                q = d[j * KMAXX + (i - 1)];
                d[j * KMAXX + (i - 1)] = dv[j];
                delta = c[j] - q;
                dv[j] = f1 * delta;
                c[j] = f2 * delta;
                vout[j] += dv[j];

            }
        }

        /* save current D's for future calls to this function in the appropriate   *
         * column of the D tableau                                                 */

        for (i = 0; i < n; i++)
            d[i * KMAXX + iest] = dv[i];

    }

    free(c);
}

/*** mmid: implements the modified midpoint method used by BuSt; it subdi- *
 *         vides a large step (htot) into nstep intervals and uses a mid-  *
 *         point method over the whole of htot, with a stepsize of 2*h     *
 *         except for the first and the last step, hence *modified* mid-   *
 *         point method); the nice thing about this method is that the     *
 *         error follows a power law depending on the stepsize, i.e. its   *
 *         error converges to zero really fast as h is diminished; this    *
 *         allows for good extrapolation to h=0 (see pzextr() above)       *
 ***************************************************************************/

void BuSt::mmid(double *vin, double *vout, double *deriv, double tin,
                double htot, int nstep, int n)
{

    int i, j; /* loop counters */

    double *vm; /* v's at beginning of 2*h */
    double *vn; /* v's at midpoint of 2*h */
    double swap; /* tmp var used for swapping vm and vn */

    double t; /* current time */

    double h; /* small (h) stepsize (equals htot / nstep) */
    double h2; /* double the small (h) stepsize */

    vm = (double *)calloc(n, sizeof(double));
    vn = (double *)calloc(n, sizeof(double));

    /* calculate h from H and n */

    h = htot / nstep;

    /* calculate the first step (according to Euler) */

    for (i = 0; i < n; i++) {
        vm[i] = vin[i];
        vn[i] = vin[i] + h * deriv[i];
    }

    /* do the intermediate steps */

    h2 = 2.0 * h; /* mmid uses 2 * h as stepsize for interm. steps */
    t = tin + h;
    zygote.p_deriv(vn, t, vout, n);

    for (i = 2; i <= nstep; i++) {
        for (j = 0; j < n; j++) {
            swap = vm[j] + h2 * vout[j];
            vm[j] = vn[j];
            vn[j] = swap;
        }

        t += h;
        zygote.p_deriv(vn, t, vout, n);
    }

    /* calculate the last step */

    for (i = 0; i < n; i++)
        vout[i] = 0.5 * (vm[i] + vn[i] + h * vout[i]);

    /* clean up */

    free(vm);
    free(vn);
}

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

void BaDe::ps(double *vin, double *vout, double tin, double tout,
              double stephint, double accuracy, int n, FILE *slog)
{

    int i; /* local loop counter */

    double t; /* current time */

    double ht; /* stepsize of next step */
    double hd; /* stepsize of step we just did */

    double *initderiv; /* array for initial derivatives */

    /* the do-nothing case; too small steps dealt with under usual */

    if (tin == tout)
        return;

    /* the usual case: steps big enough -> initialize variables */

    initderiv = (double *) calloc(n, sizeof(double));

    t = tin;

    /* initial stepsize is either stephint or tout-tin */

    ht = min(stephint, (tout - tin));

    /* we're saving old v's in vin and propagate vout below */

    for (i = 0; i < n; i++)
        vout[i] = vin[i];

    /* this loop steps through the whole interval tout-tin */

    do {

        zygote.p_deriv(vout, t, initderiv, n);
        stifbs(vout, initderiv, n, &t, ht, accuracy, &hd, &ht);

        if (t + ht > tout)
            ht = tout - t;

    } while (t < tout);

    /* clean up and go home */

    free(initderiv);

    return;
}

/*** stifbs: this is the Bader-Deuflhard stepper function; it propagates ***
 *           the equations at a given overall stepsize and order, then     *
 *           evaluates the error and repeats the step with increased order *
 *           if necessary, until the required accuracy is achieved; its    *
 *           arguments are:                                                *
 *           - v:        input and output v's over a give total stepsize   *
 *           - deriv:    array of initial derivatives (at tin)             *
 *           - n:        size of v and deriv                               *
 *           - t:        pointer to the current time                       *
 *           - htry:     initial stepsize to be tried                      *
 *           - accuracy: relative accuracy to be achieved                  *
 *           - hdid:     returns the stepsize of the step we just did      *
 *           - hnext:    returns the stepsize of the next step to be done  *
 ***************************************************************************/
void BaDe::stifbs(double *v, double *deriv, int n, double *t, double htry,
                  double accuracy, double *hdid, double *hnext)
{

    /*** constants *************************************************************/

    const int KMAXX = 7; /* this is how many successive stepsize */
    const int IMAXX = KMAXX + 1; /* refinements we're going to try */

    const double SAFE1 = 0.25; /* safety margins for errors and h reduction */
    const double SAFE2 = 0.7;

    const double REDMAX = 1.0e-5; /* boundaries for stepsize reduction */
    const double REDMIN = 0.7;

    const double SCALMX = 0.1; /* maximum scaling value for stepsize */

    /*** variables *************************************************************/

    int i, k, km; /* loop counters */

    double *vsav; /* v's at beginning of interval */
    double *vseq; /* v's at end of interval for a given stepsize */

    double *verror; /* error estimates */
    double verror_max; /* the maximum error in verror */
    double *err; /* normalized error estimates */

    double *dfdt;
    double **jac; /* the Jacobian matrix */

    double h; /* stepsize for solver method used below */
    double hest; /* square of h passed to extrapolation function */
    double red; /* reduction factor for reducing stepsize */
    double scale; /* scaling factor for next stepsize */

    double work; /* temp vars for calculating row for convergence */
    double wrkmin;
    double fact;

    static double old_accuracy = -1.0; /* used to save old accuracy */
    static double tnew; /* used to save old start time */
    static int nold = -1; /* for saving old value of n */

    /* the following two arrays are used for Deuflhard's error estimation; a   *
     * contains the work coefficients and alf (alpha) the correction factors   */

    static double *a = NULL;
    static double **alf = NULL;

    double accuracy1; /* error (< accuracy) used to calculate alphas */

    static int kmax; /* used for finding kopt */
    static int kopt; /* optimal row number for convergence */

    /* sequence of separate attempts to cross interval htot with increasing    *
     * values of nsteps as suggested by Deuflhard (Num Rec, p. 726)            */

    static int nseq[] = { 0, 2, 6, 10, 14, 22, 34, 50, 70 };

    /* miscellaneous flags */

    static int first = 1; /* is this the first try for a given step? */
    int reduct; /* flag indicating if we have reduced stepsize yet */
    int exitflag = 0; /* exitflag: when set, we exit (!) */

    /* static global arrays */

//    extern double *d; /* D's used for extrapolation in pzextr */
//    extern double *hpoints; /* stepsizes h (=H/n) which we have tried */

    /* allocate arrays */

    if (!(d = (double *) calloc(KMAXX * n, sizeof(double))))
        error("BaDe: error allocating d.\n");

    if (!(dfdt = (double *) calloc(n, sizeof(double))))
        error("BaDe: error allocating dfdt.\n");

    if (!(jac = (double **) calloc(n, sizeof(double *))))
        error("BaDe: error allocating jac.\n");

    for (i = 0; i < n; i++)
        if (!(jac[i] = (double *) calloc(n, sizeof(double))))
            error("BaDe: error allocating jac's second dimension.\n");

    if (!(err = (double *) calloc(KMAXX, sizeof(double))))
        error("BaDe: error allocating err.\n");

    if (!(verror = (double *) calloc(n, sizeof(double))))
        error("BaDe: error allocating verror.\n");

    if (!(vsav = (double *) calloc(n, sizeof(double))))
        error("BaDe: error allocating vsav.\n");

    if (!(vseq = (double *) calloc(n, sizeof(double))))
        error("BaDe: error allocating vseq.\n");

    if (!(hpoints = (double *) calloc(KMAXX, sizeof(double))))
        error("BaDe: error allocating hpoints.\n");

    /* set initial stepsize to try to htry */

    h = htry;

    /* allocate memory if necessary */

    if (!a) {
        a = (double *) calloc(IMAXX + 1, sizeof(double));
        alf = (double **) calloc(KMAXX + 1, sizeof(double *));
        for (i = 0; i < KMAXX + 1; i++)
            alf[i] = (double *) calloc(KMAXX + 1, sizeof(double));
    }

    /* new accuracy? -> initialize (here, this only applies at start) */

    if ((accuracy != old_accuracy) || (n != nold)) {

        *hnext = tnew = -1.0e29; /* init these to impossible value */

        /* initialize the work coefficients */

        a[1] = nseq[1] + 1;
        for (k = 1; k <= KMAXX; k++)
            a[k + 1] = a[k] + nseq[k + 1];

        /* initialize the correction factors (alpha) */

        accuracy1 = SAFE1 * accuracy; /* accuracy1 used to calculate alphas */
        for (i = 2; i <= KMAXX; i++)
            for (k = 1; k < i; k++)
                alf[k][i] = pow(
                        accuracy1,
                        ((a[k + 1] - a[i + 1]) / ((a[i + 1] - a[1] + 1.0)
                                * (2 * k + 1))));
        old_accuracy = accuracy;
        nold = n;

        /* add cost of Jacobian evaluations to work coefficients */

        a[1] += n;
        for (k = 1; k <= KMAXX; k++)
            a[k + 1] = a[k] + nseq[k + 1];

        /* determine optimal row number for convergence */

        for (kopt = 2; kopt < KMAXX; kopt++)
            if (a[kopt + 1] > a[kopt] * alf[kopt - 1][kopt])
                break;
        kmax = kopt;

    }

    /* actual stepping starts here: first save starting values of v[] in vsav  */

    for (i = 0; i < n; i++)
        vsav[i] = v[i];

    /* evaluate jacobian */

    zygote.p_jacobn(*t, v, dfdt, jac, n);

    /* new integration or stepsize? -> re-establish the order window */

    if (*t != tnew || h != (*hnext)) {
        first = 1;
        kopt = kmax;
    }

    reduct = 0; /* we have not reduced stepsize yet */

    for (;;) {

        /* try increasing numbers of steps over the interval h */

        for (k = 1; k <= kmax; k++) {

            tnew = (*t) + h;
            if (tnew == (*t))
                error("BaDe: stepsize underflow in stifbs.\n");

            /* propagate the equations from t to t+h with nseq[k] steps using vsav and *
             * deriv as input; dfdt are the function derivs and jac is the Jacobian    *
             * trix; simpr returns vseq                                                */

            simpr(vsav, vseq, deriv, dfdt, jac, *t, h, nseq[k], n);

            /* extrapolate v's for h->0 based on the different stepsizes (hest's) we   *
             * have tried already; the result is returned in v (errors in verror);     *
             * hest is squared since the error series is even                          */

            hest = DSQR(h / nseq[k]);
            pzextr(k - 1, hest, vseq, v, verror, n, KMAXX);

            if (k > 1) {

                /* find the maximum error */

                verror_max = 0.;
                for (i = 0; i < n; i++)
                    if (v[i] != 0.)
                        verror_max = max(fabs(verror[i] / v[i]), verror_max);
                    else
                        verror_max = max(fabs(verror[i] / DBL_EPSILON),
                                         verror_max);

                /* scale error according to desired accuracy */

                verror_max /= accuracy;

                /* compute normalized error estimates (epsilons) */

                km = k - 1;
                err[km - 1] = pow(verror_max / SAFE1, 1.0 / (2 * km + 1));

            }

            /* are we in the order window? -> converged */

            if (k > 1 && (k >= kopt - 1 || first)) {

                if (verror_max < 1.0) { /* exit if accuracy good enough */
                    exitflag = 1;
                    break;
                }

                if (k == kmax || k == kopt + 1) { /* stepsize reduction possible? */
                    red = SAFE2 / err[km - 1];
                    break;
                } else if (k == kopt && alf[kopt - 1][kopt] < err[km - 1]) {
                    red = 1.0 / err[km - 1];
                    break;
                } else if (kopt == kmax && alf[km][kmax - 1] < err[km - 1]) {
                    red = alf[km][kmax - 1] * SAFE2 / err[km - 1];
                    break;
                } else if (alf[km][kopt] < err[km - 1]) {
                    red = alf[km][kopt - 1] / err[km - 1];
                    break;
                }
            }
        }

        if (exitflag)
            break;

        /* reduce stepsize by at least REDMIN and at most REDMAX, then try again */

        red = min(red, REDMIN);
        red = max(red, REDMAX);

        h *= red;
        reduct = 1;
    }

    /* we've taken a successful step */

    *t = tnew;
    *hdid = h;
    first = 0;
    wrkmin = 1.0e35;

    /* compute optimal row for convergence and corresponding stepsize */

    for (i = 1; i <= km; i++) {

        fact = max(err[i - 1], SCALMX);
        work = fact * a[i + 1];

        if (work < wrkmin) {

            scale = fact;
            wrkmin = work;
            kopt = i + 1;

        }
    }

    *hnext = h / scale;

    /* check for possible order increase, but not if stepsize was just reduced */

    if (kopt >= k && kopt != kmax && !reduct) {

        fact = max(scale / alf[kopt - 1][kopt], SCALMX);

        if (a[kopt + 1] * fact <= wrkmin) {
            *hnext = h / fact;
            kopt++;
        }
    }

    /* clean up */

    for (i = 0; i < n; i++)
        free(jac[i]);
    free(jac);

    free(d);
    free(dfdt);
    free(hpoints);
    free(vseq);
    free(verror);
    free(err);
    free(vsav);
}

/*** simpr: implements the semi-implicit midpoint rule used by BaDe; it ****
 *          subdivides a large step (htot) into nstep intervals and uses a *
 *          semi-implicit midpoint rule over the whole of htot, with a     *
 *          stepsize of 2*h except for the first and the last step; the    *
 *          nice thing about this method is that the error follows a power *
 *          law depending on the stepsize, i.e. its error converges to 0   *
 *          really fast as h is diminished; this allows for good extrapo-  *
 *          lation to h=0 (see pzextr() above)                             *
 ***************************************************************************/

void BaDe::simpr(double *vin, double *vout, double *deriv, double *dfdt,
                 double **jac, double tin, double htot, int nstep, int n)
{

    /*** variables *************************************************************/

    int i, j; /* loop counters */

    double t; /* current time */
    double h; /* small (h) stepsize (equals htot / nstep) */
    double d; /* d indicates if eve/odd rows were switched in ludcmp */

    int *indx; /* array for permutated row indices of matrix a */

    double *del; /* delta: difference in v between steps */
    double *vtemp; /* array for temp derivs and v's */

    double **a; /* matrix [1 - hf'] */

    /* allocate arrays */

    indx = (int *) calloc(n, sizeof(int));

    del = (double *) calloc(n, sizeof(double));
    vtemp = (double *) calloc(n, sizeof(double));

    a = (double **) calloc(n, sizeof(double *));
    for (i = 0; i < n; i++)
        a[i] = (double *) calloc(n, sizeof(double));

    /* calculate h from H and n */

    h = htot / nstep;

    /* set up the matrix [1 - hf'] */

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            a[i][j] = -h * jac[i][j];
        ++a[i][i];
    }

    /* the following is slightly bizarre */

    /* do LU decomposition of matrix [1 - hf']; this is needed for all steps   *
     * below, which use linearization of the equations to get estimates of the *
     * derivatives at the endpoint of a step (which is done in lubksb)         */

    ludcmp(a, n, indx, &d);

    /* do the first step */

    /* set up the right-hand side for the first step; use vout for temp sto-   *
     * rage; note that since our equations are autonomous, all dfdt's are 0.0; *
     * lubksb solves the linear system given by a and vout and then returns    *
     * the solution vector in vout; this vector contains the difference be-    *
     * tween this step's v's and the next steps' v's which are stored in del   */

    for (i = 0; i < n; i++)
        vout[i] = h * (deriv[i] + h * dfdt[i]);
    lubksb(a, n, indx, vout);

    for (i = 0; i < n; i++)
        vtemp[i] = vin[i] + (del[i] = vout[i]);

    /* do the intermediate steps */

    t = tin + h;
    zygote.p_deriv(vtemp, t, vout, n);

    for (i = 2; i <= nstep; i++) {

        /* set up the right hand side */

        for (j = 0; j < n; j++)
            vout[j] = h * vout[j] - del[j];
        lubksb(a, n, indx, vout);

        /* take the step */

        for (j = 0; j < n; j++)
            vtemp[j] += (del[j] += 2.0 * vout[j]);

        /* go to next step */

        t += h;
        zygote.p_deriv(vtemp, t, vout, n);
    }

    /* do the last step */

    /* set up the right-hand side for the last step */

    for (i = 0; i < n; i++)
        vout[i] = h * vout[i] - del[i];
    lubksb(a, n, indx, vout);

    /* take the last step */

    for (i = 0; i < n; i++)
        vout[i] += vtemp[i];

    /* clean up */

    for (i = 0; i < n; i++)
        free(a[i]);
    free(a);

    free(vtemp);
    free(del);
    free(indx);

}

/*** ludcmp: does an LU decomposition of matrix a, which is needed for in- *
 *           verting the matrix; LU decomposition dissects the matrix into *
 *           a lower and an upper triangular matrix that when multiplied,  *
 *           equal matrix a again; this function uses Crout's algorithm    *
 *           for the decomposition; the result is returned in a and con-   *
 *           tains the LU decomposition of a rowwise permutation of a;     *
 *           indx returns the order of the permutated rows; n is the size  *
 *           of matrix a (as in n x n); d returns +1 if number of row in-  *
 *           terchanges was even, -1 if odd.                               *
 ***************************************************************************/

void BaDe::ludcmp(double **a, int n, int *indx, double *d) {
    int i, j, k; /* loop counters */
    int imax; /* index of row of largest element in a column */

    double big; /* largest value, used for pivoting */
    double dum; /* dummy used for pivoting */
    double sum; /* temp var for summing */
    double temp; /* temp var for absolute values */
    double *vv; /* stores the implicit scaling information */

    const double TINY = 1.0e-20; /* this is a tiny number */

    vv = (double *) calloc(n, sizeof(double));

    *d = 1.0; /* no interchanges yet */

    /* loop over rows to get the implicit scaling information */

    for (i = 0; i < n; i++) {

        big = 0.0;
        for (j = 0; j < n; j++)
            if ((temp = fabs(a[i][j])) > big)
                big = temp;
        /* whole row == 0 -> singular matrix */
        if (big == 0.0)
            error("ludcmp: cannot decompose singular matrix");

        vv[i] = 1.0 / big; /* save the scaling */
    }

    /* loop over columns of Crout's method; this is done in two loops which    *
     * reflect the two different formulas of Crout's method, except that the   *
     * first formula is used for the diagonal elements in the second loop; in  *
     * the second loop we're also finding the largest element and its row num- *
     * ber for pivoting                                                        */

    for (j = 0; j < n; j++) {

        for (i = 0; i < j; i++) {
            sum = a[i][j];
            for (k = 0; k < i; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }

        big = 0.0; /* used for search of largest pivot element */

        for (i = j; i < n; i++) {

            sum = a[i][j]; /* summing happens here */
            for (k = 0; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;

            if ((dum = vv[i] * fabs(sum)) >= big) { /* pivoting stuff here */
                big = dum;
                imax = i;
            }

        }

        /* pivoting: do we need to interchange rows? if yes -> do so! */

        if (j != imax) {
            for (k = 0; k < n; k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }

            *d = -(*d); /* change parity of d */
            vv[imax] = vv[j]; /* rearrange scale factors */
        }

        indx[j] = imax; /* store index of permutated row */

        /* if the pivot element is zero, the matrix is singular; TINY avoids div/0 *
         * below (some applications can deal with singular matrices, we don't)     */

        if (a[j][j] == 0.0)
            a[j][j] = TINY;

        /* divide by the pivot element */

        if (j != n - 1) {
            dum = 1.0 / (a[j][j]);
            for (i = j + 1; i < n; i++)
                a[i][j] *= dum;
        } /* end of loop: next column */

    }

    free(vv);
}
/*** lubksb: does forward and backsubstitution of matrix a; in fact, a is **
 *           not passed in its original form but as the LU decomposition   *
 *           of a rowwise permutation of a as returned by the function     *
 *           ludcmp(); the right hand side vector is called b, which also  *
 *           returns the solution vector to the calling function; indx is  *
 *           a vector that indicates the order of rows in the permutated   *
 *           matrix; n is the dimension of the matrix a (as in n x n).     *
 ***************************************************************************/

void BaDe::lubksb(double **a, int n, int *indx, double *b) {
    int i, j; /* counter variables */
    int ii = -1; /* index of first non-zero element of b */
    int ip; /* index of permutated matrix a */

    double sum; /* temp var for summing */

    /* first loop does the forward substitution; we don't worry about the dia- *
     * gonal elements of the a matrix since here they're all 1 by definition   *
     * (see description of Crout's algorithm in Num Rec, chapter 2.3)          */

    for (i = 0; i < n; i++) {

        ip = indx[i]; /* get the index of the (real) first row */
        sum = b[ip]; /* get the according b */
        b[ip] = b[i]; /* un-permutate the b-vector */
        /* start summing at first non-zero element of b */
        if (ii >= 0)
            for (j = ii; j <= i - 1; j++)
                sum -= a[i][j] * b[j];
        else if (sum)
            ii = i;

        b[i] = sum; /* b's are now in the right (un-permutated) order */

    }

    /* second loop does the backsubstitution */

    for (i = n - 1; i >= 0; i--) {
        sum = b[i];
        for (j = i + 1; j <= n - 1; j++)
            sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

SoDe::SoDe(zygotic &in_zy, int in_debug) :
        solver(in_zy, in_debug), TheMaternal(zygote.get_Maternal()),
        defs(TheMaternal.getProblem())
{
    gridpos = -1;
    gridstart = 0;
    tdone = NULL;
    vdonne = NULL;
    derivv1 = NULL;
    derivv2 = NULL;
    derivv3 = NULL;
    derivv4 = NULL;
    lp = zygote.GetMutParameters();
    maxdel = 0.;
    mindel = 1000.;

    for (int j = 0; j < defs.ngenes; j++) {

        if (lp->tau[j] > maxdel)
            maxdel = lp->tau[j];
        if (lp->tau[j] < mindel)
            mindel = lp->tau[j];

    }
    numdel = defs.ngenes;
    delay = lp->tau;
}

SoDe::~SoDe()
{
    int j;
    free(fact_discons);

    for (j=0; j<=gridpos;j++)
    {
        if (derivv1 && derivv1[j])
            free(derivv1[j]);
        if (derivv2 && derivv2[j])
            free(derivv2[j]);
        if (derivv3 && derivv3[j])
            free(derivv3[j]);
        if (derivv4 && derivv4[j])
            free(derivv4[j]);
        if (vdonne && vdonne[j])
            free(vdonne[j]);
    }

    if (derivv1)
        free(derivv1);
    if (derivv2)
        free(derivv2);
    if (derivv3)
        free(derivv3);
    if (derivv4)
        free(derivv4);
    if (vdonne)
        free(vdonne);
    if (tdone)
        free(tdone);
    derivv1 = NULL;
    derivv2 = NULL;
    derivv3 = NULL;
    derivv4 = NULL;
    tdone = NULL;
    vdonne = NULL;
}

void SoDe::resetSolver()
{
    int j;
    // This is not thread safe!! Do more check before port to OpenMP
    for (j=0; j<=gridpos;j++)
    {
        free(derivv1[j]);
        free(derivv2[j]);
        free(derivv3[j]);
        free(derivv4[j]);
        free(vdonne[j]);
    }

    free(derivv1);
    free(derivv2);
    free(derivv3);
    free(derivv4);
    free(vdonne);
    free(tdone);

    derivv1 = NULL;
    derivv2 = NULL;
    derivv3 = NULL;
    derivv4 = NULL;
    tdone = NULL;
    vdonne = NULL;
}

void SoDe::ps(double *vin, double *vout, double tin, double tout,
              double stephint, double accuracy, int n, FILE *slog)
{

    int i, j;


    double *shift_discs;
    double *discon_array;
    int num_discons;

    double **vees, *tees;


    /*  printf("%d %f %f\n",fact_discons_size,fact_discons[0],fact_discons[2]);*/



    /*  printf("maxdel:%f, mindel:%f and numdel:%d\n",
     maxdel,mindel,numdel);*/

    j = -1;
    do {
        j++;

    } while ((j < fact_discons_size) && (tout > fact_discons[j]));

    shift_discs = (double *) calloc(j + 1, sizeof(double));
    shift_discs[0] = 0.;

    /*  printf("The %d discons for [%f,%f] are %f %f\n",j,tin,tout,tin,shift_discs[0]);*/

    for (i = 0; i < j; i++) {

        shift_discs[i + 1] = fact_discons[i] - tin;
        /*      printf("The discons for [%f,%f] are %.16f %f\n",tin,tout,fact_discons[i], shift_discs[i+1]);*/

    }

    discon_array = Construct_Discont_Array(tout - tin, lp->tau,
                                           defs.ngenes, shift_discs,
                                           j + 1, &num_discons);

    for (i = 0; i < num_discons; i++) {
        discon_array[i] += tin;
        /*      printf("Disconarray[%d] in [%f,%f]=%f\n",i,tin,tout,discon_array[i]);*/
    }

    /*  printf("Size_Discon_Array:%d\n",num_discons);*/

    vees = (double **) calloc(2, sizeof(double *));
    tees = (double *) calloc(2, sizeof(double));

    tees[0] = tin;
    tees[1] = tout;
    vees[0] = vin;
    vees[1] = vout;

    /*  printf("nput and output:\n");

     for (i=0; i < n; i++)
     printf("v[%d](%f)=%f,%d,%d\n",i,tees[0],vees[0][i],vees[0],vees[1]);
     */

    DCERk32(vees, n, tees, 2, discon_array, num_discons, stephint, accuracy);

    /*  for (i=0; i < n; i++)
     printf("v[%d](%f)=%f\n",i,tees[1],vees[1][i]);*/

    free(shift_discs);
    if (discon_array)
        free(discon_array);
    free(vees);
    free(tees);

    return;

}

/*** DCERk3(2): propagates v[0] (of size n) according to tarray  by  **
 the Runge-Kutta 3(2) pair with continuous extension storing the     **
 result in vatt. Initial conditions are specified in vatt[0],        **
 corresponding to tarray[0]                                              **
 *********************************************************************/

void SoDe::DCERk32(double **vatt, int n, double *tarray, int tpoints,
                   double *darray, int dpoints, double stephint,
                   double accuracy)
{
    int i, j, dc, vc; /* local loop counters */

    int tpos = 1; /* where in the t array we are, starting
     at the value right after the starting time*/

    int dpos = 0; /* where is the discontinuities array we
     are (darray) */
    int DISCON = 0; /* if we just hit a discontinuity */

    double **v; /** used for storing intermediate steps */
    int toggle = 0; /* used to toggle between v[0] and v[1] */

    double *vtemp; /* guessed intermediate v's */

    double *vnow; /* ptr to v at current time */
    double *vnext; /* ptr to v at current time + stepsize */
    double **v_delayed[4]; /* ptr to v at t-delay */

    double *verror; /* error estimate */
    double verror_max; /* the maximum error in verror */

    double *tempswap;

    double t; /* the current time */
    double tms[4]; /* array for the rk step times */
    double h = stephint; /* initial stepsize */
    double hnext; /* used to calculate next stepsize */
    const double SAFETY = 0.9; /* safety margin for decrsg stepsize */
    const double ERRCON = 5.832e-3; /* to prevent huge jumps, see
     num recipes pg. 719 */

    /* declare and initialize parameters */
    /* parameters that are zero have been added as comments for clarity */

    static double step_exponent = -1.0 / 3.0;

    static double
    /*  a1  =      0.0, */
    a2 = 0.5, a3 = 0.75;
    /*  a4  =     1.0, */

    static double b21 = 0.5,
    /*  b31 =     0.0, */
    b32 = 0.75;

    static double c1 = 2.0 / 9.0, c2 = 1.0 / 3.0, c3 = 4.0 / 9.0;

    double dc1 = c1 - 7.0 / 24.0, dc2 = c2 - 0.25, dc3 = c3 - c2, dc4 = -0.125;

    /* the do-nothing case; too small steps dealt with under usual */

    if (tarray[0] >= tarray[tpoints - 1])
        return;

    /* the usual case: steps big enough */

    v = (double **) calloc(2, sizeof(double *));

    vtemp = (double *) calloc(n, sizeof(double));
    verror = (double *) calloc(n, sizeof(double));
    v[0] = (double *) calloc(n, sizeof(double));
    v[1] = (double *) calloc(n, sizeof(double));

    for (vc = 0; vc < 4; vc++) {
        v_delayed[vc] = (double **) calloc(numdel, sizeof(double *));
        for (dc = 0; dc < numdel; dc++)
            v_delayed[vc][dc] = (double *) calloc(n, sizeof(double));
    }

    t = tarray[0];
    vnow = vatt[0];
    vnext = v[0];

    /* allocate more grid and derivs */

    gridpos++;
    if (!(tdone = (double *) realloc(tdone, (gridpos + 1) * sizeof(double)))) {
        printf("Could not allocate memory for tdone\n");
        exit(1);
    }

    if (!(vdonne = (double **) realloc(vdonne, (gridpos + 1) * sizeof(double *)))) {
        printf("COuld not allocate memory for vdonne\n");
        exit(1);
    }

    if (!(derivv1 = (double **) realloc(derivv1,
                                        (gridpos + 1) * sizeof(double *)))) {
        printf("COuld not allocate memory for derivv1\n");
        exit(1);
    }

    if (!(derivv2 = (double **) realloc(derivv2,
                                        (gridpos + 1) * sizeof(double *)))) {
        printf("COuld not allocate memory for derivv2\n");
        exit(1);
    }

    if (!(derivv3 = (double **) realloc(derivv3,
                                        (gridpos + 1) * sizeof(double *)))) {
        printf("COuld not allocate memory for derivv3 \n");
        exit(1);
    }

    if (!(derivv4 = (double **) realloc(derivv4,
                                        (gridpos + 1) * sizeof(double *)))) {
        printf("COuld not allocate memory for derivv4\n");
        exit(1);
    }

    derivv1[gridpos] = (double *) calloc(n, sizeof(double));
    derivv2[gridpos] = (double *) calloc(n, sizeof(double));
    derivv3[gridpos] = (double *) calloc(n, sizeof(double));
    derivv4[gridpos] = (double *) calloc(n, sizeof(double));
    vdonne[gridpos] = (double *) calloc(n, sizeof(double));

    tdone[gridpos] = t;
    memcpy(vdonne[gridpos], vnow, sizeof(*vnow) * n);

    /* initial stepsize cannot be bigger than total time */

    if (tarray[0] + h >= tarray[tpoints - 1])
        h = tarray[tpoints - 1] - tarray[0];

    /*  printf("%f\n",step_exponent);*/

    /*  printf("disconny:%d %d\n",dpoints,dpos);*/

    if ((dpoints) && (dpos < dpoints))
        if ((darray[dpos] > t) && (darray[dpos] <= t + h)) {

            /*          printf("\n***Hit discont at %f ***\n",darray[dpos]); */

            h = darray[dpos] - t;
            dpos++;
            /*          printf("dpos:%d\n",dpos);*/
            DISCON = 1;
        }

    if ((h > mindel) && (h <= 2. * mindel)) {

        if (DISCON) {

            dpos--;
            /*          printf("dpos:%d\n",dpos);*/
            DISCON = 0;

        }
        h = .5 * h;
    }

    tms[0] = t;
    tms[1] = t + a2 * h;
    tms[2] = t + a3 * h;
    tms[3] = t + h;

    /* we need to calculate derivv1 only the first time, since if the
     previous step was a success, we can use the last derivv4, and if it is
     was a failure, we don't have to recalculate it */

    while (y_delayed(v_delayed, n, tms, delay, tdone, vdonne, derivv1, derivv2,
                     derivv3, derivv4, gridpos + 1, accuracy)) {
        /*          printf("Rejected Iteration for [%f,%f]!\n",t,t+h);*/
        h = 0.5 * h;

        tms[0] = t;
        tms[1] = t + a2 * h;
        tms[2] = t + a3 * h;
        tms[3] = t + h;

        if (DISCON) {

            dpos--;
            /*              printf("dpos:%d\n",dpos);*/
            DISCON = 0;

        }
    }
    zygote.d_deriv(vnow, v_delayed[0], t, derivv1[gridpos], n);

    while (t < tarray[tpoints - 1]) {

        /* Take one step and evaluate the error. Repeat until the resulting error  *
         * is less than the desired accuracy                                       */

        while (1) {

            tms[0] = t;
            tms[1] = t + a2 * h;
            tms[2] = t + a3 * h;
            tms[3] = t + h;

            if (!y_delayed(v_delayed, n, tms, delay, tdone, vdonne, derivv1,
                           derivv2, derivv3, derivv4, gridpos + 1, accuracy))

                           {

                /* do the Runge-Kutta thing here: calulate intermediate
                 derivatives */

                for (i = 0; i < n; i++)
                    vtemp[i] = vnow[i] + h * (b21 * derivv1[gridpos][i]);

                zygote.d_deriv(vtemp, v_delayed[1], t + a2 * h,
                               derivv2[gridpos], n);

                for (i = 0; i < n; i++)
                    vtemp[i] = vnow[i] + h * (b32 * derivv2[gridpos][i]);

                zygote.d_deriv(vtemp, v_delayed[2], t + a3 * h,
                               derivv3[gridpos], n);

                /* ... then feed them to the Rk32 formula */

                for (i = 0; i < n; i++)
                    vnext[i] = vnow[i]
                            + h * (c1 * derivv1[gridpos][i] + c2
                                    * derivv2[gridpos][i]
                                   + c3 * derivv3[gridpos][i]);

                /* Now calculate k4 for the embedded 4-stage formula, if this step is
                 succesful, it will get used as derivv1 (k1) in the next step */

                zygote.d_deriv(vnext, v_delayed[3], t + h, derivv4[gridpos], n);

                /* calculate the error estimate using the embedded formula */

                for (i = 0; i < n; i++)
                    verror[i] = h
                            * (dc1 * derivv1[gridpos][i] + dc2
                                    * derivv2[gridpos][i]
                               + dc3 * derivv3[gridpos][i]
                               + dc4 * derivv4[gridpos][i]);

                /* find the maximum error */

                verror_max = 0.;
                for (i = 0; i < n; i++)
                    if (fabs(vnext[i]) > 10000. * DBL_EPSILON)
                        verror_max = max(fabs(verror[i] / vnext[i]),
                                         verror_max);
                    else
                        verror_max = max(
                                fabs(verror[i] / (10000. * DBL_EPSILON)),
                                verror_max);

                /*      for (i=0; i<n; i++)
                 printf("The deriv1[%d]=%.10f deriv2[%d]=%.10f deriv3[%d]=%.10f "
                 "deriv4[%d]=%.10f verror[%d]=%.10f verror_max=%.10f \n",i,
                 derivv1[gridpos][i],i, derivv2[gridpos][i],i,derivv3[gridpos][i],i,derivv4[gridpos][i],i,
                 verror[i],verror_max);
                 */
                /* scale error according to desired accuracy */

                verror_max /= accuracy;

                /* compare maximum error to the desired accuracy; if error < accuracy, we  *
                 * are done with this step; otherwise, the stepsize has to be reduced and  *
                 * the step repeated; for detailed comments on the approximation involving *
                 * SAFETY and -0.25 see 'Numerical Recipes in C', 2nd Edition, p.718       */

                if (verror_max <= 1.0)
                    break;

                /*    printf("Step size %f rejected!\n",h);*/

                /* kludge for going back to the original position in the discontinuity
                 array if the step is rejected */
                if (DISCON) {

                    dpos--;
                    /*      printf("dpos:%d\n",dpos);*/
                    DISCON = 0;

                }

                hnext = SAFETY * h * pow(verror_max, step_exponent);

                /* decrease stepsize by no more than a factor of 10; check for underflows */

                h = (hnext > 0.1 * h) ? hnext : 0.1 * h;

                if ((h > mindel) && (h <= 2. * mindel)) {

                    if (DISCON) {

                        dpos--;
                        /*          printf("dpos:%d\n",dpos);*/
                        DISCON = 0;

                    }
                    h = .5 * h;
                }

                /*    printf("Step size %f suggested!\n",h);*/

                if (h == 0)
                    error("Rkck: stepsize underflow");

            } else {

                /*      printf("Rejected Iteration for [%f,%f]!\n",t,t+h);*/
                h = 0.5 * h;

                if ((h > mindel) && (h <= 2. * mindel))
                    h = .5 * h;

                if (DISCON) {

                    dpos--;
                    /*      printf("dpos:%d\n",dpos);*/
                    DISCON = 0;

                }
            }

        }

        /* advance the current time by last stepsize */

        /*  printf("The stepsize at time %f is %f error is %f\n",
         t,h,verror_max*accuracy);*/

        /* Now we are going to do the continuous extension business, we will
         look for t < t* <= t+h, and for all such t*, we will calculate the CE
         and add the entried to the the vatt array. If t* = t+h, we just
         exchange the pointers */

        while ((tarray[tpos] < t + h) && (tpos < tpoints)) {

            CE(tarray[tpos], vatt[tpos], t, vnow, h, derivv1[gridpos],
               derivv2[gridpos], derivv3[gridpos], derivv4[gridpos], n);
            /*      printf("Vatt: %d %f %f %f %f %f\n",
             tpos,tarray[tpos],vatt[tpos][0],vnow[0],t,t+h); */

            tpos += 1;
        }

        if ((tarray[tpos] == t + h) && (tpos < tpoints)) {

            memcpy(vatt[tpos], vnext, sizeof(double) * n);
            tpos += 1;
        }

        /* toggle vnow and vnext between v[0] and v[1], vnow is just
         vnow_temp right now, as we need the old vnow in the Continuous
         Extension stuff*/

        vnow = v[toggle];
        toggle++;
        toggle %= 2;
        vnext = v[toggle];

        t += h;

        /*    if (t >= tarray[tpoints - 1])
         break;                                that was the last iteration */

        /* increase stepsize according to error (5th order) for next iteration */

        /*  printf("ERRCON:%f\n",ERRCON);*/

        if (verror_max > ERRCON) {

            h = SAFETY * h * pow(verror_max, step_exponent);

        } else
            h = 5.0 * h;

        /* make sure t does not overstep last time point */

        if (t + h >= tarray[tpoints - 1])
            h = tarray[tpoints - 1] - t;

        DISCON = 0;

        if ((dpoints) && (dpos < dpoints))
            if ((darray[dpos] > t) && (darray[dpos] <= t + h)) {

                /*          printf("\n***Hit discont at %f ***\n",darray[dpos]); */

                h = darray[dpos] - t;
                dpos++;
                /*          printf("dpos:%d\n",dpos);*/
                DISCON = 1;
            }

        if ((h > mindel) && (h <= 2. * mindel)) {

            if (DISCON) {

                dpos--;
                /*          printf("dpos:%d\n",dpos);*/
                DISCON = 0;

            }
            h = .5 * h;
        }

        /* allocate more grid and derivs */

        gridpos++;
        if (!(tdone = (double *) realloc(tdone, (gridpos + 1) * sizeof(double)))) {
            printf("Could not allocate memory for tdone\n");
            exit(1);
        }

        if (!(vdonne = (double **) realloc(vdonne,
                                           (gridpos + 1) * sizeof(double *)))) {
            printf("COuld not allocate memory for vdonne\n");
            exit(1);
        }

        if (!(derivv1 = (double **) realloc(derivv1,
                                            (gridpos + 1) * sizeof(double *)))) {
            printf("COuld not allocate memory for derivv1\n");
            exit(1);
        }

        if (!(derivv2 = (double **) realloc(derivv2,
                                            (gridpos + 1) * sizeof(double *)))) {
            printf("COuld not allocate memory for derivv2\n");
            exit(1);
        }

        if (!(derivv3 = (double **) realloc(derivv3,
                                            (gridpos + 1) * sizeof(double *)))) {
            printf("COuld not allocate memory for derivv3 \n");
            exit(1);
        }

        if (!(derivv4 = (double **) realloc(derivv4,
                                            (gridpos + 1) * sizeof(double *)))) {
            printf("COuld not allocate memory for derivv4\n");
            exit(1);
        }

        derivv1[gridpos] = (double *) calloc(n, sizeof(double));
        derivv2[gridpos] = (double *) calloc(n, sizeof(double));
        derivv3[gridpos] = (double *) calloc(n, sizeof(double));
        derivv4[gridpos] = (double *) calloc(n, sizeof(double));
        vdonne[gridpos] = (double *) calloc(n, sizeof(double));

        tdone[gridpos] = t;

        if (t - maxdel > tdone[0]) {
            while (t - maxdel >= tdone[gridstart])
                gridstart++;

            gridstart--;
        }

        memcpy(vdonne[gridpos], vnow, sizeof(*vnow) * n);

        /* put present derivv4 into future derivv1 */

        memcpy(derivv1[gridpos], derivv4[gridpos - 1], sizeof(**derivv4) * n);

    }

    /*       for (j=0; j<=gridpos;j++)
     {
     for (i=0; i<n; i++)
     printf("At time %f, The deriv1[%d][%d]=%.10f deriv2[%d][%d]"
     "=%.10f deriv3[%d][%d]=%.10f deriv4[%d][%d]=%.10f\n",tdone[j],j,i,
     derivv1[j][i],j,i, derivv2[j][i],j,i,derivv3[j][i],j,i,
     derivv4[j][i]);
     }
     */

    /*  printf("vs on grid:\n");*/

    /*      for (j=0; j<=gridpos;j++)
     {
     printf("%f %f\n",tdone[j],vdonne[j][0]);

     }
     */

    free(v[0]);
    free(v[1]);

    for (vc = 0; vc < 4; vc++) {
        for (dc = 0; dc < numdel; dc++)
            free(v_delayed[vc][dc]);
        free(v_delayed[vc]);
    }

    /* Put zeroes in derivv1[gridpos], all the rest are zero anyways */

    free(derivv1[gridpos]);
    derivv1[gridpos] = (double *) calloc(n, sizeof(double));

    free(v);
    free(vtemp);
    free(verror);

}

int SoDe::y_delayed(double ***vd, int n, double *rktimes, double *tau,
                    double *grid, double **vdone, double **deriv1,
                    double **deriv2, double **deriv3, double **deriv4,
                    int gridsize, double accu)
{

    int i, j, vc, dc, it;
    double t, ech;
    double *vtemp, *vnext, *vprev, *dummy;
    double *drv1, *drv2, *drv3, *drv4;
    double verror_max;

    static double
    /*  a1  =      0.0, */
    a2 = 0.5, a3 = 0.75;
    /*  a4  =     1.0, */

    static double b21 = 0.5,
    /*  b31 =     0.0, */
    b32 = 0.75;

    static double c1 = 2.0 / 9.0, c2 = 1.0 / 3.0, c3 = 4.0 / 9.0;

    vtemp = (double *) calloc(n, sizeof(double));
    vnext = (double *) calloc(n, sizeof(double));
    vprev = (double *) calloc(n, sizeof(double));
    drv1 = (double *) calloc(n, sizeof(double));
    drv2 = (double *) calloc(n, sizeof(double));
    drv3 = (double *) calloc(n, sizeof(double));
    drv4 = (double *) calloc(n, sizeof(double));

    ech = rktimes[3] - rktimes[0];

    /* First lets initialiaze all the the vdelays based on what is
     available to us without iteration. We will carry out iteration
     subsequently, only if t-mindel lies within the integration interval */

    for (vc = 0; vc < 4; vc++) {

        t = rktimes[vc];

        for (dc = 0; dc < numdel; dc++)
            if (tau[dc] == 0.)
                vd[vc][dc] = (double *)memcpy(vd[vc][dc], vdone[gridsize - 1],
                                    sizeof(double) * n);
            else if (t - tau[dc] <= grid[0])
                History(t - tau[dc], t, vd[vc][dc], n);
            else if (t - tau[dc] <= grid[gridsize - 1]) {

                j = 0;
                while ((j < gridsize) && (t - tau[dc] > grid[j]))
                    j++;

                if (j >= gridsize) {
                    printf("y_past:time requested not in grid,"
                           "bailing!\n");
                    exit(1);
                }

                if ((j == gridsize - 1) && (t - tau[dc] == grid[gridsize - 1])) {
                    vd[vc][dc] = (double *)memcpy(vd[vc][dc], vdone[gridsize - 1],
                                        sizeof(double) * n);
                    /*              for (i=0; i<n; i++)
                     printf("t=%f, t-del=%f, vdone[%d]=%f,"
                     "vpast[%d]=%f\n",
                     t,t-tau[dc],j,vdone[j][i],i,vd[vc][dc][i]);*/
                } else {

                    CE(t - tau[dc], vd[vc][dc], grid[j - 1], vdone[j - 1],
                       grid[j] - grid[j - 1], deriv1[j - 1], deriv2[j - 1],
                       deriv3[j - 1], deriv4[j - 1], n);
                    /*                  for (i=0; i<n; i++)
                     printf("t=%f, t-del=%f, vdone[%d]=%f,"
                     "vpast[%d]=%f\n",
                     t,t-tau[dc],j,vdone[j-1][i],i,vd[vc][dc][i]);*/
                }
            } else if (t - tau[dc] > grid[gridsize - 1]) {

                CE(t - tau[dc], vd[vc][dc], grid[gridsize - 2],
                   vdone[gridsize - 2], grid[gridsize - 1] - grid[gridsize - 2],
                   deriv1[gridsize - 2], deriv2[gridsize - 2],
                   deriv3[gridsize - 2], deriv4[gridsize - 2], n);
                /*              for (i=0; i<n; i++)
                 printf("IterInit [%.6f,%.6f], t-del=%.6f,"
                 "vd[%d][%d]=%.6f\n",
                 grid[gridsize-2],grid[gridsize-1],
                 t-tau[dc],vc,dc,vd[vc][dc][i]); */
            } else {

                printf("Requested point not to be found\n");
                exit(1);

            }
    }

    /* Now lets do the promised iteration, if required ofcourse! */

    if (rktimes[3] - mindel > grid[gridsize - 1]) {

        zygote.d_deriv(vdone[gridsize - 1], vd[0], rktimes[0], drv1, n);

        it = 0;
        verror_max = 100.;

        while ((it < 5) && (verror_max > 0.01 * accu)) {

            /*          printf("Iteration No.%d, error: %f\n",it,verror_max);*/

            for (i = 0; i < n; i++)
                vtemp[i] = vdone[gridsize - 1][i] + ech * (b21 * drv1[i]);

            zygote.d_deriv(vtemp, vd[1], rktimes[1], drv2, n);

            for (i = 0; i < n; i++)
                vtemp[i] = vdone[gridsize - 1][i] + ech * (b32 * drv2[i]);

            zygote.d_deriv(vtemp, vd[2], rktimes[2], drv3, n);

            for (i = 0; i < n; i++)
                vnext[i] = vdone[gridsize - 1][i]
                        + ech * (c1 * drv1[i] + c2 * drv2[i] + c3 * drv3[i]);

            zygote.d_deriv(vnext, vd[3], rktimes[3], drv4, n);

            for (vc = 0; vc < 4; vc++) {

                t = rktimes[vc];

                for (dc = 0; dc < numdel; dc++)
                    if (t - tau[dc] > grid[gridsize - 1]) {
                        CE(t - tau[dc], vd[vc][dc], grid[gridsize - 1],
                           vdone[gridsize - 1], ech, drv1, drv2, drv3, drv4, n);

                        /*                  for (i=0; i<n; i++)
                         printf("[%.6f,%.6f], t-del=%.6f,"
                         "vd1[%d]=%.6f,vd2[%d]=%.6f, vd3[%d]=%.6f,"
                         "vd4[%d]=%.6f, vnext[%d]=%.6f, t=%.6f\n",
                         rktimes[0],rktimes[3],
                         t-tau[dc],dc,vd[0][dc][0],
                         dc,vd[1][dc][0], dc,vd[2][dc][0],
                         dc,vd[3][dc][0],i,vnext[i],
                         grid[gridsize-1]);*/
                    }
            }

            if (it) {
                verror_max = 0.;
                for (i = 0; i < n; i++)
                    if (fabs(vprev[i]) > 10000. * DBL_EPSILON)
                        verror_max = max(fabs((vnext[i] - vprev[i]) / vprev[i]),
                                         verror_max);
                    else
                        verror_max = max(
                                fabs((vnext[i] - vprev[i]) / (10000.
                                        * DBL_EPSILON)),
                                verror_max);
            }

            dummy = vnext;
            vnext = vprev;
            vprev = dummy;

            it++;
        }

        if (it == 5)
            return 1;
    }

    free(vtemp);
    free(vnext);
    free(vprev);
    free(drv1);
    free(drv2);
    free(drv3);
    free(drv4);

    return 0;
}

int compare(const double *x, const double *y) {

    if (*x > *y)
        return 1;
    else if (*x < *y)
        return -1;
    else
        return 0;

}

double *SoDe::Construct_Discont_Array(double range, double *taus, int n,
                                      double *starts, int sn, int *disc_size)
{

    int M = 0, m;
    int i, j, k;
    double *delay_array = NULL;

    for (i = 0; i < n; i++) {

        if (taus[i] != 0.)
            m = (int) floor(range / taus[i]) + 1;
        else
            m = 1;

        if (m > 5)
            m = 5;

        /*      printf("%f %f m=%d\n",range/taus[i],floor(range/taus[i]),m);*/

        for (k = 0; k < sn; k++) {
            if (starts[k] > -taus[i]) {

                delay_array = (double *) realloc(delay_array,
                                                 (M + m) * sizeof(double));
                for (j = M; j < M + m; j++)
                    delay_array[j] = starts[k] + (j - M + 1) * taus[i];
                M += m;
            }
        }
    }

    /*  for (i=0;i<M;i++) printf("%1.14f\n",delay_array[i]);*/

    qsort((void *) delay_array, M, sizeof(double),
          (int (*)(const void*, const void*))compare);

    qsort((void *) taus, n, sizeof(double),
          (int (*)(const void*, const void*))compare);

     /* Now lets remove duplicates */

    i    = 0;

    while (i < M) {
        if (i < M - 1) {
            if ((fabs(delay_array[i] - delay_array[i + 1]) < 1E-10) || (delay_array[i]
                                                                                    > range)) {
                memmove((delay_array + i), (delay_array + i + 1),
                        (M - i - 1) * sizeof(double));
                /*              printf("Shifted %d elements to %d\n",M-i-1,i);*/
                M--;
                i--;
            }
        } else if ((i == M - 1) && (delay_array[i] > range)) {
            M--;
            i--;
        }

        i++;
        delay_array = (double *) realloc(delay_array, M * sizeof(double));
    }

    /*  for (i=0;i<M;i++) printf("Filtered:%1.14f\n",delay_array[i]);*/

    *disc_size = M;

    return delay_array;

}
void SoDe::CE(double t, double *vans, double tbegin, double *v_at_tbegin,
double ech, double *d1, double *d2, double *d3, double *d4, int n)
{
        static double

        e1 =   -4.0/3.0,
        e2 =   5.0/9.0,
        e3 =   -2.0/3.0,
        e4 =   -8.0/9.0;

        int i;
        double sigma;       /* for the continuous ext. v(n+sigma)*/
        double sigma_sq;               /* sigma's square */
        double dt,f1,f2,f3,f4;

        sigma = (t - tbegin)/ech;
        sigma_sq = sigma*sigma;
        dt = t - tbegin;

        f1 = dt*(1 + e1*sigma + e2*sigma_sq);
        f2 = dt*(sigma + e3*sigma_sq);
        f3 = dt*(e1*sigma - e4*sigma_sq);
        f4 = dt*(sigma_sq - sigma);

        for (i=0; i<n; i++)
            vans[i] = v_at_tbegin[i] + f1*d1[i];

        for (i=0; i<n; i++)
            vans[i] = vans[i] + f2*d2[i];

        for (i=0; i<n; i++)
            vans[i] = vans[i] - f3*d3[i];

        for (i=0; i<n; i++)
            vans[i] = vans[i] + f4*d4[i];

        return;

}

void SoDe::SetHistoryInterp(InterpObject interp_info)
{

    hist_interp_object.fact_discons     = interp_info.fact_discons;
    hist_interp_object.fact_discons_size= interp_info.fact_discons_size;
    hist_interp_object.func             = interp_info.func;
    hist_interp_object.slope            = interp_info.slope;
    hist_interp_object.maxsize          = interp_info.maxsize;
    hist_interp_object.maxtime          = interp_info.maxtime;

    int i;

    fact_discons_size = hist_interp_object.fact_discons_size +
                                extinp_interp_object.fact_discons_size;

    fact_discons = (double *) calloc(fact_discons_size,
                                                sizeof(double));

    for (i = 0; i < hist_interp_object.fact_discons_size; i++)
        fact_discons[i] = hist_interp_object.fact_discons[i];

    for (i = 0; i < extinp_interp_object.fact_discons_size; i++)
        fact_discons[i+hist_interp_object.fact_discons_size]
                            = extinp_interp_object.fact_discons[i];

    qsort((void *) fact_discons, fact_discons_size,
            sizeof(double),(int (*) (const void*,const void*)) compare);

}

void SoDe::SetExternalInputInterp(InterpObject interp_info)
{

    extinp_interp_object.fact_discons   = interp_info.fact_discons;
    extinp_interp_object.fact_discons_size= interp_info.fact_discons_size;
    extinp_interp_object.func               = interp_info.func;
    extinp_interp_object.slope              = interp_info.slope;
    extinp_interp_object.maxsize            = interp_info.maxsize;
    extinp_interp_object.maxtime            = interp_info.maxtime;

}

void SoDe::History(double t, double t_size, double *yd, int n)
{

    int j,k;
    double *blug;
    double t_interp, t_diff;

/*  printf("Going from %d to %d, time:%f time for size:%f\n",maxsize,n,t,t_size);*/

    blug = (double *) calloc(hist_interp_object.maxsize,
                                                sizeof(double));

    k = -1;
    do {
            k++;

    } while ( (k < hist_interp_object.slope.size ) &&
                        (t > hist_interp_object.slope.array[k].time));
    if (k == 0)
        t_interp = hist_interp_object.slope.array[k].time;
    else {

        k--;
        t_interp = t;

    }

    t_diff = t_interp - hist_interp_object.slope.array[k].time;

    for (j=0; j < hist_interp_object.maxsize; j++)
    {

        blug[j] = hist_interp_object.func.array[k].state.array[j]
                    + hist_interp_object.slope.array[k].state.array[j]
                    *t_diff;

    }

    if (n >= hist_interp_object.maxsize)
        Go_Forward(yd, blug, TheMaternal.GetStartLinIndex(t_size),
                   TheMaternal.GetStartLinIndex(hist_interp_object.maxtime),
                   defs.ngenes);
    else
        Go_Backward(yd, blug, TheMaternal.GetStartLinIndex(t_size),
                    TheMaternal.GetStartLinIndex(hist_interp_object.maxtime),
                    defs.ngenes);
    free(blug);

    return;
}

void SoDe::ExternalInputs(double t, double t_size, double *yd, int n)
{

    int j,k;
    double *blug;
    double t_interp, t_diff;

/*  printf("Going from %d to %d, time:%f time for size:%f\n",maxsize,n,t,t_size);*/

    blug = (double *) calloc(extinp_interp_object.maxsize,
                                                sizeof(double));

    k = -1;
    do {
            k++;

    } while ( (k < extinp_interp_object.slope.size ) &&
                        (t > extinp_interp_object.slope.array[k].time));
    if (k == 0)
        t_interp = extinp_interp_object.slope.array[k].time;
    else {

        k--;
        t_interp = t;

    }

    t_diff = t_interp - extinp_interp_object.slope.array[k].time;

    for (j=0; j < extinp_interp_object.maxsize; j++)
    {

        blug[j] = extinp_interp_object.func.array[k].state.array[j]
                + extinp_interp_object.slope.array[k].state.array[j]
                *t_diff;

    }

    if (n >= extinp_interp_object.maxsize)
        Go_Forward(yd, blug, TheMaternal.GetStartLinIndex(t_size),
                   TheMaternal.GetStartLinIndex(extinp_interp_object.maxtime),
                   defs.egenes);
    else
        Go_Backward(yd, blug, TheMaternal.GetStartLinIndex(t_size),
                    TheMaternal.GetStartLinIndex(extinp_interp_object.maxtime),
                    defs.egenes);
    free(blug);

    return;
}

void SoDe::DivideHistory(double t1, double t2)
{

    double *blug;
    int i,size;

    if ((size = zygote.get_NNucs(t2) * zygote.get_ngenes())
            > zygote.get_NNucs(t1) * zygote.get_ngenes())
        for (i=0; i <= gridpos; i++)
        {

            blug = (double *) calloc(size, sizeof(double));

            Go_Forward(blug, vdonne[i], TheMaternal.GetStartLinIndex(t2),
                       TheMaternal.GetStartLinIndex(t1),
                       defs.ngenes);

            free(vdonne[i]);
            vdonne[i] = blug;

            blug = (double *) calloc(size, sizeof(double));

            Go_Forward(blug, derivv1[i], TheMaternal.GetStartLinIndex(t2),
                       TheMaternal.GetStartLinIndex(t1),
                       defs.ngenes);

            free(derivv1[i]);
            derivv1[i] = blug;

            blug = (double *) calloc(size, sizeof(double));

            Go_Forward(blug, derivv2[i], TheMaternal.GetStartLinIndex(t2),
                       TheMaternal.GetStartLinIndex(t1),
                       defs.ngenes);

            free(derivv2[i]);
            derivv2[i] = blug;

            blug = (double *) calloc(size, sizeof(double));

            Go_Forward(blug, derivv3[i], TheMaternal.GetStartLinIndex(t2),
                       TheMaternal.GetStartLinIndex(t1),
                       defs.ngenes);

            free(derivv3[i]);
            derivv3[i] = blug;

            blug = (double *) calloc(size, sizeof(double));

            Go_Forward(blug, derivv4[i], TheMaternal.GetStartLinIndex(t2),
                       TheMaternal.GetStartLinIndex(t1),
                       defs.ngenes);

            free(derivv4[i]);
            derivv4[i] = blug;

        }

}

void SoDe::Go_Forward(double *output, double *input, int output_ind, int
input_ind, int num_genes)
{

    double *y;
    int output_lin, input_lin, size, newsize;
    int k, ap, i, j, ii;

    output_lin = TheMaternal.Index2StartLin(output_ind);
    newsize = TheMaternal.Index2NNuc(output_ind)*num_genes;

/*  printf("output lineage start, indices:%d, %d, %d\n",output_lin,output_ind,input_ind);*/

    if (output_ind < input_ind - 1)
    {
        size = TheMaternal.Index2NNuc(output_ind+1)*num_genes;
        y = (double *) calloc(size, sizeof(double));
/*      printf("Passing on to another fwd with targets %d %d %d\n",size,output_ind+1,input_ind);*/
        Go_Forward(y,input,output_ind+1,input_ind,num_genes);
    } else if  (output_ind == input_ind - 1)
    {
        size = TheMaternal.Index2NNuc(input_ind)*num_genes;
        y = (double *) calloc(size, sizeof(double));
/*      printf("Goin' to do the tranfer:%d %d\n",size,newsize);*/
        y = (double *)memcpy(y,input,size*sizeof(double));
    } else if (output_ind == input_ind) {
        output = (double *)memcpy(output,input,newsize*sizeof(double));
        return;
    }
    else error("You are trying to go from nnucs %d to %d!",
               TheMaternal.Index2NNuc(input_ind),
               TheMaternal.Index2NNuc(output_ind));


      for (j=0; j < size; j++) {

    k  = j % num_genes;     /* k: index of gene k in current nucleus */
    ap = j / num_genes;      /* ap: rel. nucleus position on AP axis */

/* evaluate ii: index of anterior daughter nucleus */

    if ( output_lin % 2 )
      ii = 2 * ap * num_genes + k - num_genes;
    else
      ii = 2 * ap * num_genes + k;

/* skip the first most anterior daughter nucleus in case output_lin is odd */

        if ( ii >= 0 )
      output[ii] = y[j];

/* the second daughter only exists if it is still within the region */

    if ( ii + num_genes < newsize )
      output[ii+num_genes] =  y[j];
      }

    free(y);

    return;
}


void SoDe::Go_Backward(double *output, double *input, int output_ind, int
input_ind, int num_genes)
{

    double *y;
    int output_lin, input_lin, size, newsize;
    int k, ap, i, j, ii;

    output_lin = TheMaternal.Index2StartLin(output_ind);
    input_lin  = TheMaternal.Index2StartLin(input_ind);
    newsize = TheMaternal.Index2NNuc(output_ind)*num_genes;

/*  printf("output lineage start, indices:%d, %d, %d\n",output_lin,output_ind,input_ind);*/

    if (output_ind > input_ind + 1)
    {
        size = TheMaternal.Index2NNuc(output_ind-1)*num_genes;
        y = (double *) calloc(size, sizeof(double));
/*      printf("Passing on to another bkd with targets %d %d %d\n",size,output_ind-1,input_ind);*/
        Go_Backward(y,input,output_ind-1,input_ind,num_genes);
        input_lin  = TheMaternal.Index2StartLin(output_ind-1);
    } else if  (output_ind == input_ind + 1 )
    {
        size = TheMaternal.Index2NNuc(input_ind)*num_genes;
        y = (double *) calloc(size, sizeof(double));
/*      printf("Goin' to do the tranfer\n");*/
        memcpy(y,input,size*sizeof(double));
    } else if (output_ind == input_ind){
        memcpy(output,input,newsize*sizeof(double));
        return;
    }
    else error("You are trying to go from nnucs %d to %d!",
               TheMaternal.Index2NNuc(input_ind),
               TheMaternal.Index2NNuc(output_ind));

      for (j=0; j < newsize; j++) {

    k  = j % num_genes;     /* k: index of gene k in current nucleus */
    ap = j / num_genes;      /* ap: rel. nucleus position on AP axis */

/* evaluate ii: index of anterior daughter nucleus */

    if ( input_lin % 2 )
      ii = 2 * ap * num_genes + k - num_genes;
    else
      ii = 2 * ap * num_genes + k;

/* skip the first most anterior daughter nucleus in case output_lin is odd */

    if (ii < 0)
        output[j] = y[ii + num_genes];

    if (( ii >= 0 ) && (ii + num_genes < size))
        output[j] = .5*(y[ii] + y[ii+num_genes]);

    if (ii + num_genes >= size)
        output[j] = y[ii];

      }

    free(y);

    return;
}

solver *SolverFactory(zygotic &zygote, int debug, const char *name)
{
    if (!(strcmp(name, "a")) || !(strcmp(name, "Adams")))
        return new Adams(zygote, debug);
    else if (!(strcmp(name, "bd")) || !(strcmp(name, "BaDe")))
        return new BaDe(zygote, debug);
    else if (!(strcmp(name, "bs")) || !(strcmp(name, "BuSt")))
        return new BuSt(zygote, debug);
    else if (!(strcmp(name, "e")) || !(strcmp(name, "Euler")))
        return new Euler(zygote, debug);
    else if (!(strcmp(name, "h")) || !(strcmp(name, "Heun")))
        return new Heun(zygote, debug);
    else if (!(strcmp(name, "mi")) || !(strcmp(name, "m"))
             || !(strcmp(name, "Milne")))
        return new Milne(zygote, debug);
    else if (!(strcmp(name, "me")) || !(strcmp(name, "Meuler")))
        return new Meuler(zygote, debug);
    else if (!(strcmp(name, "r4")) || !(strcmp(name, "r"))
             || !(strcmp(name, "Rk4")))
        return new Rk4(zygote, debug);
    else if (!(strcmp(name, "r2")) || ! (strcmp(name, "Rk2")))
        return new Rk2(zygote, debug);
    else if (!(strcmp(name, "rck")) || !(strcmp(name, "Rkck")))
        return new Rkck(zygote, debug);
    else if (!(strcmp(name, "rf")) || !(strcmp(name, "Rkf")))
        return new Rkf(zygote, debug);
    else if (!(strcmp(name, "sd")) || !(strcmp(name, "SoDe")))
        return new SoDe(zygote, debug);
    else
        return NULL;

}
