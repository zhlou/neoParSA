/*
 * solvers.h
 *
 *  Created on: Feb 6, 2013
 *      Author: zhlou
 */

#ifndef SOLVERS_H_
#define SOLVERS_H_

class zygotic;

class solver
{
public:
    solver(const zygotic &in_zy, int in_debug) : zygote(in_zy), debug(in_debug) {};
    virtual ~solver() = 0;
    virtual void ps(double *, double *, double, double, double, double, int,
            FILE *) = 0;
protected:
    void WriteSolvLog(char *solver, double tin, double tout, double h, int n,
              int nderivs, FILE *slog);
    int debug;
    zygotic &zygote;
};

class Euler : public solver
{
public:
    // for newer compilers support c++11, we can simply use
    // using solver::solver;
    Euler(const zygotic &in_zy, int in_debug) : solver(in_zy, in_debug) {};
    void ps(double *vin, double *vout, double tin, double tout,
            double stephint, double accuracy, int n, FILE *slog);
};

#endif /* SOLVERS_H_ */
