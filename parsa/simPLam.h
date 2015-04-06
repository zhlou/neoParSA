/*
 * simPLam.h
 *
 *  Created on: Dec 17, 2013
 *      Author: zhlou
 *
 *  Simplified Parallel Lam Schedule
 *  var & std are fixed over proc_tau iterations
 */


#ifndef SIMPLAM_H_
#define SIMPLAM_H_

#include <mpi.h>
#include <libxml/tree.h>
#include "MPIState.h"
#include "aState.h"
#include "dynDebug.h"

class simPLam {
private:
    dynDebug debugOut;
    double acc_ratio;
    double alpha;
    double mean;
    double sd;
    double sum;
    double sumsq;
    int success;
    int count;
    const int proc_tau;
    const double lambda;
    const MPIState &mpiState;

    void calcStats(int nsteps);

public:
    class Param {
    public:
        Param(xmlNode *root, debugStatus in_st=ignore, const char *name=NULL);
        int proc_tau;
        double lambda;
        debugStatus st;
        const char * outname;
    };
    simPLam(Param param, const MPIState &);
    ~simPLam();
    void initStats(const aState &);
    void initStats(double initMean, double initVar, double initAccRatio, 
            const aState &state);
    void updateInitStep(bool, const aState &);
    void resetSegmentStats();
    void updateStep(bool, const aState &);
    double updateS(const aState &state);
    bool inSegment(aState state){return !((state.step_cnt % proc_tau) == 0);}
    void updateSegment(const aState &);
    void updateStats(const aState &state);
    void setDebug(debugStatus st, const char* outname=NULL)
        {debugOut.setDebug(st, outname);}
    static const char * name;



};



#endif /* SIMPLAM_H_ */
