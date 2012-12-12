/*
 * lam.h
 *
 *  Created on: Dec 8, 2012
 *      Author: zhlou
 */

#ifndef LAM_H_
#define LAM_H_

#include "annealer.h"
#include "invLinearFit.h"
#include <mpi.h>

class lam: public annealer
{
public:
    lam(movable *, xmlNode*);
    virtual ~lam();
protected:
    int proc_tau;
    double acc_ratio;
    double vari;
    double mean;
    int success;
    bool frozen();
    void updateStep(bool accept, double delta);
    void updateS();
    bool inSegment();
    void updateSegment();
    void resetSegmentStats();
    void resetLam();
    void updateLam();
    // below are lam parameters
    invLinearFit *fit_mean, *fit_sd;
    double estimate_mean;
    double estimate_sd;
    double alpha;


};

#endif /* LAM_H_ */
