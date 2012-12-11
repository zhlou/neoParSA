/*
 * lam.h
 *
 *  Created on: Dec 8, 2012
 *      Author: zhlou
 */

#ifndef LAM_H_
#define LAM_H_

#include "annealer.h"
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
    double w_a, w_b;
    double usx, usy, usxx, usxy, usyy, usum;
    double vsx, vsy, vsxx, vsxy, vsyy, vsum;
    double estimate_mean;
    double estimate_sd;


};

#endif /* LAM_H_ */
