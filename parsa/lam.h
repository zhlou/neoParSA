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
    int init_loop;
    bool frozen();
    double freeze_crit;
    double old_energy;
    int freeze_cnt;
    int cnt_crit;
    void updateStep(bool accept);
    void updateS();
    bool inSegment();
    void updateSegment();
    void resetSegmentStats();
    void resetLam();
    void updateLam();
    // below are lam parameters
    invLinearFit *fit_mean, *fit_sd;
    double alpha;


};

#endif /* LAM_H_ */
