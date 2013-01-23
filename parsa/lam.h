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

class lam
{
public:
    lam(xmlNode *root);
    virtual ~lam();
    double getInitS();
    int getInitLoop();
    void initStats();
    bool frozen();
    void updateInitStep(bool accept, double energy);
    void resetSegmentStats();
    void updateStep(bool accept, double energy);
    double updateS(double s);
    bool inSegment();
    virtual void updateSegment();

protected:
    int proc_tau;
    double acc_ratio;
    double vari;
    double mean;
    int success;


    double freeze_crit;
    double energy;
    double s;
    double old_energy;
    int freeze_cnt;
    int cnt_crit;

    static const double UNINITIALIZED;


    void collectStats();
    virtual void collectInitStats();
    void resetLam();
    void updateLam();
    // below are lam parameters
    invLinearFit *fit_mean, *fit_sd;
    double alpha;
    double w_mean;
    double w_sd;

private:
    double lambda;
    int init_loop;
    double init_S;
    long step_cnt;

};

#endif /* LAM_H_ */
