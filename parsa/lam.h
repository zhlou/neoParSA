/*
 * lam.h
 *
 *  Created on: Dec 8, 2012
 *      Author: zhlou
 */

#ifndef LAM_H_
#define LAM_H_

//#include "annealer.h"
#include "aState.h"
#include "invLinearFit.h"
#include <libxml/tree.h>

class lam
{
public:
    lam(xmlNode *root);
    virtual ~lam();
    //double getInitS();
    //int getInitLoop();
    bool frozen(const aState state);
    void updateInitStep(bool accept, aState state);
    void initStats(aState state);
    void resetSegmentStats();
    void updateStep(bool accept, aState state);
    double updateS(aState state);
    bool inSegment(aState state);
    void updateSegment(aState state);

protected:
    int proc_tau;
    double acc_ratio;
    double vari;
    double mean;
    int success;


    double freeze_crit;
    //double energy;
    //double s;
    double old_energy;
    int freeze_cnt;
    int cnt_crit;

    static const double UNINITIALIZED;


    void collectStats();
    virtual void collectInitStats(unsigned long init_loop);
    void resetLam();
    void updateLam();
    virtual void updateEstimators(double s);
    void local_frozen(const aState& state);
    virtual bool global_frozen();

    // below are lam parameters
    invLinearFit *fit_mean, *fit_sd;
    double alpha;
    double w_mean;
    double w_sd;

private:
    double lambda;
    //int init_loop;
    //double init_S;
    //long step_cnt;

};

#endif /* LAM_H_ */
