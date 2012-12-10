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
    int accept;
    bool frozen();
    void updateStep(bool accept, double delta);
    void cool_s();
    void resetStats();
    bool inSegment();
    void updateSegment();


};

#endif /* LAM_H_ */
