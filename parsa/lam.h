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
    double nvari; //tau * sigma ^2, the sum of energy difference squared within proc_tau
    int accept;
    bool frozen();
    void updateStats(bool accept, double delta);
    void cool_s();
    virtual void updateEstimates();


};

#endif /* LAM_H_ */
