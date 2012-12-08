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
    double loop();
protected:
    int proc_tau;
    bool frozen();

};

#endif /* LAM_H_ */
