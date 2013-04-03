/*
 * pSimpleSchedule.h
 *
 *  Created on: Apr 2, 2013
 *      Author: zhlou
 */

#ifndef PSIMPLESCHEDULE_H_
#define PSIMPLESCHEDULE_H_

#include "simpleScheduler.h"
#include "MPIState.h"

class pSimpleSchedule : public simpleSchedule
{
public:
    pSimpleSchedule(xmlNode *root, const MPIState &mpiState);
    static const char *name;
    bool frozen(const aState);
protected:
    const MPIState &mpi;
    unsigned cnt;
};

#endif /* PSIMPLESCHEDULE_H_ */
