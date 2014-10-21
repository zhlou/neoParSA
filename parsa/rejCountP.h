/*
 * rejCountP.h
 *
 *  Created on: Apr 5, 2013
 *      Author: zhlou
 */

#ifndef REJCOUNTP_H_
#define REJCOUNTP_H_

#include <libxml/tree.h>
#include "rejCount.h"
#include "MPIState.h"

class rejCountP
{
public:
    class Param {
    public:
        rejCount::Param serParam;
        const MPIState &mpi;
        Param(xmlNode *root, const MPIState &mpiState, debugStatus st=ignore,
              const char * debugname=NULL);
    };
    rejCountP(const Param &param, const MPIState &mpiState);
    virtual ~rejCountP();
    void updateStep(bool accept, const aState &state)
    {serCount.updateStep(accept, state);}
    bool frozen(const aState &state) const;
private:
    const MPIState &mpi;
    rejCount serCount;
};

#endif /* REJCOUNTP_H_ */
