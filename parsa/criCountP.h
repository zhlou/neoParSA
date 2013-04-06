/*
 * criCountP.h
 *
 *  Created on: Apr 5, 2013
 *      Author: zhlou
 */

#ifndef CRICOUNTP_H_
#define CRICOUNTP_H_

#include "criCount.h"
#include "MPIState.h"

class criCountP {
public:
    class Param
    {
        criCount::Param serParam;
        const MPIState &mpi;
        Param(xmlNode *root, const MPIState &MPIState);
    };
    criCountP(const Param &param, const MPIState &mpiState);
    ~criCountP();
    void updateStep(bool accept, const aState &state){}
    bool frozen(const aState &state);
private:
    criCount serCount;
    const MPIState &mpi;
};

#endif /* CRICOUNTP_H_ */
