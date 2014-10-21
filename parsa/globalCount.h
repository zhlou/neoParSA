/*
 * globalCount.h
 *
 *  Created on: Dec 18, 2013
 *      Author: zhlou
 */

#ifndef GLOBALCOUNT_H_
#define GLOBALCOUNT_H_

#include "criCount.h"
#include "MPIState.h"

class globalCount {
private:
    criCount serCount;
    const MPIState &mpi;
    const int globalcritcnt;
public:
    class Param
    {
    public:
        criCount::Param serParam;
        Param(xmlNode *root) : serParam(root) {};
    };
    globalCount(const Param &param, const MPIState &mpiState) :
        serCount(param.serParam), mpi(mpiState),
        globalcritcnt(serCount.cnt_crit * mpi.nnodes)
    {
    }
    ~globalCount(){}
    void updateStep(bool accept, const aState &state){}
    bool frozen(const aState &state);
};

#endif /* GLOBALCOUNT_H_ */
