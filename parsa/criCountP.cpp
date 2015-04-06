/*
 * criCountP.cpp
 *
 *  Created on: Apr 5, 2013
 *      Author: zhlou
 */

#include "criCountP.h"

criCountP::criCountP(const Param& param, const MPIState& mpiState) :
    serCount(param.serParam), mpi(mpiState), interval(param.serParam.interval),
        step_cnt(0)
{
}

criCountP::~criCountP() {
    // TODO Auto-generated destructor stub
}

bool criCountP::frozen(const aState &state)
{
    step_cnt ++;
    if (interval != step_cnt){
        return false;
    }
    step_cnt = 0;
    int local_freeze = serCount.checkFrozen(state);
    int global_flag;
    MPI_Allreduce(&local_freeze, &global_flag, 1, MPI_INT, MPI_SUM, mpi.comm);
    return (bool)global_flag;
}

criCountP::Param::Param(xmlNode* root) :
        serParam(root)
{
}
