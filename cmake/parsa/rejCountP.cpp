/*
 * rejCountP.cpp
 *
 *  Created on: Apr 5, 2013
 *      Author: zhlou
 */

#include "rejCountP.h"


rejCountP::rejCountP(const rejCountP::Param& param, const MPIState& mpiState) :
    serCount(param.serParam), mpi(mpiState)
{
}

rejCountP::~rejCountP() {
    // TODO Auto-generated destructor stub
}

bool rejCountP::frozen(const aState& state) const
{
    unsigned tot_rej;
    int local_rej = serCount.reject_cnt;
    MPI_Allreduce(&local_rej, &tot_rej, 1, MPI_INT, MPI_SUM, mpi.comm);

    return (tot_rej >= serCount.max_rej * mpi.nnodes);

}

rejCountP::Param::Param(xmlNode* root,const MPIState &mpiState, debugStatus st,
                        const char* name):
        serParam(root, st, name), mpi(mpiState)
{
}
