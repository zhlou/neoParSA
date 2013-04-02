/*
 * pSimpleSchedule.cpp
 *
 *  Created on: Apr 2, 2013
 *      Author: zhlou
 */

#include "pSimpleSchedule.h"

const char *pSimpleSchedule::name = "parallel exponential";

pSimpleSchedule::pSimpleSchedule(xmlNode *root, const MPIState &mpiState) :
    simpleSchedule(root), mpi(mpiState)
{
}


bool pSimpleSchedule::frozen(aState)
{
    unsigned tot_rej;
    MPI_Allreduce(&reject_cnt, &tot_rej, 1, MPI_INT, MPI_SUM, mpi.comm);
    return (tot_rej >= max_rej * mpi.nnodes);
}
