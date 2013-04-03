/*
 * pSimpleSchedule.cpp
 *
 *  Created on: Apr 2, 2013
 *      Author: zhlou
 */

#include "pSimpleSchedule.h"
#include "utils.h"
#include <exception>

const char *pSimpleSchedule::name = "parallel exponential";

pSimpleSchedule::pSimpleSchedule(xmlNode *root, const MPIState &mpiState) :
    simpleSchedule(root), mpi(mpiState), cnt(0)
{
    check_freq = 100;
    try {
        check_freq = getPropInt(xmlsection, "check_freq");
    } catch (const std::exception &e) {

    }

}


bool pSimpleSchedule::frozen(aState)
{
    ++ cnt;
    if (cnt % check_freq == 0) {
        cnt = 0;
        unsigned tot_rej;
        MPI_Allreduce(&reject_cnt, &tot_rej, 1, MPI_INT, MPI_SUM, mpi.comm);
        return (tot_rej >= max_rej * mpi.nnodes);
    }
    return false;

}
