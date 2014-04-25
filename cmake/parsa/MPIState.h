/*
 * MPIState.h
 *
 *  Created on: Jan 31, 2013
 *      Author: zhlou
 */

#ifndef MPISTATE_H_
#define MPISTATE_H_

#include <mpi.h>

struct MPIState
{
    MPI_Comm comm;
    MPI_Group group;
    int nnodes;
    int rank;
};


#endif /* MPISTATE_H_ */
