/*
 * PStats.h
 *
 *  Created on: Jan 1, 2013
 *      Author: zhlou
 */

#ifndef PSTATS_H_
#define PSTATS_H_
#include <mpi.h>

class PStats
{
public:
    PStats(MPI_Comm thecomm);
    virtual ~PStats();
};

#endif /* PSTATS_H_ */
