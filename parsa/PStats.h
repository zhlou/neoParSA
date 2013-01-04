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
    PStats(MPI_Comm thecomm, int in_nnodes, int in_rank);
    void CommSegment(double mean, double vari, invLinearFit *fit_mean,
            invLinearFit *fit_sd);
    bool CommCheckFrozen(int freeze_cnt);

    virtual ~PStats();
private:
    MPI_Comm comm;
    int nnodes;
    int rank;
};

#endif /* PSTATS_H_ */
