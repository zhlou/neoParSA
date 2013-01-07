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
    struct StatData
    {
        double s;
        double mean;
        double vari;
        double energy;
        long success;
        long moves;
    };
    MPI_Comm comm;
    int nnodes;
    int rank;
    StatData l_stat, *local_stat_buf;
    MPI_Win stat_win;
};

#endif /* PSTATS_H_ */
