/*
 * PStats.cpp
 *
 *  Created on: Jan 1, 2013
 *      Author: zhlou
 */

#include "PStats.h"
#include <mpi.h>

PStats::PStats(MPI_Comm thecomm, int in_nnodes, int in_rank) :
    comm(thecomm), nnodes(in_nnodes), rank(in_rank)
{
    local_stat_buf = new StatData[nnodes];
    l_stat.s=-1;
    MPI_Win_create(&l_stat, sizeof(StatData), sizeof(StatData), MPI_INFO_NULL,
            comm, &stat_win);
}

void PStats::CommSegment(double mean, double vari, invLinearFit* fit_mean,
        invLinearFit* fit_sd)
{
}

bool PStats::CommCheckFrozen(int freeze_cnt)
{
}

PStats::~PStats()
{
    // TODO Auto-generated destructor stub
}

