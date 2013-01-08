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
    l_stat.moves = -1;
    MPI_Win_create(&l_stat, sizeof(StatData), sizeof(StatData), MPI_INFO_NULL,
            comm, &stat_win);
    MPI_Comm_group(comm, &group);
}

void PStats::CommSegment(double mean, double vari, long success, long steps,
        invLinearFit* fit_mean, invLinearFit* fit_sd)
{
    //pack statistics
    //l_stat.s = S;
    l_stat.mean = mean;
    l_stat.vari = vari;
    l_stat.success = success;
    l_stat.moves = steps;
    //l_stat.energy = energy;
    MPI_Win_post(group, 0, stat_win);
    MPI_Win_start(group, 0, stat_win);
    for (int j = 0; j < nnodes; j++)
      {
        if (j != rank)
        {
          MPI_Get(local_stat_buf + j, sizeof(StatData), MPI_BYTE, j, 0,
                  sizeof(StatData), MPI_BYTE, stat_win);
        }
      }
    MPI_Win_complete(stat_win);
    MPI_Win_wait(stat_win);
}

bool PStats::CommCheckFrozen(int freeze_cnt)
{
}

PStats::~PStats()
{
    MPI_Win_free(&stat_win);
    delete[] local_stat_buf;
}

