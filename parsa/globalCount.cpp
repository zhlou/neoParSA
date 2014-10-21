/*
 * globalCount.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: zhlou
 */

#include "globalCount.h"

bool globalCount::frozen(const aState& state)
{
    int total_count;
    serCount.frozen(state);
    MPI_Allreduce(&serCount.freeze_cnt, &total_count, 1, MPI_INT, MPI_SUM,
            mpi.comm);
    return (total_count>=globalcritcnt);
}
