/*
 * periodBest.hpp
 *
 *  Created on: Aug 20, 2014
 *      Author: zhlou
 */

#ifndef PERIODBEST_HPP_
#define PERIODBEST_HPP_

#include <stdexcept>
#include <string>
#include <mpi.h>
#include "xmlUtils.h"

template <class Problem>
const char * periodBest<Problem>::name = "periodBest";

template <class Problem>
periodBest<Problem>::Param::Param(xmlNode *docroot)
{
    xmlNode *section = getSectionByName(docroot, "periodBest");
    if (section == NULL)
        throw std::runtime_error(std::string("Error: cannot find section periodBest"));

    frequency = getPropInt(section, "frequency");
}

template <class Problem>
periodBest<Problem>::periodBest(Problem &problem, MPIState const&mpiState,
                                unirandom &, const Param &param) :
                                problem(problem), mpi(mpiState),
                                frequency(param.frequency), counter(0),
                                buf_size(problem.getStateSize())
{

    MPI_Alloc_mem(buf_size, MPI_INFO_NULL, &state_buf);
}


template <class Problem>
periodBest<Problem>::~periodBest()
{
    MPI_Free_mem(state_buf);
}

template <class Problem>
mixState periodBest<Problem>::Mix(aState &state)
{

    if (counter % frequency) {
        ++ counter;
        return mixState();
    }
    ++ counter;

    doubleint energyRank;
    energyRank.energy=state.energy;
    energyRank.rank=mpi.rank;
    MPI_Allreduce(MPI_IN_PLACE, &energyRank, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                  mpi.comm);
    if (mpi.rank == energyRank.rank)
        problem.serialize(state_buf);
    MPI_Bcast(state_buf, buf_size, MPI_BYTE, energyRank.rank, mpi.comm);
    debugOut << state.step_cnt << " " << state.s << " "
             << energyRank.energy << "  " << energyRank.rank;
    debugOut << std::endl;
    if (mpi.rank != energyRank.rank) {
        problem.deserialize(state_buf);
        state.energy = energyRank.energy;
        return mixState(energyRank.rank);
    }
    return mixState();


}

#endif /* PERIODBEST_HPP_ */
