/*
 * intervalBest.hpp
 *
 *  Created on: Oct 8, 2015
 *      Author: zhlou
 */

#ifndef PARSA_INTERVALBEST_HPP_
#define PARSA_INTERVALBEST_HPP_

#include <mpi.h>
#include "xmlUtils.h"

template <class Problem>
const char *intervalBest<Problem>::name = "intervalBest";

template<class Problem>
intervalBest<Problem>::Param::Param(xmlNode* root)
{
	xmlNode *section = getSectionByName(root, "intervalBest");
	interval = getPropInt(section, "interval");
}

template<class Problem>
intervalBest<Problem>::intervalBest(Problem& in_problem,
		const MPIState& mpiState, unirandom& in_rand, const Param& param) :
		problem(in_problem),interval(param.interval), mpi(mpiState),
		rnd(in_rand), count(0)
{
	buf_size = problem.getStateSize();
	MPI_Alloc_mem(buf_size, MPI_INFO_NULL, &state_buf);
}

template<class Problem>
intervalBest<Problem>::~intervalBest()
{
	MPI_Free_mem(state_buf);
}

template<class Problem>
mixState intervalBest<Problem>::Mix(aState& state)
{
	++ count;
	if (interval != count)
		return mixState();
	count = 0;
	struct {
		double energy;
		int rank;
	} doubleint;
	doubleint.energy = state.energy;
	doubleint.rank = mpi.rank;
	MPI_Allreduce(MPI_IN_PLACE, &doubleint, 1, MPI_DOUBLE_INT, MPI_MINLOC,
			mpi.comm);
	if (mpi.rank == doubleint.rank) {
		problem.serialize(state_buf);
	}
	MPI_Bcast(state_buf, buf_size, MPI_BYTE, doubleint.rank, mpi.comm);
	if (mpi.rank != doubleint.rank) {
		problem.deserialize(state_buf);
		state.energy = doubleint.energy;
		return mixState(doubleint.rank);
	}
	return mixState();

}



#endif /* SOURCE_DIRECTORY__PARSA_INTERVALBEST_HPP_ */
