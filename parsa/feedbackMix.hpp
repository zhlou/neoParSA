/*
 * feedbackMix.hpp
 *
 *  Created on: Aug 27, 2014
 *      Author: zhlou
 */

#ifndef FEEDBACKMIX_HPP_
#define FEEDBACKMIX_HPP_

#include <cmath>
#include <exception>
#include <mpi.h>
#include "xmlUtils.h"

template<class Problem>
const char * feedbackMix<Problem>::name = "feedbackMix";

template<class Problem>
feedbackMix<Problem>::Param::Param(xmlNode *root) : target(0.5)
{
    xmlNode *section = getSectionByName(root, "feedbackMix");
    interval = getPropInt(section, "interval");
    try {
        target = getPropDouble(section, "target");
    } catch(std::exception &e) {
        // ignore
    }
}

template<class Problem>
feedbackMix<Problem>::feedbackMix(Problem &problem, const MPIState &mpiState,
                                  unirandom &rand, const Param &param):
        mix(problem, mpiState, rand), tau_count(0),
        mpi(mpiState), target(param.target), interval(param.interval)
{
    adoptArray = new int[mpi.nnodes];
}

template<class Problem>
feedbackMix<Problem>::~feedbackMix()
{
    delete []adoptArray;
}

template<class Problem>
mixState feedbackMix<Problem>::Mix(aState &state)
{
    ++tau_count;
    if ((tau_count % interval) != 0)
        return mixState();
    int i,p, nadopt=0;
    const double C = 1.38629436112; // 2ln(2)
    tau_count = 0;
    mix.calProbTab(state);
    p = mix.getPartner();
    state.energy = mix.adoptState(p);
    for (i = 0; i < mpi.nnodes; ++ i)
        adoptArray[i] = 0;
    adoptArray[p] = 1;
    MPI_Allreduce(MPI_IN_PLACE, adoptArray, mpi.nnodes, MPI_INT, MPI_LOR,
                  mpi.comm);
    for (i = 0; i < mpi.nnodes; ++i)
        nadopt+=adoptArray[i];
    double adoptRate = (double)(nadopt) / mpi.nnodes;
    interval = lround(interval * std::exp(C*(adoptRate - target)));
    if (interval < 1)
        interval = 1;
    debugOut << state.step_cnt << " " << state.s << " "
             << adoptRate << " " << interval << std::endl;
    return mixState(p);

}

#endif /* FEEDBACKMIX_HPP_ */
