/* 
 * File:   rateMix.hpp
 * Author: zhlou
 *
 * Created on April 16, 2015, 1:18 PM
 */

#ifndef RATEMIX_HPP
#define	RATEMIX_HPP


#include <cmath>
#include <exception>
#include <mpi.h>
#include "xmlUtils.h"

template<class Problem>
const char * rateMix<Problem>::name = "rateMix";

template<class Problem>
rateMix<Problem>::Param::Param(xmlNode *root) : weight(0.)
{
    xmlNode *section = getSectionByName(root, "rateMix");
    factor = getPropDouble(section, "factor");
    double memLength = 0.0;
    try {
        memLength = getPropDouble(section, "memLeght");
    } catch (exception &e) {
        
    }
    if (0.0 != memLength) {
        weight = std::exp(-7.0/memLength); // exp(-7) ~= 0.001
    }
}

template<class Problem>
rateMix<Problem>::rateMix(Problem &problem, const MPIState &mpiState,
                                  unirandom &rand, const Param &param):
        mix(problem, mpiState, rand), factor(param.factor), weight(param.weight),
        interval(1), tau_count(0), w(0.), v(0.), mpi(mpiState)
{
    adoptArray = new int[mpi.nnodes];
}

template<class Problem>
rateMix<Problem>::~rateMix()
{
    delete []adoptArray;
}

template<class Problem>
mixState rateMix<Problem>::Mix(aState &state)
{
    ++tau_count;
    if ((tau_count % interval) != 0)
        return mixState();
    int i,p, nadopt=0;
//    const double C = 1.38629436112; // 2ln(2)
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
//    double adoptRate = (double)(nadopt) / mpi.nnodes;
    w = 1.0 + w * weight;
    v = (double)nadopt + v * weight;
    interval = lround(factor * v/w);
    if (interval < 1)
        interval = 1;
    debugOut << state.step_cnt << " " << state.s << " "
             << nadopt << " " << interval << std::endl;
    return mixState(p);

}


#endif	/* RATEMIX_HPP */

