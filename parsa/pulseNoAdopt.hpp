/*
 * pulseNoAdopt.hpp
 *
 *  Created on: Sep 16, 2013
 *      Author: zhlou
 */
#include <mpi.h>
#include <cmath>
#include "xmlUtils.h"

template <class Problem>
const char * pulseNoAdopt<Problem>::name = "pulseNoAdopt";

template<class Problem>
pulseNoAdopt<Problem>::Param::Param(xmlNode *docroot)
{
    xmlNode *section = getSectionByName(docroot, "pulseBcast");
    if (section == NULL) {
        throw std::runtime_error(std::string("Error: cannot find section pulseBcast"));
    }
    frequency = getPropInt(section, "frequency");

}

template<class Problem>
pulseNoAdopt<Problem>::pulseNoAdopt(Problem& in_problem, const MPIState& mpiState,
                                unirandom& in_rnd, const Param &param) :
        problem(in_problem), mpi(mpiState), rnd(in_rnd), counter(0),
        frequency(param.frequency)
{

    //MPI_Alloc_mem(buf_size,MPI_INFO_NULL, &state_buf);

}

template<class Problem>
pulseNoAdopt<Problem>::~pulseNoAdopt()
{
    //MPI_Free_mem(state_buf);
}

template<class Problem>
mixState pulseNoAdopt<Problem>::Mix(aState& state)
{
    if (counter % frequency) {
        ++ counter;
        return mixState();
    }
    int root = (counter / frequency) % mpi.nnodes;
    ++ counter;
    double energy;
    if (mpi.rank == root) {
        energy = state.energy;
        //problem.serialize(state_buf);
    }
    MPI_Bcast(&energy, 1, MPI_DOUBLE, root, mpi.comm);
    // MPI_Bcast(state_buf, buf_size, MPI_BYTE, root, mpi.comm);
    if (mpi.rank != root) {
        // The incoming state is treated as current state
        // and real current state is treated as proposed state.
        // Thus we do one step of Metropolis algorithm
        double delta = energy - state.energy;
        double crit = 1.0 - std::exp(state.s * delta); // no negative sign here
        double randval = rnd.random();
        if ((delta <= 0.0) && (crit > randval)) { // adopt the incoming state
            // problem.deserialize(state_buf);
            // state.energy = energy;
            return mixState(root);
        }
    }

    return mixState();
}
