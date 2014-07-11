/*
 * pannealer.hpp
 *
 *  Created on: Mar 5, 2013
 *      Author: zhlou
 */
#include "mixState.h"
template <class Problem, class Schedule, class FrozenCnd,
          template<class> class Move, template<class> class PopBased>
pannealer<Problem, Schedule, FrozenCnd, Move,
          PopBased>::pannealer(Problem &problem, unirandom& in_rand,
                               typename Schedule::Param scheParam,
                               typename FrozenCnd::Param frozenParam,
                               xmlNode *root, const MPIState &mpiState) :
          annealer<Problem, Schedule, FrozenCnd, Move>::annealer(in_rand, root),
          mpi(mpiState), pop(problem, mpiState, in_rand, root)
{
    // this pointer is necessary because otherwise the lookup to parent members
    // may fail depending on compilers
    this->cooling = new Schedule(scheParam, mpiState);
    this->move = new Move<Problem>(problem, in_rand, root, mpiState);
    this->state.energy = this->move->get_score();
    this->frozen = new FrozenCnd(frozenParam, mpiState);
}

template<class Problem, class Schedule, class FrozenCnd,
        template<class > class Move, template<class> class PopBased>
pannealer<Problem, Schedule, FrozenCnd, Move, PopBased>::~pannealer()
{

}

template<class Problem, class Schedule, class FrozenCnd,
        template<class> class Move, template<class> class PopBased>
int pannealer<Problem, Schedule, FrozenCnd, Move, PopBased>::getWinner()
{
    struct
    {
        double energy;
        int rank;
    } doubleint;
    doubleint.energy = (this->state).energy;
    doubleint.rank = mpi.rank;
    MPI_Allreduce(MPI_IN_PLACE, &doubleint, 1, MPI_DOUBLE_INT, MPI_MINLOC,
            mpi.comm);
    return doubleint.rank;
}

template <class Problem, class Schedule, class FrozenCnd,
          template<class> class Move, template<class> class PopBased>
void pannealer<Problem, Schedule, FrozenCnd, Move, PopBased>::updateSegment(aState &state)
{

    annealer<Problem, Schedule, FrozenCnd, Move>::updateSegment(state);
    mixState ms = pop.Mix(state);
    this->move->processMix(ms, state);

}

template <class Problem, class Schedule, class FrozenCnd,
          template<class> class Move, template<class> class PopBased>
void pannealer<Problem, Schedule, FrozenCnd, Move,
               PopBased>::writeMethodText(xmlNode* method)
{
    annealer<Problem, Schedule, FrozenCnd, Move>::writeMethodText(method);
    xmlNewProp(method, (const xmlChar*)"mixing",
               (const xmlChar*)PopBased<Problem>::name);

}
