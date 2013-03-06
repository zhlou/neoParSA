/*
 * pannealer.hpp
 *
 *  Created on: Mar 5, 2013
 *      Author: zhlou
 */

template <class Problem, class Schedule, template<class> class Move,
          template<class> class PopBased>
pannealer<Problem, Schedule, Move, PopBased>::pannealer(Problem &problem,
          unirandom * const in_rand, xmlNode *root, const MPIState &mpiState) :
          annealer<Problem, Schedule, Move>(problem, in_rand, root),
          mpi(mpiState), pop(problem, mpiState, root)
{

}

template <class Problem, class Schedule, template<class> class Move,
          template<class> class PopBased>
void pannealer<Problem, Schedule, Move, PopBased>::updateSegment(aState &state)
{
    annealer<Problem, Schedule, Move>::updateSegment(state);
    pop.Mix(state);
}
