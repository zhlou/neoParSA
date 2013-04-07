/*
 * pannealer.h
 *
 *  Created on: Mar 5, 2013
 *      Author: zhlou
 */

#ifndef PANNEALER_H_
#define PANNEALER_H_
#include "annealer.h"
#include "MPIState.h"
template <class Problem, class Schedule, class FrozenCnd, template<class> class Move,
          template<class> class PopBased>
class pannealer : public annealer<Problem, Schedule, FrozenCnd, Move>
{
public:
    pannealer(Problem &problem, unirandom * const in_rand,
              typename FrozenCnd::Param frozenParam, xmlNode *root,
              const MPIState &mpiState);
    virtual ~pannealer();
    int getWinner();
protected:
    const MPIState &mpi;
    PopBased<Problem> pop;
    virtual void updateSegment(aState &state);
    virtual void writeMethodText(xmlNode *method);

};

#include "pannealer.hpp"

#endif /* PANNEALER_H_ */
