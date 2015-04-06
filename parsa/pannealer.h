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
    pannealer(Problem &problem, unirandom& in_rand,
              typename Schedule::Param scheParam,
              typename FrozenCnd::Param frozenParam,
              const typename PopBased<Problem>::Param &popParam,
              xmlNode *root,
              const MPIState &mpiState);
    virtual ~pannealer();
    int getWinner();
    void setMixLog(debugStatus st, const char* outname=NULL)
    {pop.setDebug(st, outname);}
#ifdef USE_BOOST
    virtual void ptreeGetResult(ptree &pt);
#endif
protected:
    const MPIState &mpi;
    PopBased<Problem> pop;
    virtual void updateStats(aState &state);
    virtual void writeMethodText(xmlNode *method);

};

#include "pannealer.hpp"

#endif /* PANNEALER_H_ */
