/* 
 * File:   coupledAnnealer.h
 * Author: zhlou
 *
 * Created on January 27, 2015, 1:21 PM
 */

#ifndef COUPLEDANNEALER_H
#define	COUPLEDANNEALER_H

#include "annealer.h"
#include "MPIState.h"

template <class Problem, class Schedule, class FrozenCnd, template<class> class CoupledMove>
class coupledAnnealer : public annealer<Problem, Schedule, FrozenCnd, CoupledMove>
{
public:
    coupledAnnealer(Problem &problem, 
            unirandom& in_rand,
            typename Schedule::Param scheParam,
            typename FrozenCnd::Param frozenParam,
            const typename CoupledMove<Problem>::Param &moveParam,
            xmlNode *root,
            const MPIState &mpiState);
    virtual ~coupledAnnealer();
    int getWinner();
    void setMixLog(debugStatus st, const char* outname=NULL) {
        this->move->setMixLog(st, outname);
    }

protected:
    const MPIState &mpi;
    virtual void updateStats(aState &state);

        
        
};

#include "coupledAnnealer.hpp"

#endif	/* COUPLEDANNEALER_H */

