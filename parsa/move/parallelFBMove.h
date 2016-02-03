/*
 * parallelFBMove.h
 *
 *  Created on: Jan 28, 2013
 *      Author: zhlou
 */

#ifndef PARALLELFBMOVE_H_
#define PARALLELFBMOVE_H_

#include <mpi.h>
#include "move/feedbackMove.h"
#include "MPIState.h"
#include "mix/mixState.h"

template<class Problem>
class parallelFBMove: public feedbackMove<Problem>
{
public:
    parallelFBMove(Problem &in_problem, unirandom& in_rnd, xmlNode *root,
                   const MPIState &mpiState);
    parallelFBMove(Problem &in_problem, unirandom& in_rnd, const ptree &root,
                const MPIState &mpiState);
    ~parallelFBMove();
    int getWinner();
    // void doMix(aState &state);
    void processMix(const mixState &ms, const aState &state)
    {if (ms.doesMix()) this->energy = state.energy;}
    static const char *name;
    void readState(xmlNodePtr docroot); // it overrides the base class. bite me
    void readState(const ptree &root); // and you may want to bite me twice
protected:
    void collectMoveStats();
private:
    const MPIState &mpi;

};


#include "parallelFBMove.hpp"

#endif /* PARALLELFBMOVE_H_ */
