/*
 * parallelFBMove.h
 *
 *  Created on: Jan 28, 2013
 *      Author: zhlou
 */

#ifndef PARALLELFBMOVE_H_
#define PARALLELFBMOVE_H_

#include <mpi.h>
#include "feedbackMove.h"
#include "MPIState.h"
#include "mixState.h"

template<class Problem>
class parallelFBMove: public feedbackMove<Problem>
{
public:
    parallelFBMove(Problem &in_problem, unirandom * const in_rnd, xmlNode *root,
                   const MPIState &mpiState);
    ~parallelFBMove();
    int getWinner();
    // void doMix(aState &state);
    void processMix(const mixState &ms, const aState &state)
    {if (ms.doesMix()) this->energy = state.energy;}
    static const char *name;
protected:
    void collectMoveStats();
private:
    const MPIState &mpi;

};


#include "parallelFBMove.hpp"

#endif /* PARALLELFBMOVE_H_ */
