/*
 * parallelFBMove.h
 *
 *  Created on: Jan 28, 2013
 *      Author: zhlou
 */

#ifndef PARALLELCAUCHYMOVE_H_
#define PARALLELCAUCHYMOVE_H_

#include <mpi.h>
#include "move/cauchyMove.h"
#include "MPIState.h"
#include "mix/mixState.h"

template<class Problem>
class parallelCauchyMove: public cauchyMove<Problem>
{
public:
    parallelCauchyMove(Problem &in_problem, unirandom& in_rnd, const ptree &root,
                const MPIState &mpiState);
    ~parallelCauchyMove();
    int getWinner();
    // void doMix(aState &state);
    void processMix(const mixState &ms, const aState &state)
    {if (ms.doesMix()) this->energy = state.energy;}
    static const char *name;
    void readState(const ptree &root); // and you may want to bite me twice
protected:
    void collectMoveStats();
private:
    const MPIState &mpi;

};


#include "parallelCauchyMove.hpp"

#endif /* PARALLELCAUCHYMOVE_H_ */
