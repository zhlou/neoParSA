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

template<class Problem, class Debug, template <class> class PopBased>
class parallelFBMove: public feedbackMove<Problem, Debug>
{
public:
    parallelFBMove(Problem &in_problem, xmlNode *root, const MPIState &mpiState);
    ~parallelFBMove();
    int getWinner();
    void DoMix(aState &state);
protected:
    void collectMoveStats();
private:
    const MPIState &mpi;

    PopBased<Problem> pop;



};


#include "parallelFBMove.hpp"

#endif /* PARALLELFBMOVE_H_ */
