/*
 * parallelFBMove.h
 *
 *  Created on: Jan 28, 2013
 *      Author: zhlou
 */

#ifndef PARALLELFBMOVE_H_
#define PARALLELFBMOVE_H_

#include "feedbackMove.h"

template<class Problem, class Debug>
class parallelFBMove: public feedbackMove<Problem, Debug>
{
public:
    parallelFBMove(Problem &in_problem, xmlNode *root, MPI_Comm thecomm,
            int in_nnodes, int in_rank);
    int getWinner();
protected:
    void collectMoveStats();
private:
    MPI_Comm comm;
    int nnodes;
    int rank;

};


#include "parallelFBMove.hpp"

#endif /* PARALLELFBMOVE_H_ */
