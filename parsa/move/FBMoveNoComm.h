/*
 * FBMoveNoComm.h
 *
 *  Created on: Jun 26, 2013
 *      Author: zhlou
 */

#ifndef FBMOVENOCOMM_H_
#define FBMOVENOCOMM_H_

/*
 * A parallel move control with no communication except for after adopting other
 * processor's state, which then adopt that processor's move control as well.
 */

#include <mpi.h>
#include "move/feedbackMove.h"
#include "MPIState.h"
#include "mixState.h"

template<class Problem>
class FBMoveNoComm: public feedbackMove<Problem>
{
public:
    FBMoveNoComm(Problem &in_problem, unirandom& in_rnd, xmlNode *root,
                 const MPIState &miState);
    ~FBMoveNoComm();
    int getWinner();
    void processMix(const mixState &ms, const aState &state);
    static const char *name;
private:
    const MPIState &mpi;
    double *theta_buf;
    MPI_Win theta_win;
};



#include "FBMoveNoComm.hpp"

#endif /* FBMOVENOCOMM_H_ */
