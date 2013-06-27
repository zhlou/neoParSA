/*
 * FBMoveNoComm.hpp
 *
 *  Created on: Jun 27, 2013
 *      Author: zhlou
 */

template <class Problem>
const char *FBMoveNoComm<Problem>::name = "FeedbackMoveNoComm";

template<class Problem>
FBMoveNoComm<Problem>::FBMoveNoComm(Problem& in_problem,
        unirandom* const in_rnd, xmlNode* root, const MPIState& miState) {
}

template<class Problem>
void FBMoveNoComm<Problem>::processMix(const mixState& ms,
        const aState& state) {
}

template<class Problem>
int FBMoveNoComm<Problem>::getWinner()
{
    struct
    {
        double energy;
        int rank;
    } doubleint;
    doubleint.energy = feedbackMove<Problem>::energy;
    doubleint.rank = mpi.rank;
    MPI_Allreduce(MPI_IN_PLACE, &doubleint, 1, MPI_DOUBLE_INT, MPI_MINLOC,
            mpi.comm);
    return doubleint.rank;
}
