/*
 * parallelFBMove.hpp
 *
 *  Created on: Jan 28, 2013
 *      Author: zhlou
 */

template<class Problem>
parallelFBMove<Problem>::parallelFBMove(Problem& in_problem, xmlNode* root,
        MPI_Comm thecomm, int in_nnodes, int in_rank) :
        feedbackMove<Problem>(in_problem, root), comm(thecomm), nnodes(
                in_nnodes), rank(in_rank)
{
}

template<class Problem>
void parallelFBMove<Problem>::collectMoveStats()
{
    MPI_Allreduce(MPI_IN_PLACE, success, nparams, MPI_LONG, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, moves, nparams, MPI_LONG, MPI_SUM, comm);
}
