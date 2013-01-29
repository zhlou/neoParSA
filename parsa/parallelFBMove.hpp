/*
 * parallelFBMove.hpp
 *
 *  Created on: Jan 28, 2013
 *      Author: zhlou
 */

template<class Problem, class Debug>
parallelFBMove<Problem, Debug>::parallelFBMove(Problem& in_problem, xmlNode* root,
        MPI_Comm thecomm, int in_nnodes, int in_rank) :
        feedbackMove<Problem, Debug>(in_problem, root), comm(thecomm), nnodes(
                in_nnodes), rank(in_rank)
{
}

template<class Problem, class Debug>
int parallelFBMove<Problem, Debug>::getWinner()
{
    struct {
        double energy;
        int rank;
    } doubleint;
    doubleint.energy = feedbackMove<Problem, Debug>::energy;
    doubleint.rank = rank;
    MPI_Allreduce(MPI_IN_PLACE, &doubleint, 1, MPI_DOUBLE_INT, MPI_MINLOC,
            comm);
    return doubleint.rank;
}

template<class Problem, class Debug>
void parallelFBMove<Problem, Debug>::collectMoveStats()
{
    MPI_Allreduce(MPI_IN_PLACE, feedbackMove<Problem, Debug>::success,
            feedbackMove<Problem, Debug>::nparams, MPI_LONG, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, feedbackMove<Problem, Debug>::moves,
            feedbackMove<Problem, Debug>::nparams, MPI_LONG, MPI_SUM, comm);
}
