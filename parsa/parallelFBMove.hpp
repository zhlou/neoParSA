/*
 * parallelFBMove.hpp
 *
 *  Created on: Jan 28, 2013
 *      Author: zhlou
 */

template<class Problem, template<class > class PopBased>
parallelFBMove<Problem, PopBased>::parallelFBMove(Problem& in_problem,
        unirandom &in_rnd, xmlNode* root,
        const MPIState &mpiState) :
        feedbackMove<Problem>(in_problem, in_rnd, root), mpi(mpiState),
        pop(in_problem, mpiState, root)
{
}

template<class Problem, template<class > class PopBased>
int parallelFBMove<Problem, PopBased>::getWinner()
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

template<class Problem, template<class > class PopBased>
void parallelFBMove<Problem, PopBased>::doMix(aState &state)
{
    //cout << "Run adapt mix" << endl;
    feedbackMove<Problem>::energy = pop.Mix(state);
}

template<class Problem, template<class > class PopBased>
parallelFBMove<Problem, PopBased>::~parallelFBMove()
{

}

template<class Problem, template<class > class PopBased>
void parallelFBMove<Problem, PopBased>::collectMoveStats()
{
    MPI_Allreduce(MPI_IN_PLACE, feedbackMove<Problem>::success,
            feedbackMove<Problem>::nparams, MPI_LONG, MPI_SUM, mpi.comm);
    MPI_Allreduce(MPI_IN_PLACE, feedbackMove<Problem>::moves,
            feedbackMove<Problem>::nparams, MPI_LONG, MPI_SUM, mpi.comm);
}
