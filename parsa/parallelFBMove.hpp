/*
 * parallelFBMove.hpp
 *
 *  Created on: Jan 28, 2013
 *      Author: zhlou
 */

template<class Problem, class Debug, template<class > class PopBased>
parallelFBMove<Problem, Debug, PopBased>::parallelFBMove(Problem& in_problem,
        unirandom &in_rnd, xmlNode* root, const MPIState &mpiState) :
        feedbackMove<Problem, Debug>(in_problem, in_rnd, root), mpi(mpiState),
        pop(in_problem, mpiState, root)
{
}

template<class Problem, class Debug, template<class > class PopBased>
int parallelFBMove<Problem, Debug, PopBased>::getWinner()
{
    struct
    {
        double energy;
        int rank;
    } doubleint;
    doubleint.energy = feedbackMove<Problem, Debug>::energy;
    doubleint.rank = mpi.rank;
    MPI_Allreduce(MPI_IN_PLACE, &doubleint, 1, MPI_DOUBLE_INT, MPI_MINLOC,
            mpi.comm);
    return doubleint.rank;
}

template<class Problem, class Debug, template<class > class PopBased>
void parallelFBMove<Problem, Debug, PopBased>::doMix(aState &state)
{
    //cout << "Run adapt mix" << endl;
    feedbackMove<Problem, Debug>::energy = pop.Mix(state);
}

template<class Problem, class Debug, template<class > class PopBased>
parallelFBMove<Problem, Debug, PopBased>::~parallelFBMove()
{

}

template<class Problem, class Debug, template<class > class PopBased>
void parallelFBMove<Problem, Debug, PopBased>::collectMoveStats()
{
    MPI_Allreduce(MPI_IN_PLACE, feedbackMove<Problem, Debug>::success,
            feedbackMove<Problem, Debug>::nparams, MPI_LONG, MPI_SUM, mpi.comm);
    MPI_Allreduce(MPI_IN_PLACE, feedbackMove<Problem, Debug>::moves,
            feedbackMove<Problem, Debug>::nparams, MPI_LONG, MPI_SUM, mpi.comm);
}
