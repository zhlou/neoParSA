/*
 * parallelCauchyMove.hpp
 *
 *  Created on: Jun 27, 2016
 *      Author: zhlou
 */
template <class Problem>
const char *parallelCauchyMove<Problem>::name = "parallelCauchyMove";

template<class Problem>
parallelCauchyMove<Problem>::parallelCauchyMove(Problem &in_problem,
        unirandom &in_rnd, const ptree &root, const MPIState &mpiState) :
        cauchyMove<Problem>(in_problem, in_rnd, root), mpi(mpiState)
{
}

template<class Problem>
int parallelCauchyMove<Problem>::getWinner()
{
    struct
    {
        double energy;
        int rank;
    } doubleint;
    doubleint.energy = cauchyMove<Problem>::energy;
    doubleint.rank = mpi.rank;
    MPI_Allreduce(MPI_IN_PLACE, &doubleint, 1, MPI_DOUBLE_INT, MPI_MINLOC,
            mpi.comm);
    return doubleint.rank;
}

template<class Problem>
parallelCauchyMove<Problem>::~parallelCauchyMove()
{

}

template<class Problem>
void parallelCauchyMove<Problem>::collectMoveStats()
{
    MPI_Allreduce(MPI_IN_PLACE, cauchyMove<Problem>::success,
            cauchyMove<Problem>::nparams, MPI_LONG, MPI_SUM, mpi.comm);
    MPI_Allreduce(MPI_IN_PLACE, cauchyMove<Problem>::moves,
            cauchyMove<Problem>::nparams, MPI_LONG, MPI_SUM, mpi.comm);
}

template<class Problem>
void parallelCauchyMove<Problem>::readState(const ptree &root)
{
    cauchyMove<Problem>::readState(root);
    MPI_Allreduce(MPI_IN_PLACE, cauchyMove<Problem>::theta_bars,
            cauchyMove<Problem>::nparams,
            MPI_DOUBLE, MPI_SUM, mpi.comm);
    for (int i = 0; i < cauchyMove<Problem>::nparams; ++i) {
        cauchyMove<Problem>::theta_bars[i] /= mpi.nnodes;
    }
}
