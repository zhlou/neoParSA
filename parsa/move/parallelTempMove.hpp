#ifndef PARALLELTEMPMOVE_HPP
#define PARALLELTEMPMOVE_HPP

template<class Problem>
const char *parallelTempMove<Problem>::name = "parallelTempMove";

template<class Problem>
parallelTempMove<Problem>::parallelTempMove(Problem &in_problem,
        unirandom &in_rnd, const ptree &root, const MPIState &mpiState) :
        tempFeedback<Problem>(in_problem, in_rnd, root), mpi(mpiState)
{
    bool use_rank = root.get<bool>("tempMove.<xmlattr>.use_rank", false);
    if (use_rank)
        this->index = (mpi.rank % (this->nparams));

}

template<class Problem>
parallelTempMove<Problem>::~parallelTempMove()
{

}

template<class Problem>
int parallelTempMove<Problem>::getWinner()
{
    struct
    {
        double energy;
        int rank;
    } doubleint;
    doubleint.energy = tempFeedback<Problem>::energy;
    doubleint.rank = mpi.rank;
    MPI_Allreduce(MPI_IN_PLACE, &doubleint, 1, MPI_DOUBLE_INT, MPI_MINLOC,
            mpi.comm);
    return doubleint.rank;
}

template<class Problem>
void parallelTempMove<Problem>::collectMoveStats()
{
    MPI_Allreduce(MPI_IN_PLACE, &tempFeedback<Problem>::success[0],
            tempFeedback<Problem>::nparams, MPI_LONG, MPI_SUM, mpi.comm);
    MPI_Allreduce(MPI_IN_PLACE, &tempFeedback<Problem>::moves[0],
            tempFeedback<Problem>::nparams, MPI_LONG, MPI_SUM, mpi.comm);
}

template<class Problem>
void parallelTempMove<Problem>::readState(const ptree &root)
{
    tempFeedback<Problem>::readState(root);
    MPI_Allreduce(MPI_IN_PLACE, &tempFeedback<Problem>::theta_bars[0],
            tempFeedback<Problem>::nparams,
            MPI_DOUBLE, MPI_SUM, mpi.comm);
    for (int i = 0; i < tempFeedback<Problem>::nparams; ++i) {
        tempFeedback<Problem>::theta_bars[i] /= mpi.nnodes;
    }
}


#endif
