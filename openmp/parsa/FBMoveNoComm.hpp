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
        unirandom& in_rnd, xmlNode* root, const MPIState& mpiState) :
        feedbackMove<Problem>(in_problem, in_rnd, root), mpi(mpiState)
{
    theta_buf = new double[this->nparams];
    int win_size = this->nparams*sizeof(double);
    MPI_Win_create(theta_buf, win_size, win_size, MPI_INFO_NULL, mpi.comm,
                   &theta_win);
}

template<class Problem>
FBMoveNoComm<Problem>::~FBMoveNoComm()
{
    MPI_Win_free(&theta_win);
    delete []theta_buf;
}

template<class Problem>
void FBMoveNoComm<Problem>::processMix(const mixState& ms,
        const aState& state)
{
    for (int i = 0; i < this->nparams; ++ i) {
        theta_buf[i] = this->theta_bars[i];
    }
    MPI_Win_post(mpi.group, MPI_MODE_NOPUT, theta_win);
    MPI_Win_start(mpi.group, 0, theta_win);
    if (ms.doesMix()) {
        this->energy = state.energy;
        MPI_Get(this->theta_bars, this->nparams,MPI_DOUBLE, ms.adoptList[0], 0,
                this->nparams, MPI_DOUBLE, theta_win);
    }
    MPI_Win_complete(theta_win);
    MPI_Win_wait(theta_win);

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
