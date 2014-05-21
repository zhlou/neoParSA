/*
 * FBNCTheta.hpp
 *
 *  Created on: Sep 17, 2013
 *      Author: zhlou
 */

template <class Problem>
const char *FBNCTheta<Problem>::name = "FBNCTheta";

template <class Problem>
FBNCTheta<Problem>::FBNCTheta(Problem& in_problem,
unirandom& in_rnd, xmlNode* root, const MPIState& mpiState) :
        feedbackMove<Problem>(in_problem, in_rnd, root), mpi(mpiState),
        buf_size(in_problem.getStateSize())
{
    theta_buf = new double[this->nparams];
    MPI_Alloc_mem(buf_size,MPI_INFO_NULL, &state_buf);
    MPI_Win_create(state_buf, buf_size, buf_size, MPI_INFO_NULL, mpi.comm,
            &state_win);
    local_buf = malloc(buf_size);
}

template <class Problem>
FBNCTheta<Problem>::~FBNCTheta()
{
    MPI_Win_free(&state_win);
    MPI_Free_mem(state_buf);
    free(local_buf);
}

template <class Problem>
void FBNCTheta<Problem>::processMix(const mixState& ms, const aState& state)
{
    this->problem.serialize(state_buf);
    MPI_Win_post(mpi.group, MPI_MODE_NOPUT, state_win);
    MPI_Win_start(mpi.group, 0, state_win);
    if(ms.doesMix()) {
        MPI_Get(local_buf, buf_size, MPI_BYTE, ms.adoptList[0],0,
                buf_size, MPI_BYTE, state_win);
        this->problem.state2theta(local_buf,theta_buf);
        this->energy = state.energy;
        for (int i=0;i<(this->nparams);++i) {
            this->theta_bars[i] = theta_buf[i];
        }
    }
    MPI_Win_complete(state_win);
    MPI_Win_wait(state_win);
}

template<class Problem>
int FBNCTheta<Problem>::getWinner()
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
