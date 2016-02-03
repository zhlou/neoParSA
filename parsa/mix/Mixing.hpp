/*
 * Mixing.hpp
 *
 *  Created on: Mar 22, 2013
 *      Author: zhlou
 */
#include <cmath>
#include <limits>
#include <mpi.h>
using namespace std;

template<typename Problem>
Mixing<Problem>::Mixing(Problem & in_problem, const MPIState &mpiState,
                  unirandom& in_rnd) :
        problem(in_problem), mpi(mpiState), rnd(in_rnd)
{
    energy_tab = new double[mpi.nnodes];
    prob_tab = new double[mpi.nnodes];
    buf_size = problem.getStateSize();
    MPI_Info info_no_locks;
    MPI_Info_create(&info_no_locks);
    MPI_Info_set(info_no_locks, (char *)"no_locks",(char *)"true");
    MPI_Alloc_mem(buf_size, MPI_INFO_NULL, &state_buf);
    MPI_Win_create(state_buf, buf_size, buf_size, info_no_locks, mpi.comm,
            &state_win);
    recv_buf = calloc(1, buf_size);
    norm = 0;
    adoptArray = new int[mpi.nnodes];
}

template<typename Problem>
Mixing<Problem>::~Mixing()
{
    MPI_Win_free(&state_win);
    MPI_Free_mem(state_buf);
    free(recv_buf);
    delete[] energy_tab;
    delete[] prob_tab;
    delete[] adoptArray;
}

template<typename Problem>
void Mixing<Problem>::calProbTab(const aState &state)
{
    double prob;
    double energy = state.energy;
    int i;
    norm = 0;
    MPI_Allgather(&energy, 1, MPI_DOUBLE, energy_tab, 1, MPI_DOUBLE,
            mpi.comm);
    debugOut << mpi.rank << " @ " << state.step_cnt << " " << state.s << " ";
    for (i = 0; i < mpi.nnodes; ++i) {
        prob = exp((energy - energy_tab[i]) * state.s);
        if (prob < numeric_limits<double>::min())
            prob = numeric_limits<double>::min();
        if (prob > numeric_limits<double>::max()/mpi.nnodes)
            prob = numeric_limits<double>::max()/mpi.nnodes;
        norm += prob_tab[i] = prob;
        //debugOut << " " << prob_tab[i];
        //debugOut << " " << energy_tab[i];
    }
    // debugOut << endl;

}

template<typename Problem>
double Mixing<Problem>::adoptState(int Id)
{
    problem.serialize(state_buf);
    MPI_Win_post(mpi.group, MPI_MODE_NOPUT, state_win);
    MPI_Win_start(mpi.group, 0, state_win);
    if (Id < 0 || Id >= mpi.nnodes) { // should not happen
        std::cerr<< "adoptState " << mpi.rank <<": Invalid Id " << Id <<endl;
        Id = mpi.rank;
    }
    if (Id != mpi.rank)
        MPI_Get(recv_buf, buf_size, MPI_BYTE, Id, 0, buf_size, MPI_BYTE,
                state_win);
    MPI_Win_complete(state_win);
    MPI_Win_wait(state_win);
    if (Id != mpi.rank) {

        problem.deserialize(recv_buf);
    }
    return problem.get_score();
}

template<typename Problem>
int Mixing<Problem>::getPartner() const
{
    double rand = rnd.random();
    double psum = 0.;
    int i;
    for (i = 0; i < mpi.nnodes; ++i) {
        psum += (prob_tab[i]/norm);
        if (psum > rand)
            break;
    }
    if (i == mpi.nnodes) {
        std::cerr << "getPartner " << mpi.rank <<": psum = " << psum
                  << " rand = " << rand << endl;
        i = mpi.rank;
    }
    debugOut << " Adopt " << i << endl;
    return i;
}

template<typename Problem>
int Mixing<Problem>::getBest() const
{
    int curBest = 0;
    double bestEnergy = energy_tab[0];
    for (int i = 1; i < mpi.nnodes; ++i) {
        if (energy_tab[i] < bestEnergy) {
            bestEnergy = energy_tab[i];
            curBest = i;
        }
    }
    return curBest;
}

template<typename Problem>
double Mixing<Problem>::getEnergyVar(size_t ddof) const
{
    // Calculating variance using two-pass with correction a la A. Bjorck
    const int N = mpi.nnodes;
    if (N < 2)
        return 0.;

    int i = 0;
    
    double m=0.;
    for (i = 0; i < N; ++i) {
        m += energy_tab[i];
    }
    m = m/N;
    double s=0.;
    double e=0.;
    double delta;
    for (i = 0; i < N; ++i) {
        delta = energy_tab[i] - m;
        s+= delta*delta;
        e+= delta;
    }
    s -= (e*e)/N;
    return s/(N-ddof);
}

template<typename Problem>
int Mixing<Problem>::getNAdopt(int id)
{
    int i,nadopt = 0;
    for (i = 0; i < mpi.nnodes; ++i)
        adoptArray[i] = 0;
    adoptArray[id] = 1;
    MPI_Allreduce(MPI_IN_PLACE, adoptArray, mpi.nnodes, MPI_INT, MPI_LOR, mpi.comm);
    for (i = 0; i < mpi.nnodes; ++i)
        nadopt += adoptArray[i];
    return nadopt;
}