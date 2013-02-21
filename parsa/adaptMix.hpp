/*
 * adaptMix.hpp
 *
 *  Created on: Jan 31, 2013
 *      Author: zhlou
 */
#include <limits>
#include <cmath>
#include <iostream>
using namespace std;
template<class Problem>
adaptMix<Problem>::adaptMix(Problem& in_problem, const MPIState& mpiState,
        xmlNode *docroot) :
        problem(in_problem), mpi(mpiState), rnd(mpi.rank), root(docroot)
{
    xmlNode *section = getSectionByName(root, "mix");
    adaptCoef = getPropDouble(section, "adaptcoef");

    energy_tab = new double[mpi.nnodes];
    prob_tab = new double[mpi.nnodes];
    buf_size = problem.getStateSize();
    MPI_Alloc_mem(buf_size, MPI_INFO_NULL, &state_buf);
    MPI_Win_create(state_buf, buf_size, buf_size, MPI_INFO_NULL, mpi.comm,
            &state_win);
}

template<class Problem>
adaptMix<Problem>::~adaptMix()
{
    MPI_Win_free(&state_win);
    MPI_Free_mem(state_buf);
    delete[] energy_tab;
}

template<class Problem>
double adaptMix<Problem>::Mix(aState& state)
{
    double energy = state.energy;
    int i;
    double prob, norm = 0;
    MPI_Allgather(&energy, 1, MPI_DOUBLE, energy_tab, 1, MPI_DOUBLE,
            mpi.comm);
    //cout << "energy_tab@" << state.step_cnt;
    for (i = 0; i < mpi.nnodes; ++i) {
        prob = exp((energy - energy_tab[i]) * state.s);
        if (prob < numeric_limits<double>::min())
            prob = numeric_limits<double>::min();
        if (prob > numeric_limits<double>::max()/mpi.nnodes)
            prob = numeric_limits<double>::max()/mpi.nnodes;
        norm += prob_tab[i] = prob;
        //cout << " " << prob_tab[i];
    }
    //cout << " norm = " << norm << endl;

    double rand = rnd.random();
    if (rand > adaptCoef * mpi.nnodes / norm) {
        rand = rnd.random();
        double psum = 0.;
        for (i = 0; i < mpi.nnodes; ++i) {
            psum += (prob_tab[i]/norm);
            if (psum > rand)
                break;
        }
        cout << "rank " << mpi.rank << " adopts state from rank " << i << endl;
        adoptState(i);
        return (state.energy = problem.get_score());
    } else {
        adoptState(mpi.rank);
        return state.energy;
    }
}
template<class Problem>
void adaptMix<Problem>::adoptState(int Id)
{
    problem.serialize(state_buf);
    MPI_Win_post(mpi.group, MPI_MODE_NOPUT, state_win);
    MPI_Win_start(mpi.group, 0, state_win);
    if (Id != mpi.rank)
        MPI_Get(state_buf, buf_size, MPI_BYTE, Id, 0, buf_size, MPI_BYTE,
                state_win);
    MPI_Win_complete(state_win);
    MPI_Win_wait(state_win);
    if (Id != mpi.rank) {

        problem.deserialize(state_buf);
    }

}
