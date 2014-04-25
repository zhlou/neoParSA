/*
 * Mixing.h
 *
 *  Created on: Mar 22, 2013
 *      Author: zhlou
 */

#ifndef MIX_H_
#define MIX_H_

#include "unirandom.h"
#include "MPIState.h"

template<typename Problem>
class Mixing
{
private:
    Problem &problem;
    const MPIState &mpi;
    unirandom * const rnd;

    MPI_Win state_win;
    void *state_buf;
    int buf_size;
    double *energy_tab;
    double *prob_tab;
    double norm;
public:
    Mixing(Problem &problem, const MPIState &mpiState, unirandom * const in_rnd);
    ~Mixing();
    void calProbTab(const aState &state);
    double getNorm() const {return norm;}
    double adoptState(int id);
    double adoptState(){return adoptState(mpi.rank);}
    int getPartner() const;
    mutable dynDebug debugOut;
};

#include "Mixing.hpp"

#endif /* MIX_H_ */
