/*
 * FBNCTheta.h
 *
 *  Created on: Sep 16, 2013
 *      Author: zhlou
 */

#ifndef FBNCTHETA_H_
#define FBNCTHETA_H_
// Feedback move, no communication,
// with theta_bar adjustment
template<class Problem>
class FBNCTheta: public feedbackMove<Problem>
{
private:
    const MPIState &mpi;
    double *theta_buf;
    MPI_Win state_win;
    void *state_buf;
    void *local_buf;
    const int buf_size;
public:
    FBNCTheta(Problem &in_problem, unirandom * const in_rnd, xmlNode *root,
              const MPIState &mpiState);
    ~FBNCTheta();
    int getWinner();
    void processMix(const mixState &ms, const aState &state);
    static const char * name;
};

#include "FBNCTheta.hpp"

#endif /* FBNCTHETA_H_ */
