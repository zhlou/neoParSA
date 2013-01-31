/*
 * adaptMix.h
 *
 *  Created on: Jan 31, 2013
 *      Author: zhlou
 */

#ifndef ADAPTMIX_H_
#define ADAPTMIX_H_

#include "MPIState.h"
#include "unirandom.h"
#include <libxml/tree.h>


template <class Problem>
class adaptMix
{
public:
    adaptMix(Problem &in_problem, const MPIState &mpiState, xmlNode *docroot);
    ~adaptMix();
    double Mix(aState &state);
private:
    Problem &problem;
    unirandom rnd;
    const MPIState &mpi;
    MPI_Win state_win;
    void *state_buf;
    int buf_size;
    double *energy_tab;
    double *prob_tab;
    double adaptCoef;
    xmlNode *root;
    void adoptState(int Id);

};



#include "adaptMix.hpp"

#endif /* ADAPTMIX_H_ */
