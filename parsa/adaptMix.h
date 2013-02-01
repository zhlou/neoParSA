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

/*
 * In addition to the requirements listed in feedbackMove.h, using adaptMix
 * requires the Problem class additional interfaces
 * int getStateSize(); // returns the byte count of the state
 * void serialize(void *buf); // serialize its state to buf
 * void deserialize(void *buf); // inflate buf to a new state. calculation
 *                              // of new score is not required here.
 */

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
