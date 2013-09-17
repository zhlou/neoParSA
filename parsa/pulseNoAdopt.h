/*
 * pulseNoAdopt.h
 *
 *  Created on: Sep 16, 2013
 *      Author: zhlou
 */

#ifndef PULSENOADOPT_H_
#define PULSENOADOPT_H_

#include "MPIState.h"
#include "unirandom.h"
#include "mixState.h"
#include <libxml/tree.h>

template <class Problem>
class pulseNoAdopt
{
private:
    Problem &problem;
    const MPIState &mpi;
    unirandom * const rnd;
    //void *state_buf;
    //const int buf_size;
    unsigned counter;
    unsigned frequency;
    dynDebug debugOut;
public:
    pulseNoAdopt(Problem &problem, const MPIState &mpiState,
                 unirandom * const in_rnd, xmlNode *docroot);
    ~pulseNoAdopt();
    mixState Mix(aState &state);
    static const char * name;

};

#include "pulseNoAdopt.hpp"

#endif /* PULSENOADOPT_H_ */
