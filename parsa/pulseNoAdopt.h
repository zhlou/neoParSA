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
    unirandom& rnd;
    //void *state_buf;
    //const int buf_size;
    unsigned counter;
    const unsigned frequency;
    dynDebug debugOut;
public:
    class Param {
    public:
        unsigned frequency;
        Param(xmlNode *root);
    };
    pulseNoAdopt(Problem &problem, const MPIState &mpiState,
                 unirandom& in_rnd, const Param &param);
    ~pulseNoAdopt();
    mixState Mix(aState &state);
    static const char * name;

};

#include "pulseNoAdopt.hpp"

#endif /* PULSENOADOPT_H_ */
