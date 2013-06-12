/*
 * pulseBcast.h
 *
 *  Created on: May 15, 2013
 *      Author: zhlou
 */

#ifndef PULSEBCAST_H_
#define PULSEBCAST_H_

#include "MPIState.h"
#include "unirandom.h"
#include "mixState.h"
#include <libxml/tree.h>

template <class Problem>
class pulseBcast
{
private:
    Problem &problem;
    const MPIState &mpi;
    unirandom * const rnd;
    void *state_buf;
    const int buf_size;
    unsigned counter;
    unsigned frequency;
    dynDebug debugOut;
public:
    pulseBcast(Problem &problem, const MPIState &mpiState,
               unirandom * const in_rnd, xmlNode *docroot);
    ~pulseBcast();
    mixState Mix(aState &state);
    static const char * name;

};

#include "pulseBcast.hpp"

#endif /* PULSEBCAST_H_ */
