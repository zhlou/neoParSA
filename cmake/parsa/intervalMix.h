/*
 * intervalMix.h
 *
 *  Created on: Mar 21, 2013
 *      Author: zhlou
 */

#ifndef INTERVALMIX_H_
#define INTERVALMIX_H_

#include <libxml/tree.h>
#include "dynDebug.h"
#include "MPIState.h"
#include "mixState.h"
#include "Mixing.h"

template <class Problem>
class intervalMix {
public:
    intervalMix(Problem &in_problem, const MPIState &mpiState,
                unirandom * const in_rand, xmlNode *docroot);
    ~intervalMix();
    mixState Mix(aState &state);
    static const char * name;
    void setDebug(debugStatus st, const char* outname=NULL)
    {
        mix.debugOut.setDebug(st, outname);
    }
private:
    Mixing<Problem> mix;
    xmlNode *root;
    const int interval;
    int tau_count;
};

#include "intervalMix.hpp"


#endif /* INTERVALMIX_H_ */
