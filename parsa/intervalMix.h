/*
 * intervalMix.h
 *
 *  Created on: Mar 21, 2013
 *      Author: zhlou
 */

#ifndef INTERVALMIX_H_
#define INTERVALMIX_H_

#include "MPIState.h"
#include "mixState.h"
#include <libxml/tree.h>
template <class Problem>
class intervalMix {
public:
    intervalMix(Problem &in_problem, const MPIState &mpiState,
                unirandom * const in_rand, xmlNode *docroot);
    ~intervalMix();
    mixState Mix(aState &state);
    static const char * name;
private:
    Mixing<Problem> mix;
    xmlNode *root;
    const int interval;
    int tau_count;
};

#include "intervalMix.hpp"


#endif /* INTERVALMIX_H_ */
