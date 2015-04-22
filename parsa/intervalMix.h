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
    class Param{
    public:
        int interval;
        int reportNAdopt;
        Param(xmlNode *root);
    };
    //intervalMix(Problem &in_problem, const MPIState &mpiState,
    //            unirandom& in_rand, xmlNode *docroot);
    intervalMix(Problem &in_problem, const MPIState &mpiState,
                unirandom& in_rand, const Param &param);
    ~intervalMix();
    mixState Mix(aState &state);
    static const char * name;
    void setDebug(debugStatus st, const char* outname=NULL)
    {
        if (reportNAdopt) {
            mix.debugOut.setDebug(ignore, NULL);
            debugOut.setDebug(st, outname);
        } else {
            mix.debugOut.setDebug(st, outname);
            debugOut.setDebug(ignore, NULL);
        }
    }
private:
    Mixing<Problem> mix;
    //xmlNode *root;
    const int interval;
    int tau_count;
    const int reportNAdopt;
    dynDebug debugOut;
};

#include "intervalMix.hpp"


#endif /* INTERVALMIX_H_ */
