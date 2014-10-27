/*
 * twoTierMix.h
 *
 *  Created on: Oct 23, 2014
 *      Author: zhlou
 */

#ifndef TWOTIERMIX_H_
#define TWOTIERMIX_H_

#include <libxml/tree.h>
#include "dynDebug.h"
#include "MPIState.h"
#include "mixState.h"
#include "Mixing.h"


template <class Problem>
class twoTierMix {
public:
    class Param{
    public:
        unsigned partSize;
        unsigned mixFreq;
        unsigned globalFreq;
        Param(xmlNode *docroot);
    };
    twoTierMix(Problem &in_problem, const MPIState &globalState,
               unirandom& in_rand, const Param &param);
    ~twoTierMix();
    mixState Mix(aState &state);
    static const char * name;
    void setDebug(debugStatus st, const char * outname=NULL){}
    // I don't know what to write for debug info
private:
    MPIState localState;
    const MPIState &globalState;
    Mixing<Problem> globalMix;
    Mixing<Problem> *localMix;
    unsigned mixFreq;
    unsigned globalFreq;
    unsigned tauCount;
    unsigned mixCount;
    int *local2Global;
};

#include "twoTierMix.hpp"
#endif /* TWOTIERMIX_H_ */
