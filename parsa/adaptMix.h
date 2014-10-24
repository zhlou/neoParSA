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
#include "mixState.h"
#include "Mixing.h"
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
    class Param{
    public:
        double adaptCoef;
        Param(xmlNode *root);
    };
    adaptMix(Problem &in_problem, const MPIState &mpiState,
             unirandom& in_rnd, const Param &param);
    ~adaptMix();
    mixState Mix(aState &state);
    static const char *name;
    void setDebug(debugStatus st, const char* outname=NULL)
    {
        mix.debugOut.setDebug(st, outname);
    }
private:
    Mixing<Problem> mix;
    unirandom& rnd;
    //xmlNode *root;
    int nnodes;


    const double adaptCoef;
};



#include "adaptMix.hpp"

#endif /* ADAPTMIX_H_ */
