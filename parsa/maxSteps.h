/*
 * maxSteps.h
 *
 *  Created on: Apr 22, 2013
 *      Author: zhlou
 */

#ifndef MAXSTEPS_H_
#define MAXSTEPS_H_
#include "aState.h"
#include "dynDebug.h"
#include <libxml/tree.h>
class maxSteps {
private:
    int max_steps;
    mutable dynDebug debugOut;
public:
    class Param {
    public:
        int max_steps;
        debugStatus st;
        const char * debugname;
        Param(xmlNode *root, debugStatus st=ignore,
              const char *debugname = NULL);
    };
    maxSteps(const Param &param): max_steps(param.max_steps),
            debugOut(param.st, param.debugname){};
    ~maxSteps();
    void updateStep(bool, const aState &state){}
    bool frozen(const aState &state) const {return (state.step_cnt >= max_steps);}
};

#endif /* MAXSTEPS_H_ */
