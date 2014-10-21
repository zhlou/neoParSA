/*
 * tempCount.h
 *
 *  Created on: May 5, 2013
 *      Author: zhlou
 */

#ifndef TEMPCOUNT_H_
#define TEMPCOUNT_H_

#include "dynDebug.h"
#include "aState.h"
#include <libxml/tree.h>

class tempCount {
private:
    mutable dynDebug debugOut;
    const int max_steps;
    int step_cnt;
    double target_s;
public:
    class Param {
    public:
        int max_steps;
        double target_s;
        debugStatus st;
        const char * debugname;
        Param(xmlNode *root, debugStatus st=ignore,
              const char *debugname = NULL);
    };
    tempCount(const Param &param) : max_steps(param.max_steps),
            step_cnt(0), target_s(param.target_s),
            debugOut(param.st, param.debugname) {}
    virtual ~tempCount();
    void updateStep(bool, const aState &state);
    bool frozen(const aState &state);
};

#endif /* TEMPCOUNT_H_ */
