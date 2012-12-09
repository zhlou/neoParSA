/*
 * simpleAnnealer.h
 *
 *  Created on: Dec 9, 2012
 *      Author: zhlou
 */

#ifndef SIMPLEANNEALER_H_
#define SIMPLEANNEALER_H_

#include "annealer.h"

class simpleAnnealer: public annealer
{
public:
    simpleAnnealer(movable *theproblem, xmlNode *root);
    virtual ~simpleAnnealer();
protected:
    unsigned reject_cnt;
    void updateStats(bool accept, double delta);
    bool frozen();
    void cool_s();
};

#endif /* SIMPLEANNEALER_H_ */
