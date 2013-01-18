/*
 * simpleAnnealer.h
 *
 *  Created on: Dec 9, 2012
 *      Author: zhlou
 */

#ifndef SIMPLEANNEALER_H_
#define SIMPLEANNEALER_H_

#include "annealer.h"
#include <libxml/tree.h>

class simpleSchedule
{
public:
    simpleSchedule(movable *theproblem, xmlNode *root);
    virtual ~simpleSchedule();
    void initStats();
    double getInitS();
    int getInitLoop();
    void updateInitStep(bool accept, double energy);
    bool frozen();
    void resetSegmentStats();
    void updateStep(bool accept, double energy);
    double updateS(double s);
    bool inSegment();
    void updateSegment();
private:
    unsigned reject_cnt;
    double lambda;
    int init_loop;
    double init_S;
};

#endif /* SIMPLEANNEALER_H_ */
