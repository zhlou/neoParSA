/*
 * simpleAnnealer.h
 *
 *  Created on: Dec 9, 2012
 *      Author: zhlou
 */

#ifndef SIMPLEANNEALER_H_
#define SIMPLEANNEALER_H_

#include "annealer.h"
#include "aState.h"
#include <libxml/tree.h>

class simpleSchedule
{
public:
    simpleSchedule(xmlNode *root);
    virtual ~simpleSchedule();
    void initStats(aState state);
    //double getInitS();
    //int getInitLoop();
    void updateInitStep(bool accept, aState state);
    bool frozen(aState state);
    void resetSegmentStats(){};
    void updateStep(bool accept, aState state);
    double updateS(aState state);
    bool inSegment(aState state) {return false;}
    void updateSegment(aState state) {};
    bool needMix(){return false;}
private:
    unsigned reject_cnt;
    double lambda;
    //int init_loop;
    //double init_S;
};

#endif /* SIMPLEANNEALER_H_ */
