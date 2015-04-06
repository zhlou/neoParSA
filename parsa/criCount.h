/*
 * criCount.h
 *
 *  Created on: Apr 5, 2013
 *      Author: zhlou
 */

#ifndef CRICOUNT_H_
#define CRICOUNT_H_

#include <libxml/tree.h>
#include "aState.h"

class criCountP;
class criCount {
public:
    class Param
    {
    public:
        double freeze_crit;
        int cnt_crit;
        int interval;
        Param(xmlNode *root);
    };
    criCount(const Param &param);
    ~criCount();
    void updateStep(bool, const aState &){}
    bool frozen(const aState &state);
    bool checkFrozen(const aState &state);
    friend class criCountP;
    friend class globalCount;
private:
    const double freeze_crit;
    double old_energy;
    const int cnt_crit;
    const int interval;
    int step_cnt;
    int freeze_cnt;
};

#endif /* CRICOUNT_H_ */
 