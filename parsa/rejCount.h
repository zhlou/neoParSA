/*
 * rejCount.h
 *
 *  Created on: Apr 3, 2013
 *      Author: zhlou
 */

#ifndef REJCOUNT_H_
#define REJCOUNT_H_
#include "aState.h"
class rejCount {
public:
    class Param {
    public:
        int max_rej;
        int output_freq;
    };
    rejCount(const Param &param);
    virtual ~rejCount();
    void updateState(bool accept, const aState &state);
    bool frozen(const aState &state);
protected:
    dynDebug debugOut;
    int max_rej;
    int output_freq;
    int reject_cnt;
};

#endif /* REJCOUNT_H_ */
