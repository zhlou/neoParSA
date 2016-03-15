/*
 * rejCount.h
 *
 *  Created on: Apr 3, 2013
 *      Author: zhlou
 */

#ifndef REJCOUNT_H_
#define REJCOUNT_H_

#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

#include "aState.h"
#include "dynDebug.h"
class rejCountP;
class rejCount {
public:
    class Param {
    public:
        int max_rej;
        //int output_freq;
        debugStatus st;
        const char * debugname;
        Param(const ptree &root, debugStatus st=ignore, const char *debugname=NULL);
    };
    rejCount(const Param &param);
    ~rejCount(){};
    void updateStep(bool accept, const aState &state);
    bool frozen(const aState &state) const;
    friend class rejCountP;
private:
    mutable dynDebug debugOut;
    const int max_rej;
    //int output_freq;
    int reject_cnt;
};

#endif /* REJCOUNT_H_ */
