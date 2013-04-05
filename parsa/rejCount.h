/*
 * rejCount.h
 *
 *  Created on: Apr 3, 2013
 *      Author: zhlou
 */

#ifndef REJCOUNT_H_
#define REJCOUNT_H_

#include <libxml/tree.h>
#include "aState.h"
#include "dynDebug.h"

class rejCount {
public:
    class Param {
    public:
        int max_rej;
        //int output_freq;
        debugStatus st;
        const char * debugname;
    };
    rejCount(const Param &param);
    virtual ~rejCount();
    void updateState(bool accept, const aState &state);
    bool frozen(const aState &state) const;
protected:
    mutable dynDebug debugOut;
    int max_rej;
    //int output_freq;
    int reject_cnt;
};
rejCount::Param rejCountParamXML(xmlNode *root, debugStatus st=ignore,
                                 const char * debugname=NULL);
#endif /* REJCOUNT_H_ */
