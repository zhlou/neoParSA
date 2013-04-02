/*
 * pSimpleSchedule.h
 *
 *  Created on: Apr 2, 2013
 *      Author: zhlou
 */

#ifndef PSIMPLESCHEDULE_H_
#define PSIMPLESCHEDULE_H_

#include "simpleScheduler.h"

class pSimpleSchedule : public simpleSchedule
{
public:
    pSimpleSchedule(xmlNode *root, const MPIState &) : simpleSchedule(root) {}
    static const char *name;
};
const char *pSimpleSchedule::name = "parallel exponential";
#endif /* PSIMPLESCHEDULE_H_ */
