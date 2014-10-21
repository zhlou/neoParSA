/*
 * aState.h
 *
 *  Created on: Jan 25, 2013
 *      Author: zhlou
 */

#ifndef ASTATE_H_
#define ASTATE_H_

struct aState
{
    unsigned long step_cnt;
    double s;
    double energy;
    double proposed; // proposed energy
};


#endif /* ASTATE_H_ */
