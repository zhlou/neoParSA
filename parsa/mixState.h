/*
 * mixState.h
 *
 *  Created on: Mar 6, 2013
 *      Author: zhlou
 */

#ifndef MIXSTATE_H_
#define MIXSTATE_H_

#include <vector>
class mixState
{
public:
    std::vector<int> adoptList;
    mixState() {}
    mixState(int i) {adoptList.push_back(i);}
    bool doesMix() const {return !adoptList.empty();}

};


#endif /* MIXSTATE_H_ */
