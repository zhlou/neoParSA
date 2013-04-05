/*
 * criCount.cpp
 *
 *  Created on: Apr 5, 2013
 *      Author: zhlou
 */

#include "criCount.h"
#include <limits>

criCount::criCount(const criCount::Param &param) :
    freeze_crit(param.freeze_crit), cnt_crit(param.cnt_crit)
{
    old_energy = std::numeric_limits<double>::max();
    freeze_cnt = 0;
}

criCount::~criCount() {
    // TODO Auto-generated destructor stub
}

bool criCount::frozen(const aState& state)
{
    if ((old_energy - state.energy) < freeze_crit)
        freeze_cnt++;
    else
        freeze_cnt = 0;
    old_energy = state.energy;
    return (freeze_cnt >= cnt_crit);
}
