/*
 * criCount.cpp
 *
 *  Created on: Apr 5, 2013
 *      Author: zhlou
 */

#include "criCount.h"
#include <stdexcept>
#include <string>
#include <cmath>
#include <limits>

criCount::criCount(const criCount::Param &param) :
    freeze_crit(param.freeze_crit), cnt_crit(param.cnt_crit),
        interval(param.interval), step_cnt(0)
{
    old_energy = std::numeric_limits<double>::max();
    freeze_cnt = 0;
}

criCount::~criCount() {
    // TODO Auto-generated destructor stub
}

bool criCount::frozen(const aState& state)
{
    step_cnt ++;
    if (interval != step_cnt){
        return false;
    }
    step_cnt = 0;

    return checkFrozen(state);
}

bool criCount::checkFrozen(const aState& state)
{
    if (std::abs(old_energy - state.energy) < freeze_crit)
        freeze_cnt++;
    else
        freeze_cnt = 0;
    old_energy = state.energy;
    return (freeze_cnt >= cnt_crit);
}

criCount::Param::Param(const ptree &root)
{
    const ptree &section = root.get_child("count_criterion");
    freeze_crit = section.get<double>("<xmlattr>.freeze_crit");
    cnt_crit = section.get<int>("<xmlattr>.freeze_cnt");
    interval = section.get<int>("<xmlattr>.interval", 100);
}
