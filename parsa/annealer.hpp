/*
 * annealer.hpp
 *
 *  Created on: Jan 16, 2013
 *      Author: zhlou
 */

#ifndef ANNEALER_HPP_
#define ANNEALER_HPP_

#include <iostream>

template<class Schedule, class Move>
inline double annealer<Schedule, Move>::loop()
{
    if (!is_init)
        initMoves();
    bool accepted;
    while (!frozen()) {
        cooling.resetSegmentStats();
        do {
            accepted = move(); // TODO move may be expanded here
            cooling.updateStep(accepted, energy);
            s = cooling.updateS(s);
            step_cnt++;
        } while (cooling.inSegment());
        cooling.updateSegment();

    }
    std::cout << "Annealing stopped at s = " << s << std::endl << "Total steps is "
            << step_cnt << std::endl;
    return problem->get_score();
}

template<class Schedule, class Move>
inline double annealer<Schedule, Move>::initMoves()
{

}


#endif /* ANNEALER_HPP_ */
