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
double annealer<Schedule, Move>::loop()
{
    if (!is_init)
        initMoves();
    bool accepted;
    while (!cooling.frozen()) {
        cooling.resetSegmentStats();
        do {
            accepted = step();
            cooling.updateStep(accepted, energy);
            s = cooling.updateS(s);
            step_cnt++;
        } while (cooling.inSegment());
        cooling.updateSegment();

    }
    std::cout << "Annealing stopped at s = " << s << std::endl
            << "Total steps is " << step_cnt << std::endl;
    return move.get_score();
}

template<class Schedule, class Move>
double annealer<Schedule, Move>::initMoves()
{
    bool accepted;
    for (int i = 0; i < cooling.getInitLoop(); i++) {
        accepted = step();
        cooling.updateInitStep(accepted, energy);
    }
    cooling.initStats();
    is_init = true;
    return move.get_score();

}

template<class Schedule, class Move>
annealer<Schedule, Move>::annealer(Schedule& in_cool, Move& in_move) :
cooling(in_cool), move(in_move)
{
    s = cooling.getInitS();
    is_init = false;
    step_cnt = 0;
    rnd_seed = time(NULL);
    energy = move.get_score();
}

template<class Schedule, class Move>
bool annealer<Schedule, Move>::step()
{
    double delta, crit, ran_n;
    delta = move.propose();
    crit = exp(-s * delta);
    ran_n = (double) rand_r(&rnd_seed) / RAND_MAX;
    if ((delta <= 0.0) || crit > ran_n) {
        move.accept();
        energy += delta;
        return true;
    } else {
        move.reject();
        return false;
    }
}
#endif /* ANNEALER_HPP_ */
