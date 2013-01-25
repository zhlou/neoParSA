/*
 * annealer.hpp
 *
 *  Created on: Jan 16, 2013
 *      Author: zhlou
 *  The implementation of the annealer template
 */

#include <iostream>
#include <cmath>
#include <stdexcept>

using namespace std;

template<class Schedule, class Move, class Random>
double annealer<Schedule, Move, Random>::loop()
{
    if (!is_init)
        initMoves();
    bool accepted;
    do { // while not frozen
        cooling.resetSegmentStats();
        do { // while in segment
            accepted = step();
            cooling.updateStep(accepted, energy);
            s = cooling.updateS(s);
            step_cnt++;
        } while (cooling.inSegment());
        cooling.updateSegment();
    } while (!cooling.frozen());

    std::cout << "Annealing stopped at s = " << s << std::endl
            << "Total steps is " << step_cnt << std::endl;
    return move.get_score();
}

template<class Schedule, class Move, class Random>
double annealer<Schedule, Move, Random>::initMoves()
{
    bool accepted;
    for (int i = 0; i < initLoop; i++) {
        accepted = step();
        cooling.updateInitStep(accepted, energy);
    }
    cooling.initStats();
    is_init = true;
    return move.get_score();

}

/*
 * When the constructor is called, all other parts should already be set up.
 */
template<class Schedule, class Move, class Random>
annealer<Schedule, Move, Random>::annealer(Schedule& in_cool, Move& in_move,
        Random &in_rand, xmlNode *root) :
        cooling(in_cool), move(in_move), rand(in_rand), xmlroot(root)
{
    xmlNode *xmlsection = getSectionByName(root, (char *)"annealer_input");

    if (xmlsection == NULL) {
        throw runtime_error(string("Error: fail to find section annealer_input"));
    }
    initS = 1.0 / getPropDouble(xmlsection, (char *)"init_T");
    initLoop = getPropInt(xmlsection, (char *)"init_loop");
    s = initS;
    is_init = false;
    step_cnt = 0;
    energy = move.get_score();
}

template<class Schedule, class Move, class Random>
annealer<Schedule, Move, Random>::~annealer()
{
}

template<class Schedule, class Move, class Random>
bool annealer<Schedule, Move, Random>::step()
{
    double delta, crit, ran_n;
    delta = move.propose();
    crit = std::exp(-s * delta);
    ran_n = rand.random();
    if ((delta <= 0.0) || crit > ran_n) {
        move.accept();
        energy += delta;
        return true;
    } else {
        move.reject();
        return false;
    }
}
