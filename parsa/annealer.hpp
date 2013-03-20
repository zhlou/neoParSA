/*
 * annealer.hpp
 *
 *  Created on: Jan 16, 2013
 *      Author: zhlou
 *  The implementation of the annealer template
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <stdexcept>
#include "utils.h"

using namespace std;

template<class Problem, class Schedule, template<class> class Move>
double annealer<Problem, Schedule, Move>::loop()
{
    clock_t start = clock();
    if (!is_init)
        initMoves();
    bool accepted;
    state.step_cnt = 0;
    do { // while not frozen
        cooling->resetSegmentStats();

        do { // while in segment
            accepted = step();
            cooling->updateStep(accepted, state);
            state.s = cooling->updateS(state);
            ++ (state.step_cnt);
        } while (cooling->inSegment(state));
        debugOut << state.step_cnt << " " << state.s << " "
                 << state.energy << endl;
        updateSegment(state);
    } while (!cooling->frozen(state));

    clock_t end = clock();
    tlaps = (double)(end - start)/CLOCKS_PER_SEC;

    std::cout << "Annealing stopped at s = " << state.s << std::endl
            << "Total steps is " << state.step_cnt << std::endl;
    return move->get_score();
}

template<class Problem, class Schedule, template<class> class Move>
double annealer<Problem, Schedule, Move>::initMoves()
{
    bool accepted;
    for (state.step_cnt = 0; state.step_cnt < initLoop; state.step_cnt++) {
        accepted = step();
        cooling->updateInitStep(accepted, state);
    }
    cooling->initStats(state);
    is_init = true;
    return move->get_score();

}

template<class Problem, class Schedule, template<class > class Move>
void annealer<Problem, Schedule, Move>::writeResult()
{
    xmlNode *result = xmlNewChild(xmlroot, NULL,
                                  (const xmlChar *)"annealing_result", NULL);
    std::ostringstream s_energy, s_step, s_time;
    s_energy << state.energy;
    s_step << state.step_cnt;
    s_time << tlaps;
    xmlNewProp(result, (const xmlChar*)"final_energy",
               (const xmlChar*)s_energy.str().c_str());
    xmlNewProp(result, (const xmlChar*)"max_count",
               (const xmlChar*)s_step.str().c_str());
    xmlNewProp(result, (const xmlChar*)"time",
                   (const xmlChar*)s_time.str().c_str());
}

template<class Problem, class Schedule, template<class > class Move>
void annealer<Problem, Schedule, Move>::initState(xmlNode* root)
{
    xmlNode* xmlsection = getSectionByName(root, "annealer_input");
    if (xmlsection == NULL) {
        throw runtime_error(
                string("Error: fail to find section annealer_input"));
    }
    initS = 1.0 / getPropDouble(xmlsection, "init_T");
    initLoop = getPropInt(xmlsection, "init_loop");
    state.s = initS;
    is_init = false;
    state.step_cnt = 0;
}

/*
 * When the constructor is called, all other parts should already be set up.
 */
template<class Problem, class Schedule, template<class> class Move>
annealer<Problem, Schedule, Move>::annealer(Problem &problem,
        unirandom * const in_rand, xmlNode *root) :
        rand(in_rand), xmlroot(root)
{
    initState(root);
    cooling = new Schedule(root);
    move = new Move<Problem>(problem, in_rand, root);
    state.energy = move->get_score();
}
/*
 * This construtor is used only by derived class that have there own way
 * of generating cooling schedule and move scheme
 */
template<class Problem, class Schedule, template<class > class Move>
annealer<Problem, Schedule, Move>::annealer(unirandom * const in_rand,
                                            xmlNode *root) :
        rand(in_rand), xmlroot(root)
{
    initState(root);
    cooling = NULL;
    move = NULL;
    tlaps = -1;
}

template<class Problem, class Schedule, template<class> class Move>
annealer<Problem, Schedule, Move>::~annealer()
{
    delete cooling;
    delete move;
}

template<class Problem, class Schedule, template<class> class Move>
bool annealer<Problem, Schedule, Move>::step()
{
    double delta, crit, ran_n;
    delta = move->propose();
    // cout << state.energy + delta << "@" << state.step_cnt << endl;
    crit = std::exp(-state.s * delta);
    ran_n = rand->random();
    if ((delta <= 0.0) || crit > ran_n) {
        move->accept();
        state.energy += delta;
        return true;
    } else {
        move->reject();
        return false;
    }
}
