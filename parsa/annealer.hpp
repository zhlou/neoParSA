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
#include "xmlUtils.h"

using namespace std;
/*
 * When the constructor is called, all other parts should already be set up.
 */
template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
annealer<Problem, Schedule, FrozenCnd, Move>::annealer(Problem &problem,
        unirandom& in_rand, typename Schedule::Param scheParam,
        typename FrozenCnd::Param frozenParam,
        xmlNode *root) :
        rand(in_rand), xmlroot(root)
{
    initState(root);
    cooling = new Schedule(scheParam);
    move = new Move<Problem>(problem, in_rand, root);
    frozen = new FrozenCnd(frozenParam);
    state.energy = move->get_score();
    tlaps = -1;
}

/*
 * This construtor is used only by derived class that have there own way
 * of generating cooling schedule and move scheme
 */
template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
annealer<Problem, Schedule, FrozenCnd, Move>::annealer(unirandom& in_rand,
        xmlNode *root) :
        rand(in_rand), xmlroot(root)
{
    initState(root);
    cooling = NULL;
    move = NULL;
    frozen = NULL;
    tlaps = -1;
}

template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
annealer<Problem, Schedule, FrozenCnd, Move>::~annealer()
{
    delete cooling;
    delete move;
    delete frozen;
}

template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
double annealer<Problem, Schedule, FrozenCnd, Move>::loop()
{
    clock_t start = clock();
    if (!is_init)
        initMoves();
    bool accepted;
    // state.step_cnt = 0;
    do { // while not frozen
        cooling->resetSegmentStats();

        do { // while in segment
            accepted = step();

            cooling->updateStep(accepted, state);
            frozen->updateStep(accepted, state);
            state.s = cooling->updateS(state);
            ++ (state.step_cnt);
        } while (cooling->inSegment(state));

        updateSegment(state);
    } while (!frozen->frozen(state));

    clock_t end = clock();
    tlaps = (double)(end - start)/CLOCKS_PER_SEC;

    std::cout << "Annealing stopped at s = " << state.s << std::endl
            << "Total steps is " << state.step_cnt << std::endl;
    return move->get_score();
}

template<class Problem, class Schedule, class FrozenCnd, template<class > class Move>
double annealer<Problem, Schedule, FrozenCnd, Move>::fixedTMoves(double S, long steps)
{
    double prev_s = state.s;
    state.s = S;
    for (int i = 0; i < steps; ++i) {
        step();
    }
    state.s = prev_s;
    return move->get_score();
}

template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
double annealer<Problem, Schedule, FrozenCnd, Move>::initMoves()
{
    fixedTMoves(state.s, initLoop);
    bool accepted;
    for (state.step_cnt = 0; state.step_cnt < initLoop; ++ state.step_cnt) {
        accepted = step();
        cooling->updateInitStep(accepted, state);
    }
    cooling->initStats(state);
    is_init = true;
    return move->get_score();

}

template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
double annealer<Problem, Schedule, FrozenCnd, Move>::initMovesOnly()
{
    fixedTMoves(state.s, initLoop);
    return move->get_score();

}

template<class Problem, class Schedule, class FrozenCnd, template<class > class Move>
void annealer<Problem, Schedule, FrozenCnd, Move>::writeResultData(xmlNode* result) {
    std::ostringstream s_energy, s_step, s_time;
    s_energy << state.energy;
    s_step << state.step_cnt;
    s_time << tlaps;
    xmlNewProp(result, (const xmlChar*) ("final_energy"),
               (const xmlChar*) (s_energy.str().c_str()));
    xmlNewProp(result, (const xmlChar*) ("max_count"),
               (const xmlChar*) (s_step.str().c_str()));
    xmlNewProp(result, (const xmlChar*) ("time"),
               (const xmlChar*) (s_time.str().c_str()));
}

template<class Problem, class Schedule, class FrozenCnd, template<class > class Move>
void annealer<Problem, Schedule, FrozenCnd, Move>::writeMethodText(xmlNode *method)
{
    xmlNewProp(method, (const xmlChar*)"cooling-schedule",
               (const xmlChar*)Schedule::name);
    xmlNewProp(method, (const xmlChar*)"move-generation",
               (const xmlChar*)Move<Problem>::name);
}
template<class Problem, class Schedule, class FrozenCnd, template<class > class Move>
void annealer<Problem, Schedule, FrozenCnd, Move>::writeResult()
{
    xmlNode *result = getSectionByName(xmlroot, "annealing_result");
    if (result != NULL) {
        xmlUnlinkNode(result);
        xmlFreeNode(result);
    }
    result = xmlNewChild(xmlroot, NULL, (const xmlChar *)"annealing_result",
                         NULL);
    writeResultData(result);
    xmlNode *method = getSectionByName(xmlroot, "method");
    if (method != NULL) {
        xmlUnlinkNode(method);
        xmlFreeNode(method);
    }
    method = xmlNewChild(xmlroot, NULL, (const xmlChar *)"method", NULL);
    writeMethodText(method);
}

template<class Problem, class Schedule, class FrozenCnd, template<class > class Move>
void annealer<Problem, Schedule, FrozenCnd, Move>::initState(xmlNode* root)
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


template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
bool annealer<Problem, Schedule, FrozenCnd, Move>::step()
{
    double delta, crit, ran_n;
    bool flag;
    delta = move->propose();
    state.proposed = state.energy + delta;
    // cout << state.energy + delta << "@" << state.step_cnt << endl;
    crit = std::exp(-state.s * delta);
    ran_n = rand.random();
    if ((delta <= 0.0) || crit > ran_n) {
        move->accept();
        state.energy += delta;

        flag = true;
    } else {
        move->reject();
        flag = false;
    }
    debugOut << state.step_cnt << " " << state.s << " "
             << state.energy << " " << state.proposed << endl;
    return flag;
}
