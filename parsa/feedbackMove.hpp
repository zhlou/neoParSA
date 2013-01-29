/*
 * feedbackMove.hpp
 *
 *  Created on: Jan 21, 2013
 *      Author: zhlou
 */
#include <limits>
using namespace std;

template<class Problem, class Debug>
const double feedbackMove<Problem, Debug>::theta_min = 0.;

template<class Problem, class Debug>
feedbackMove<Problem, Debug>::feedbackMove(Problem& in_problem, xmlNode* root) :
        problem(in_problem)
{
    nparams = problem.getDimension();
    index = -1;
    success = new long[nparams];
    moves = new long[nparams];
    theta_bars = new double[nparams];
    for (int i = 0; i < nparams; ++i) {
        success[i] = 0;
        moves[i] = 0;
        theta_bars[i] = 1.0;
    }
    energy = problem.get_score();
    prev_energy = energy;
    sweep = 0;

    xmlNode *section = NULL;
    if (root != NULL) {
        section = root->children;
        while (section != NULL) {
            if (!xmlStrcmp(section->name, (xmlChar *)"move"))
                break;
            section = section->next;
        }
    }
    if (section == NULL) {
        move_gain = 0.03;
        move_interval = 100;
    } else {
        xmlChar *prop = NULL;
        if ((prop = xmlGetProp(section, (xmlChar *)"gain")) != NULL) {
            move_gain = strtod((char *)prop, NULL);
            xmlFree(prop);
            prop = NULL;
        } else
            move_gain = 0.03;
        if ((prop = xmlGetProp(section, (xmlChar *)"interval")) != NULL) {
            move_interval = atoi((char *)prop);
            xmlFree(prop);
            prop = NULL;
        } else
            move_interval = 100;

    }
}

template<class Problem, class Debug>
inline double feedbackMove<Problem, Debug>::get_score()
{
    return energy;
}

template<class Problem, class Debug>
double feedbackMove<Problem, Debug>::propose()
{
    index++;
    index %= nparams;
    if (index == 0) {
        ++sweep;
        if (sweep % move_interval == 0) {
            collectMoveStats();
            move_control();
        }

    }

    prev_energy = energy;
    problem.generateMove(index, theta_bars[index]);
    ++ moves[index]; // prefix increment has lower precedence than [], right?
    energy = problem.get_score();
    return (energy - prev_energy);
}

template<class Problem, class Debug>
void feedbackMove<Problem, Debug>::accept()
{
    ++ success[index];
}

template<class Problem, class Debug>
feedbackMove<Problem, Debug>::~feedbackMove()
{
    delete[] success;
    delete[] moves;
    delete[] theta_bars;
}

template<class Problem, class Debug>
void feedbackMove<Problem, Debug>::reject()
{
    problem.restoreMove(index);
    energy = prev_energy;
}

template<class Problem, class Debug>
void feedbackMove<Problem, Debug>::collectMoveStats()
{
    // Do nothing in base class.
}

template<class Problem, class Debug>
void feedbackMove<Problem, Debug>::move_control()
{
    for (int i = 0; i < nparams; ++i) {
        double acc_ratio = (double) success[i] / (double) moves[i];
        double x = log(theta_bars[i]);
        Debug::debugOut << i << "\t" << acc_ratio;
        x += move_gain * (acc_ratio - 0.44);
        theta_bars[i] = exp(x);
        if (theta_bars[i] < theta_min)
            theta_bars[i] = theta_min;
        Debug::debugOut << "\t" << theta_bars[i];
        success[i] = 0;
        moves[i] = 0;
    }
    Debug::debugOut << endl;
}
