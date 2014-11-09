/*
 * fixedTheta.hpp
 *
 *  Created on: Jul 26, 2014
 *      Author: zhlou
 */

#include "xmlUtils.h"
#include <exception>
#include <stdexcept>

template<class Problem>
const char *fixedThetaMove<Problem>::name = "fixedThetaMove";



template<class Problem>
fixedThetaMove<Problem>::fixedThetaMove(Problem & problem,
                                        unirandom& rnd,
                                        xmlNode* root):
        problem(problem), rnd(rnd), nparams(problem.getDimension()),
        index(0),counter(0)
{
    xmlNode *section = getSectionByName(root, "fixedThetaMove");
    if (section == NULL)
        throw std::runtime_error("Failed to find section fixedThetaMove");
    target = getPropDouble(section, "target");
    try {
        logInterval = getPropInt(section, "interval");
    } catch (std::exception &e) {
        logInterval = 100;
    }
    accumTheta = 0.0;
    accumAccept = 0;
    energy = problem.get_score();
    prev_energy = energy;
}

template<class Problem>
fixedThetaMove<Problem>::~fixedThetaMove()
{

}

template<class Problem>
double fixedThetaMove<Problem>::propose()
{
    prev_energy = energy;
    double uniform = 2.0 * rnd.random() - 1.0;
    double theta;
    if (uniform >= 0.)
        theta = -1 * target * log(abs(uniform));
    else
        theta = target * log(abs(uniform));

    problem.generateMove(index, theta);
    accumTheta += theta;
    energy = problem.get_score();

    return (energy - prev_energy);
}

template<class Problem>
void fixedThetaMove<Problem>::accept()
{
    accumAccept++;
    index ++;
    if (index == nparams)
        index = 0;
    counter ++;
    if (counter == logInterval) {
        counter = 0;
        procStats();
    }
}

template<class Problem>
void fixedThetaMove<Problem>::reject()
{
    problem.restoreMove(index);
    energy = prev_energy;
    index ++;
    if (index == nparams)
        index = 0;
    counter ++;
    if (counter == logInterval) {
        counter = 0;
        procStats();
    }
}

template<class Problem>
void fixedThetaMove<Problem>::procStats()
{
    double accRatio = (double)accumAccept / (double)logInterval;
    double avgTheta = accumTheta / logInterval;
    accumAccept = 0;
    accumTheta = 0.0;

    debugOut << accRatio << "\t" << avgTheta << std::endl;
}
