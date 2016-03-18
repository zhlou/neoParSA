// Created 2/13/2016
#ifndef SCHEDULETHETAMOVE_HPP
#define SCHEDULETHETAMOVE_HPP

#include <string>
#include "utils/vectorUtils.h"

template<class Problem>
const char *scheduleThetaMove<Problem>::name = "scheduleThetaMove";

template<class Problem>
scheduleThetaMove<Problem>::scheduleThetaMove(Problem &in_problem, unirandom &in_rnd, const ptree &root) :
        problem(in_problem),
        rnd(in_rnd),
        index(std::numeric_limits<unsigned>::max()), //UINT_MAX + 1 == 0
        currentTabRow(0),
        sweep(0)
{
    // initialize state
    nparams = problem.getDimension();
    energy = problem.get_score();
    prev_energy = energy;
    moves.assign(nparams, 0);
    success.assign(nparams, 0);

    // read config from root
    const ptree &sec_attr = root.get_child("scheduleThetaMove.<xmlattr>");
    std::string filename = sec_attr.get<std::string>("thetaFile");
    readArray(thetaTab, filename.c_str());
    interval = sec_attr.get<unsigned>("interval", 0);

}

template<class Problem>
double scheduleThetaMove<Problem>::propose(const aState &state)
{
    const double s = state.s;
    ++ index;
    if (index == nparams) {
        index = 0;
        ++ sweep;
        if (sweep == interval) {
            sweep = 0;
            // output prolix statistics
            debugOut << s;
            for (int i = 0; i < nparams; ++i) {
                double acc_ratio = (double) success[i] / (double) moves[i];
                debugOut << "\t" << acc_ratio;
                success[i] = 0;
                moves[i] = 0;
            }
            for (int i = 0; i < nparams; ++i) {
                debugOut << "\t" << getThetaByS(s, i);
            }
            debugOut << std::endl;
        }
    }
    prev_energy = energy;
    double theta_bar = getThetaByS(s, index);
    double theta = rnd.laplace(theta_bar);
    problem.generateMove(index, theta);
    ++ moves[index];
    energy = problem.get_score();
    return (energy - prev_energy);
}

template<class Problem>
double scheduleThetaMove<Problem>::getThetaByS(double s, unsigned index)
{
    while (currentTabRow < thetaTab.size()) {
        if (s < thetaTab[currentTabRow][0]) {
            break;
        }
        ++ currentTabRow;
    }
    if (0 == currentTabRow) {
        return thetaTab[0][index + 1];
    } else if (thetaTab.size() > currentTabRow) {
        // linear interpolation
        double sc = thetaTab[currentTabRow][0];
        double scp = thetaTab[currentTabRow - 1][0];
        double t = (s - scp)/(sc - scp);
        return (1-t) * thetaTab[currentTabRow -1][index + 1]
                + t * thetaTab[currentTabRow][index + 1];
    } else { // thetaTab.size() == currentTabRow
        return thetaTab[currentTabRow - 1][index + 1];
    }
}

template<class Problem>
void scheduleThetaMove<Problem>::accept()
{
    ++ success[index];
}

template<class Problem>
void scheduleThetaMove<Problem>::reject()
{
    problem.restoreMove(index);
    energy = prev_energy;
}

#endif
