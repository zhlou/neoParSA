/*
 * cauchyMove.hpp
 *
 *  Created on: Jun 22, 2016
 *      Author: zhlou
 */
#include <limits>
#include <cstdlib>
#include <cmath>

//template<class Problem>
//const double feedbackMove<Problem>::theta_min = 0.;

template<class Problem>
const char *cauchyMove<Problem>::name = "cauchyMove";


template<class Problem>
cauchyMove<Problem>::cauchyMove(Problem &in_problem, unirandom &in_rand,
                                    const ptree &root) :
        problem(in_problem), rnd(in_rand), index(-1), sweep(0)
{
    nparams = problem.getDimension();
    energy = problem.get_score();
    prev_energy = energy;
    const ptree &sec_attr = root.get_child("cauchyMove.<xmlattr>");
    move_gain = sec_attr.get<double>("gain", 0.03);
    move_interval = sec_attr.get<int>("interval", 100);
    target = sec_attr.get<double>("target", 0.44);
    double initTheta = sec_attr.get<double>("init_theta_bar", 1.0);
    success = new long[nparams];
    moves = new long[nparams];
    theta_bars = new double[nparams];
    theta_mins = new double[nparams];
    theta_maxs = new double[nparams];
    for (int i = 0; i < nparams; ++i) {
        success[i] = 0;
        moves[i] = 0;
        theta_bars[i] = initTheta;
        theta_mins[i] = 0.;
        theta_maxs[i] = numeric_limits<double>::max();
    }
    // now initialize theta_mins and theta_maxs if there're any in the input file
    std::pair <ptree::const_assoc_iterator, ptree::const_assoc_iterator> bounds
            = root.get_child("cauchyMove").equal_range("theta");
    for (ptree::const_assoc_iterator it = bounds.first;
            it != bounds.second; ++it) {
        const ptree &node = it->second;
        int i = node.get_value<int>(0);
        if (i < 0 || i >= nparams)
            continue;
        theta_mins[i] = node.get<double>("<xmlattr>.min", theta_mins[i]);
        theta_maxs[i] = node.get<double>("<xmlattr>.max", theta_maxs[i]);
    }
}


template<class Problem>
inline double cauchyMove<Problem>::get_score()
{
    return energy;
}

template<class Problem>
double cauchyMove<Problem>::propose(const aState &)
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
    // generate theta from theta bar here and pass only theta
    // to the problem so problem doesn't have to have random
    // number generator
    double theta = rnd.cauchy(theta_bars[index]);
    problem.generateMove(index, theta);
    ++ moves[index]; // prefix increment has lower precedence than [], right?
    energy = problem.get_score();
    return (energy - prev_energy);
}

template<class Problem>
void cauchyMove<Problem>::accept()
{
    ++ success[index];
}

template<class Problem>
cauchyMove<Problem>::~cauchyMove()
{
    delete[] success;
    delete[] moves;
    delete[] theta_bars;
    delete[] theta_mins;
    delete[] theta_maxs;
}

template<class Problem>
void cauchyMove<Problem>::reject()
{
    problem.restoreMove(index);
    energy = prev_energy;
}

template<class Problem>
void cauchyMove<Problem>::move_control()
{
    debugOut << sweep;
    for (int i = 0; i < nparams; ++i) {
        double acc_ratio = (double) success[i] / (double) moves[i];
        double x = log(theta_bars[i]);
        debugOut << "\t" << acc_ratio;
        x += move_gain * (acc_ratio - target);
        theta_bars[i] = exp(x);
        if (theta_bars[i] < theta_mins[i])
            theta_bars[i] = theta_mins[i];
        if (theta_bars[i] > theta_maxs[i])
            theta_bars[i] = theta_maxs[i];
        success[i] = 0;
        moves[i] = 0;
    }
    if (!debugOut.isIgnore()){
        for (int i = 0; i < nparams; ++i) {
            debugOut << "\t" << theta_bars[i];
        }
    }

    debugOut << endl;
}

template<class Problem>
void cauchyMove<Problem>::writeState(ptree &root) const
{
    ptree node;
    for (int i = 0; i < nparams; ++i) {
        ptree &param = node.add("param", i);
        param.put("<xmlattr>.theta", theta_bars[i]);
    }
    root.put_child("moveSize", node);
}


template<class Problem>
void cauchyMove<Problem>::readState(const ptree &root)
{
    boost::optional<const ptree &> section = root.get_child_optional("moveSize");
    if (! section) {
        cerr << "Warning: fail to read moveSize from state file" << endl;
        return;
    }
    std::pair <ptree::const_assoc_iterator, ptree::const_assoc_iterator> bounds
            = (*section).equal_range("param");
    for (ptree::const_assoc_iterator it = bounds.first;
            it != bounds.second; ++it) {
        const ptree &node = it->second;
        int i = node.get_value<int>(0);
        if (i < 0 || i >= nparams)
            continue;
        theta_bars[i] = node.get<double>("<xmlattr>.theta", theta_bars[i]);
    }
}
