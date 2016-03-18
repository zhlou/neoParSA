/*
 * fixedTheta.hpp
 *
 *  Created on: Jul 26, 2014
 *      Author: zhlou
 */

#include <exception>
#include <stdexcept>

template<class Problem>
const char *fixedThetaMove<Problem>::name = "fixedThetaMove";




template<class Problem>
fixedThetaMove<Problem>::fixedThetaMove(Problem &in_problem, unirandom &in_rnd,
                                        const ptree &root):
        problem(in_problem), rnd(in_rnd), index(0), counter(0),
        accumTheta(0.0), accumAccept(0), nparams(in_problem.getDimension())
{
    energy = problem.get_score();
    prev_energy = energy;

    const ptree &sec_attr = root.get_child("fixedThetaMove.<xmlattr>");
    target = sec_attr.get<double>("target");
    logInterval = sec_attr.get<int>("interval", 100);

}

template<class Problem>
fixedThetaMove<Problem>::~fixedThetaMove()
{

}

template<class Problem>
double fixedThetaMove<Problem>::propose(const aState &)
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


template<class Problem>
void fixedThetaMove<Problem>::writeState(ptree &root) const
{
    ptree node;
    /* for now fixedThetaMove only has one theta
    int i;
    for (i = 0; i < nparams; ++i) {
        asprintf(&paramNumString, "%d", i);
        paramIter = xmlNewChild(moveNode, NULL, BAD_CAST "param",
                BAD_CAST paramNumString);
        free(paramNumString);
        asprintf(&thetaString, "%.15g", theta_bars[i]);
        xmlNewProp(paramIter, BAD_CAST "theta", BAD_CAST thetaString);
        free(thetaString);

    }
    */
    ptree &param = node.put("param", 0);
    param.put("<xmlattr>.theta", target);
    root.put_child("moveSize", node);
}

template<class Problem>
void fixedThetaMove<Problem>::readState(const ptree &)
{
    // why oh why yet another function that does nothing
}
