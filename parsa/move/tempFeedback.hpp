#ifndef TEMPFEEDBACK_HPP
#define TEMPFEEDBACK_HPP

#include <limits>

template<class Problem>
const char *tempFeedback::name = "tempFeedback";

template<class Problem>
tempFeedback<Problem>::tempFeedback(Problem &in_problem, unirandom &in_rnd,
        const ptree &root) :
        nparams(in_problem.getDimension()),
        success(nparams, 0),
        moves(nparams, 0),
        theta_bars(nparams, 1.0),
        theta_mins(nparams, 0.),
        theta_maxs(nparams, std::numeric_limits<double>::max()),
        problem(in_problem),
        rnd(in_rnd), index(-1), sweep(0)
{
    energy = problem.get_score();
    prev_energy = energy;
    const ptree &sec_attr = root.get_child("tempMove.<xmlattr>");
    proportion_gain = sec_attr.get<double>("proportion_gain"); // TODO: find a good default value
    integral_gain = sec_attr.get<double>("integral_gain"); // TODO: find a good default value
    boost::optional<double> initTheta = sec_attr.get_optional<double>("init_theta_bar");
    if (initTheta) {
        for (int i = 0;  i < nparams; ++i) {
            theta_bars[i] = *initTheta;
        }
    }
    std::pair <tree::const_assoc_iterator, ptree::const_assoc_iterator> bonds
            = root.get_child("tempMove").equal_range("theta");
    for (ptree::const_assoc_iterator it = bounds.first; it != bounds.second; ++it) {
        const ptree &node = it->second;
        int i = node.get_value<int>(0);
        if (i < 0 || i >= nparams)
            continue;
        theta_bars[i] = node.get<double>("<xmalattr>.val", theta_bars[i]);
        theta_mins[i] = node.get<double>("<xmlattr>.min", theta_mins[i]);
        theta_maxs[i] = node.get<double>("<xmlattr>.max", theta_maxs[i]);
    }
}

#endif
