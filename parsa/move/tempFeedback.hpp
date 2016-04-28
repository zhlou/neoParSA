#ifndef TEMPFEEDBACK_HPP
#define TEMPFEEDBACK_HPP

#include <limits>
#include <cmath>

template<class Problem>
const char *tempFeedback<Problem>::name = "tempFeedback";

template<class Problem>
tempFeedback<Problem>::tempFeedback(Problem &in_problem, unirandom &in_rnd,
        const ptree &root) :
        nparams(in_problem.getDimension()),
        success(nparams, 0),
        moves(nparams, 0),
        theta_bars(nparams, 1.0),
        temp_coef(nparams, 0.0),
        last_s(nparams, 0.0),
        theta_mins(nparams, 0.),
        theta_maxs(nparams, std::numeric_limits<double>::max()),
        problem(in_problem),
        rnd(in_rnd), index(0), step_cnt(0), sweep(0)
{
    energy = problem.get_score();
    prev_energy = energy;
    const ptree &sec_attr = root.get_child("tempMove.<xmlattr>");
    proportion_gain = sec_attr.get<double>("proportion_gain"); // TODO: find a good default value
    integral_gain = sec_attr.get<double>("integral_gain"); // TODO: find a good default value
    target = sec_attr.get<double>("target", 0.44);
    interval = sec_attr.get<unsigned long>("interval", 100);
    boost::optional<double> initTheta = sec_attr.get_optional<double>("init_theta_bar");
    if (initTheta) {
        for (int i = 0;  i < nparams; ++i) {
            theta_bars[i] = *initTheta;
        }
    }
    std::pair <ptree::const_assoc_iterator, ptree::const_assoc_iterator> bounds
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

template <class Problem>
double tempFeedback<Problem>::propose(const aState &state)
{
    if (index == nparams) {
        index = 0;
    }
    if (step_cnt == nparams) {
        step_cnt = 0;
        ++sweep;
        if (sweep == interval) {
            sweep = 0;
            collectMoveStats();
            move_control(state.s);
        }
    }
    prev_energy = energy;
    double theta = rnd.laplace(genThetaBar(index, state.s));
    problem.generateMove(index, theta);
    ++ moves[index];
    energy = problem.get_score();
    return (energy - prev_energy);
}

template<class Problem>
void tempFeedback<Problem>::accept()
{
    ++ success[index];
    index ++;
    step_cnt ++;
}

template<class Problem>
void tempFeedback<Problem>::reject()
{
    problem.restoreMove(index);
    energy = prev_energy;
    index ++;
    step_cnt ++;
}

template<class Problem>
void tempFeedback<Problem>::move_control(double s)
{
    debugOut << s;
    double lnDs = 0.0;
    if (0 != last_control_s) {
        lnDs = std::log(s) - std::log(last_control_s);
    }
    for (int i = 0; i < nparams; ++i) {
        double acc_ratio = (double)success[i] / (double) moves[i];
        debugOut << "\t" << acc_ratio;
        double delta_rho = acc_ratio - target;
        temp_coef[i] += integral_gain * delta_rho * lnDs;
        double lt = std::log(theta_bars[i]);
        lt += proportion_gain * delta_rho;
        double t = std::exp(lt);
        if (t < theta_mins[i]) {
            t = theta_mins[i];
        }
        if (t > theta_maxs[i]) {
            t = theta_maxs[i];
        }
        theta_bars[i] = t;
        success[i] = 0;
        moves[i] = 0;
    }
    last_control_s = s;
    if (!debugOut.isIgnore()){
        for (int i = 0; i < nparams; ++i) {
            debugOut << "\t" << theta_bars[i];
        }
    }

    debugOut << endl;
}

template<class Problem>
double tempFeedback<Problem>::genThetaBar(size_t i, double s)
{
    if (0 != last_s[i]) {
        double lt = std::log(theta_bars[i]);
        lt += temp_coef[i] * (std::log(s) - std::log(last_s[i]));
        double t = std::exp(lt);
        if (t < theta_mins[i]) {
            t = theta_mins[i];
        }
        if (t > theta_maxs[i]) {
            t = theta_maxs[i];
        }
        theta_bars[i] = t;
    }
    last_s[i] = s;
    return theta_bars[i];
}

template<class Problem>
void tempFeedback<Problem>::writeState(ptree &root) const
{
    ptree node;
    for (int i = 0; i < nparams; ++i) {
        ptree &param = node.add("param", i);
        param.put("<xmlattr>.theta", theta_bars[i]);
        param.put("<xmlattr>.temp_coef", temp_coef[i]);
    }
    root.put_child("moveSize", node);
}
template<class Problem>
void tempFeedback<Problem>::readState(const ptree &root)
{
    boost::optional<const ptree &> section = root.get_child_optional("moveSize");
    if (! section) {
        std::cerr << "Warning: fail to read moveSize from state file\n";
        return;
    }
    std::pair <ptree::const_assoc_iterator, ptree::const_assoc_iterator>
            bounds = section->equal_range("param");
    for (ptree::const_assoc_iterator it = bounds.first;
            it != bounds.second; ++it) {
        const ptree &node = it->second;
        int i = node.get_value<int>(0);
        if (i < 0 || i >= nparams)
            continue;
        theta_bars[i] = node.get<double>("<xmlattr>.theta", theta_bars[i]);
        temp_coef[i] = node.get<double>("<xmlattr>.temp_coef", temp_coef[i]);
    }
}

#endif
