#ifndef TEMPFEEDBACK_H
#define TEMPFEEDBACK_H
#include <vector>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

#include "unirandom.h"
#include "dynDebug.h"
#include "aState.h"

// This is a templerature depedent move control

template <class Problem>
class tempFeedback
{
public:
    tempFeedback(Problem &in_problem, unirandom &in_rnd, const ptree &root);
    // now every array become a vector, I don't need to do anything in dtor
    virtual ~tempFeedback() {}
    double get_score() {return energy;}
    double propose(const aState &);
    void accept();
    void reject();
    static const char *name;
    void setDebug(debugStatus st, const char*outname=NULL)
    {debugOut.setDebug(st, outname);}
    void writeState(ptree &root) const;
    void readState(const ptree &root);
    double forceUpdateEnergy(){return energy = problem.get_score();}
protected:
    virtual void collectMoveStats(){}; // do nothing in base class
    void move_control(double);
    size_t nparams;
    std::vector<unsigned long> success;
    std::vector<unsigned long> moves;
    std::vector<double> theta_bars;
    std::vector<double> temp_coef;
    std::vector<double> last_s;
    std::vector<double> theta_mins;
    std::vector<double> theta_maxs;
    double energy;
    double prev_energy;
    size_t index;
    size_t step_cnt;
private:
    Problem &problem;
    unirandom &rnd;
    dynDebug debugOut;
    double prev_s;
    double proportion_gain; // should be const
    double integral_gain; // should be const
    unsigned long sweep;
    unsigned long interval; // should be const
    double target; // should be const
    double genThetaBar(size_t i, double s);
};

#include "tempFeedback.hpp"
#endif
