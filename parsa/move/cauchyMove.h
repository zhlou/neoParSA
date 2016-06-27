#ifndef CAUCHYMOVE_H
#define CAUCHYMOVE_H

#include "unirandom.h"
#include "dynDebug.h"
#include "aState.h"

#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;


template <class Problem>
class cauchyMove
{
public:
    cauchyMove(Problem &in_problem, unirandom& in_rnd, const ptree &root);
    virtual ~cauchyMove();
    double get_score();
    double propose(const aState &);
    void accept();
    void reject();
    static const char *name;
    void setDebug(debugStatus st, const char*outname=NULL)
    {debugOut.setDebug(st, outname);}
    void writeState(ptree &root) const;
    void readState(const ptree &root);
    double forceUpdateEnergy() {return energy = problem.get_score();}
protected:
    int nparams;
    Problem &problem;
    dynDebug debugOut;
    virtual void collectMoveStats(){}; // do nothing in base class
    void move_control();
    long *success;
    long *moves;
    double *theta_bars;
    double *theta_mins;
    double *theta_maxs;
    double energy;
private:
    int index;
    unirandom& rnd;
    double prev_energy;
    double move_gain;
    long sweep;
    int move_interval;
    //static const double theta_min;
    double target;
};



#include "cauchyMove.hpp"

#endif
