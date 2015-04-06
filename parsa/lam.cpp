/*
 * lam.cpp
 *
 *  Created on: Dec 8, 2012
 *      Author: zhlou
 */

#include "lam.h"
#include "xmlUtils.h"
#include <string>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <iostream>
using namespace std;

const double lam::UNINITIALIZED = numeric_limits<double>::max();
const char * lam::name = "Lam";
// in newer compilers, string constants cast to char * will generate
// warning messages unless done explicitly.
lam::Param::Param(xmlNode* root, debugStatus in_st, const char* name) :
        st(in_st), outname(name)
{
    xmlNode *section = getSectionByName(root, "annealer_input");

    if (section == NULL) {
        throw runtime_error(string("Error: fail to find section annealer_input"));
    }
    lambda = getPropDouble(section, "lambda");
    //init_S = 1.0 / getPropDouble(section, "init_T");
    //init_loop = getPropInt(section, "init_loop");

    section = getSectionByName(root, "lam");
    if (section == NULL)
        throw runtime_error(string("Error: fail to find section lam"));
    proc_tau = getPropInt(section, "tau");

    double memlength_mean = getPropDouble(section, "memLength_mean");
    double memlength_sd = getPropDouble(section, "memLength_sd");
    w_mean = 1.0 - proc_tau / (memlength_mean / lambda);
    if (w_mean < 0.)
        w_mean = 0.;
    w_sd = 1.0 - proc_tau / (memlength_sd / lambda);
    if (w_sd < 0.)
        w_sd = 0.;
    //freeze_crit = getPropDouble(section, "criterion");
    //cnt_crit = getPropInt(section, "freeze_cnt");
}

lam::lam(Param param) : proc_tau(param.proc_tau),
        //freeze_crit(param.freeze_crit),
        //cnt_crit(param.cnt_crit),
        debugOut(param.st, param.outname),
        lambda(param.lambda), w_mean(param.w_mean), w_sd(param.w_sd),
        acc_ratio(0.), success(0), vari(0.), mean(0.), fit_mean(NULL),
        fit_sd(NULL), old_energy(UNINITIALIZED), alpha(UNINITIALIZED),
        freeze_cnt(0), tau_count(0)
{
    debugOut << "iterations S\t dS/S\t meanE\tsdE\t(e)meanE\t(e)sdE\tacc\talpha"
             << endl;
}

lam::~lam()
{
    if (fit_mean != NULL)
        delete fit_mean;
    if (fit_sd != NULL)
        delete fit_sd;
}
/*
void lam::local_frozen(const aState& state)
{
    if (abs(state.energy - old_energy) < freeze_crit)
        freeze_cnt++;
    else
        freeze_cnt = 0;

    old_energy = state.energy;

}


 double lam::loop()
 {
 long unsigned step_cnt = 0;
 int accept = 0, i;
 double delta, vari;
 while (!frozen()) {
 accept = 0;
 vari = 0.;
 for (i = 0; i < proc_tau; i++) {
 if ((delta = move()) != 0.)
 accept ++;

 }

 }

 return problem->get_score();
 }


bool lam::frozen(aState state)
{
    local_frozen(state);
    return global_frozen();

}
 */
void lam::updateStep(bool accept, aState state)
{
    //++ step_cnt;
    if (accept)
        ++ success;
    double estimate_mean = fit_mean->getEstimate(state.s);
    double d = state.energy - estimate_mean;
    mean += state.energy;
    vari += d * d;
}

double lam::updateS(aState state)
{
    double estimate_sd = fit_sd->getEstimate(state.s);
    double d = state.s * estimate_sd;
    return state.s + lambda * alpha / (d * d * estimate_sd);
}

void lam::resetSegmentStats()
{
    acc_ratio = 0.;
    success = 0;
    vari = 0.;
    mean = 0.;
    //old_energy = energy;
}

void lam::updateEstimators(double s)
{
    fit_mean->fullUpdate(1.0 / mean, s);
    fit_sd->fullUpdate(1.0 / sqrt(vari), s);
}

void lam::updateSegment(aState state)
{
    collectStats();
    updateEstimators(state.s);
    updateLam();
    double estimate_mean = fit_mean->getEstimate(state.s);
    double estimate_sd = fit_sd->getEstimate(state.s);
    double d = state.s *estimate_sd;
    debugOut << state.step_cnt << " " << state.s << " "
             << lambda * alpha / (d * d * estimate_sd) << " "
             << mean << " " << sqrt(vari) << " "
             << estimate_mean << " " << estimate_sd << " "
             << acc_ratio << " " << alpha << endl;
//    resetSegmentStats();
}

void lam::resetLam()
{
    fit_mean->reset();
    fit_sd->reset();
}

void lam::initStats(aState state)
{
    collectInitStats(state.step_cnt);

    initStatsCore(state);
    tau_count = 0;


    //cout << state.step_cnt << ": " << mean << "\t" << sd << "\t" << alpha
    //        <<endl;
    //resetSegmentStats(); // we don't need this since it is called automatically
                           // before every segment
}

void lam::initStats(double initMean, double initVar, double initAccRatio, aState state)
{

    collectInitStats(initMean, initVar, initAccRatio);
    initStatsCore(state);
    tau_count = 0;
}

void lam::initStatsCore(const aState& state) 
{
    double sd = sqrt(vari);
    fit_mean = new invLinearFit(w_mean, mean, state.s, vari / (mean * mean));
    fit_sd = new invLinearFit(w_sd, sd, state.s, sd / mean);

    double d = (1.0 - acc_ratio) / (2.0 - acc_ratio);
    alpha = 4.0 * acc_ratio * d * d;
    old_energy = state.energy;
}

void lam::updateInitStep(bool accept, aState state)
{
    if (accept)
        success++;
    mean += state.energy;
    vari += state.energy * state.energy;
}

void lam::collectStats()
{
    mean /= proc_tau;
    vari /= proc_tau;
    acc_ratio = (double) success / proc_tau;
}

void lam::collectInitStats(double initMean, double initVar, double initAccRatio) 
{
    mean = initMean;
    vari = initVar;
    acc_ratio = initAccRatio;
}


void lam::collectInitStats(unsigned long init_loop)
{
    mean /= (double) init_loop;
    vari = vari / (double) init_loop - mean * mean;
    acc_ratio = (double) success / (double) init_loop;
}

void lam::updateLam()
{

    double d = (1.0 - acc_ratio) / (2.0 - acc_ratio);
    alpha = 4.0 * acc_ratio * d * d;
}

/*
bool lam::global_frozen()
{
    return (freeze_cnt >= cnt_crit);
}
 */


void lam::updateStats(aState state) 
{
    tau_count ++;
    if (proc_tau == tau_count) {
        tau_count = 0;
        updateSegment(state);
        resetSegmentStats();
    }
}
