/*
 * lam.cpp
 *
 *  Created on: Dec 8, 2012
 *      Author: zhlou
 */

#include "lam.h"
#include "utils.h"
#include <string>
#include <stdexcept>
#include <cmath>
#include <limits>
using namespace std;

const double lam::UNINITIALIZED = numeric_limits<double>::max();

// in newer compilers, string constants cast to char * will generate
// warning messages unless done explicitly.

lam::lam(xmlNode *root)
{
    xmlNode *section = getSectionByName(root, (char *)"annealer_input");

    if (section == NULL) {
        throw runtime_error(string("Error: fail to find section annealer_input"));
    }
    init_S = 1.0 / getPropDouble(section, (char *)"init_T");
    lambda = getPropDouble(section, (char *)"lambda");
    init_loop = getPropInt(section, (char *)"init_loop");

    section = getSectionByName(root, (char *)"lam");
    if (section == NULL)
        throw runtime_error(string("Error: fail to find section lam"));
    proc_tau = getPropInt(section, (char *)"tau");

    double memlength_mean = getPropDouble(section, (char *)"memLength_mean");
    double memlength_sd = getPropDouble(section, (char *)"memLength_mean");
    w_mean = 1.0 - proc_tau / (memlength_mean / lambda);
    if (w_mean < 0.)
        w_mean = 0.;
    w_sd = 1.0 - proc_tau / (memlength_sd / lambda);
    if (w_sd < 0.)
        w_sd = 0.;
    freeze_crit = getPropDouble(section, (char *)"criterion");
    cnt_crit = getPropInt(section, (char *)"freeze_cnt");
    freeze_cnt = 0;
    fit_mean = NULL;
    fit_sd = NULL;

    old_energy = energy = UNINITIALIZED;
    alpha = UNINITIALIZED;
    s = init_S;
    step_cnt = 0;
    acc_ratio = 0.;
    success = 0;
    vari = 0.;
    mean = 0.;
    //resetSegmentStats();
}

lam::~lam()
{
    if (fit_mean != NULL)
        delete fit_mean;
    if (fit_sd != NULL)
        delete fit_sd;
}
/*
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
 */

bool lam::frozen()
{
    if (abs(energy - old_energy) < freeze_crit)
        freeze_cnt++;
    else
        freeze_cnt = 0;

    return (freeze_cnt >= cnt_crit);
}

void lam::updateStep(bool accept, double in_energy)
{
    energy = in_energy;
    ++ step_cnt;
    if (accept)
        ++ success;
    double estimate_mean = fit_mean->getEstimate(s);
    double d = energy - estimate_mean;
    mean += energy;
    vari += d * d;
}

double lam::updateS(double)
{
    double estimate_sd = fit_sd->getEstimate(s);
    double d = s * estimate_sd;
    s += lambda * alpha / (d * d * estimate_sd);
    return s;
}

void lam::resetSegmentStats()
{
    acc_ratio = 0.;
    success = 0;
    vari = 0.;
    mean = 0.;
    old_energy = energy;
}

bool lam::inSegment()
{
    if ((step_cnt % proc_tau) == 0) {
        return false;
    } else {
        return true;
    }
}

void lam::updateEstimators()
{
    fit_mean->fullUpdate(1.0 / mean, s);
    fit_sd->fullUpdate(1.0 / sqrt(vari), s);
}

void lam::updateSegment()
{
    collectStats();
    updateEstimators();
    updateLam();
//    resetSegmentStats();
}

void lam::resetLam()
{
    fit_mean->reset();
    fit_sd->reset();
}

void lam::initStats()
{
    collectInitStats();

    double sd = sqrt(vari);

    fit_mean = new invLinearFit(w_mean, mean, s, vari / (mean * mean));
    fit_sd = new invLinearFit(w_sd, sd, s, sd / mean);
    acc_ratio = (double) success / (double) init_loop;
    double d = (1.0 - acc_ratio) / (2.0 - acc_ratio);
    alpha = 4.0 * acc_ratio * d * d;
    //resetSegmentStats(); // we don't need this since it is called automatically
                           // before every segment
}

void lam::updateInitStep(bool accept, double in_energy)
{
    energy = in_energy;
    if (accept)
        success++;
    mean += energy;
    vari += energy * energy;
}

void lam::collectStats()
{
    mean /= proc_tau;
    vari /= proc_tau;
    acc_ratio = (double) success / proc_tau;
}

void lam::collectInitStats()
{
    mean /= (double) init_loop;
    vari = vari / (double) init_loop - mean * mean;
}

double lam::getInitS()
{
    return init_S;
}

int lam::getInitLoop()
{
    return init_loop;
}

void lam::updateLam()
{

    double d = (1.0 - acc_ratio) / (2.0 - acc_ratio);
    alpha = 4.0 * acc_ratio * d * d;
}
