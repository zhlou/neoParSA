/*
 * lam.cpp
 *
 *  Created on: Dec 8, 2012
 *      Author: zhlou
 */

#include "lam.h"
#include "utils.h"
#include <cmath>
using namespace std;

lam::lam(movable *theproblem, xmlNode *root) :
        annealer(theproblem, root)
{
    xmlNode *section = getSectionByName(root, "lam");
    if (section == NULL)
        throw runtime_error(string("Error: fail to find section lam"));
    proc_tau = getPropInt(section, "tau");
    init_loop = getPropInt(section, "init_loop");


    proc_tau = 100; // TODO for now. should be an input from xml later
    init_loop = 1000;
    resetSegmentStats();
    for (int i = 0; i < init_loop; i++) {
        if (move())
            success ++;
        mean += energy;
        vari += energy * energy;
    }
    mean /= (double)init_loop;
    vari = vari / (double)init_loop - mean * mean;
    double sd = sqrt(vari);
    fit_mean = new invLinearFit(0.995, mean, s, vari/(mean * mean));
    fit_sd = new invLinearFit(0.99, sd, s, sd/mean);
    acc_ratio = (double)success / (double)init_loop;
    double d = (1.0 - acc_ratio) / (2.0 - acc_ratio);
    alpha = 4.0 * acc_ratio * d * d;
    resetSegmentStats();
    freeze_crit = 0.0001;
    old_energy = energy;
    freeze_cnt = 0;
    cnt_crit = 3;

}

lam::~lam()
{
    delete fit_mean;
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
        freeze_cnt ++;
    else
        freeze_cnt = 0;
    old_energy = energy;
    return (freeze_cnt >= cnt_crit);
}

void lam::updateStep(bool is_accept)
{
    if (is_accept)
        success++;
    double estimate_mean = fit_mean->getEstimate(s);
    double d = energy - estimate_mean;
    mean += energy;
    vari += d * d;
}

void lam::updateS()
{
    double estimate_sd = fit_sd->getEstimate(s);
    double d = s * estimate_sd;
    s += lambda * alpha / (d * d * estimate_sd);
}

void lam::resetSegmentStats()
{
    acc_ratio = 0.;
    success = 0;
    vari = 0.;
    mean = 0.;
}

bool lam::inSegment()
{
    if ((step_cnt % proc_tau) == 0) {
        return false;
    } else {
        return true;
    }
}

void lam::updateSegment()
{
    mean /= proc_tau;
    vari /= proc_tau;
    acc_ratio = (double)success / proc_tau;
    updateLam();
    resetSegmentStats();
}

void lam::resetLam()
{
    fit_mean->reset();
    fit_sd->reset();
}

void lam::updateLam()
{
    fit_mean->update(1.0/mean, s);
    fit_sd->update(1.0/sqrt(vari),s);
    double d = (1.0 - acc_ratio) / (2.0 - acc_ratio);
    alpha = 4.0 * acc_ratio * d * d;
}
