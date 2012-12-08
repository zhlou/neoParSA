/*
 * lam.cpp
 *
 *  Created on: Dec 8, 2012
 *      Author: zhlou
 */

#include "lam.h"

lam::lam(movable *theproblem, xmlNode *root) :
        annealer(theproblem, root)
{
    // TODO Auto-generated constructor stub
    proc_tau = 100; // for now

}

lam::~lam()
{
    // TODO Auto-generated destructor stub
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
