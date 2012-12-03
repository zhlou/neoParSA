/*
 * rastrigin_problem.cpp
 *
 *  Created on: Dec 3, 2012
 *      Author: zhlou
 */
#include "rastrigin_problem.h"



void variable::init(unsigned int *in_seed)
{
	seed = in_seed;
	x = VAR_MAX * ( (rand_r(seed)*2.0/RAND_MAX) - 1.0);
	is_restorable = false;
}

void variable::restore_tweak()
{
	if (is_restorable){
		x = prev_x;
		is_restorable = false;
	}
}

void variable::generate_tweak(double theta_bar)
{
    prev_x = x;
	double uniform = rand_r(seed)*2.0/RAND_MAX - 1.0;
    if (uniform >= 0)
        x -= theta_bar * log(abs(uniform));
    else
        x += theta_bar * log(abs(uniform));
    if (x > VAR_MAX)
        x = VAR_MAX;
    else if (x < VAR_MIN)
        x = VAR_MIN;
    is_restorable = true;
}

rastrigin_problem::rastrigin_problem(rastrigin *rst_problem):
		movable(rst_problem->get_dimension()),therst(rst_problem)
{
	vars = new variable[nparams];
	for (int i = 0; i < nparams; i ++) {
		vars[i].set_idx(i);
		params[i] = vars + i;
	}
}
