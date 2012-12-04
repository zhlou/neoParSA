/*
 * rastrigin_problem.cpp
 *
 *  Created on: Dec 3, 2012
 *      Author: zhlou
 */
#include "rastrigin_problem.h"

variable::variable(rastrigin *rst, int in_idx) :
		therst(rst), idx(in_idx)
{
	seed = rst->get_seed();
	is_restorable = false;
	prev_x = 0;
}


void variable::restore_tweak()
{
	if (is_restorable) {
		therst->set_param(idx, prev_x);
		is_restorable = false;
	}
}

void variable::generate_tweak(double theta_bar)
{
	double x;
	x = prev_x = therst->get_param(idx);
	double uniform = rand_r(seed) * 2.0 / RAND_MAX - 1.0;
	if (uniform >= 0)
		x -= theta_bar * log(abs(uniform));
	else
		x += theta_bar * log(abs(uniform));
	if (x > rastrigin::VAR_MAX)
		x = rastrigin::VAR_MAX;
	else if (x < rastrigin::VAR_MIN)
		x = rastrigin::VAR_MIN;
	therst->set_param(idx, x);
	is_restorable = true;
}

rastrigin_problem::rastrigin_problem(rastrigin *rst_problem) :
		movable(rst_problem->get_dimension()), therst(rst_problem)
{
	for (int i = 0; i < nparams; i++) {

		params[i] = new variable(therst, i);
	}
}

double rastrigin_problem::get_score(){
	return therst->value();
}

rastrigin_problem::~rastrigin_problem()
{
	for (int i = 0; i < nparams; i++) {
		delete params[i];
	}
}
