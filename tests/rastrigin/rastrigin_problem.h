/*
 * rastrigin_problem.h
 *
 *  Created on: Dec 3, 2012
 *      Author: zhlou
 */

#ifndef RASTRIGIN_PROBLEM_H_
#define RASTRIGIN_PROBLEM_H_

#include <libxml/parser.h>
#include <cstdlib>
#include "movable.h"

class rastrigin;
class unirandom;
class variable: public abstract_param
{
public:
	variable(rastrigin *rst, int in_idx);
	void generate_tweak(double theta_bar);
	void restore_tweak();
private:
	rastrigin *therst;
	unirandom &rnd;
	double prev_x;
	bool is_restorable;
	int idx;

};

class rastrigin_problem: public movable
{
public:
	rastrigin_problem(rastrigin *rst, xmlNode *root=NULL);
	double get_score();
	~rastrigin_problem();
private:
	rastrigin *therst;

};

#endif /* RASTRIGIN_PROBLEM_H_ */
