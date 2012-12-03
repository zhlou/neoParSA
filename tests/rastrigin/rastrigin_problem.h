/*
 * rastrigin_problem.h
 *
 *  Created on: Dec 3, 2012
 *      Author: zhlou
 */

#ifndef RASTRIGIN_PROBLEM_H_
#define RASTRIGIN_PROBLEM_H_

#include <libxml/parser.h>

class variable : public abstract_param {
    public:
        double x;
        void init(unsigned int *in_seed);
        void generate_tweak(double theta_bar);
        void restore_tweak();
        void set_idx(int i);
    private:
        unsigned int *seed;
        double prev_x;
        bool is_restorable;
        int idx;
};

class rastrigin_problem : public movable
{
public:
	rastrigin_problem(rastrigin *rst_problem);
	double get_score();
	~rastrigin_problem();
private:
	rastrigin *therst;
	variable *vars;

};

#endif /* RASTRIGIN_PROBLEM_H_ */
