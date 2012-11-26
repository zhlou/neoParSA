#include <cmath>
#include <cstdlib>
#include <ctime>
#include "rastrigin.h"

using namespace std;
void variable::init(unsigned int *in_seed)
{
	seed = in_seed;
	x = VAR_MAX * ( (rand_r(seed)*2/RAND_MAX) - 1.0);
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



rastrigin::rastrigin (int dimension)
{
    int i;
	nparams = dimension;
    seed = time(NULL);
    params = new abstract_param*[nparams];
    vars = new variable[nparams];
    for (i = 0; i < nparams; i++) {
    	vars[i].init(&seed);
    	params[i] = vars + i;
    }
    init_stats();
    for (i = 0; i < nparams; i++) {
    	set_theta(i, 1.0);
    }

}

void rastrigin::print_solution(ostream& o) const
{
	o << "{" << endl;
	for (int i = 0; i < nparams; i++) {
		o << vars[i].x << endl;
	}
	o << "}" << endl;
}

rastrigin::~rastrigin()
{
	if (vars != NULL) {
		delete []vars;
	}
	if (params != NULL) {
		delete []params;
	}
}
double rastrigin::get_score()
{
    double tot = 0;
    for (int i = 0; i < nparams; i ++) {
        tot += vars[i].x * vars[i].x - 10.0 * cos(2 * M_PI * vars[i].x);
    }
    return (10 * nparams + tot);
}



