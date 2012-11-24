#include <cmath>
#include <stdlib>
#include "rastrigin.h"

using namespace std;

void variable::generate_tweek(double theta_bar)
{
    double uniform = rand_r(seed)*2.0/RAND_MAX - 1.0;
    if (uniform >= 0)
        x -= theta_bar * log(abs(uniform));
    else
        x += theta_bar * log(abs(uniform));
    if (x > VAR_MAX)
        x = VAR_MAX;
    else if (x < VAR_MIN)
        x = VAR_MIN;
}

rastrigin::rastrigin (int dimension)
{
    nparams = dimension;
    vars = new variables[nparams];
}

void rastrigin::init_vars(unsigned int *seed)
{
    int i;
    for (i = 0; i < ndim; i++) {
        vars[i] = VAR_MAX * ( (rand_r(seed)*2/RAND_MAX) - 1.0);
    }
}

double value() 
{
    double tot = 0;
    for (i = 0; i < ndim; i ++) {
        tot += vars[i] * vars[i] - 10.0 * cos(2 * M_PI * vars[i]);
    }
    return (10 * ndim + tot);
}
