#include "tsp_problem.h"
#include <cstdlib>
#include <ctime>

using namespace std;

tsp_reorder::tsp_reorder(tsp *from_problem)
{
    thetsp = from_problem;
    ncities = thetsp->get_ncities();
    seed = time(NULL);
}

tsp_reorder::~tsp_reorder()
{
    return;
}

void tsp_reorder::generate_tweak(double theta_bar)
{
    int c1, c2;
    c1 = rand_r(&seed) % ncities;
    c2 = rand_r(&seed) % ncities;

    thetsp->step(c1, c2);
}

void tsp_reorder::restore_tweak()
{
    thetsp->roll_back();
}


double tsp_problem::get_score()
{
    return thetsp->cost();
}

tsp_problem::tsp_problem(tsp *ext_tsp)
{
    thetsp = ext_tsp;
    nparams = 1;
    params = new abstract_param*;
    params[0] = new tsp_reorder(thetsp);
    init_stats();

}

tsp_problem::~tsp_problem()
{
    delete params[0];
    delete params;
}
