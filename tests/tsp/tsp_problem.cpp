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
    do {
        c2 = rand_r(&seed) % ncities;
    } while (c2 == c1);

    thetsp->step(c1, c2);
}

void tsp_reorder::restore_tweak()
{
    thetsp->roll_back();
}


double tsp_problem::get_score()
{
    return thetsp->cost();
	//return thetsp->calc_route();
}

tsp_problem::tsp_problem(tsp *ext_tsp) :movable(1)
{
    thetsp = ext_tsp;

    params = new abstract_param*;
    params[0] = new tsp_reorder(thetsp);


}

tsp_problem::~tsp_problem()
{
    delete params[0];
    delete params;
}
