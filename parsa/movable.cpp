#include "movable.h"
#include <cstdlib>
#include <iostream>
#include <limits>
using namespace std;
static const double UNINITIALIZED = numeric_limits<double>::max() ;
abstract_param::~abstract_param() // to make compiler happy
{
    return;
}

void movable::set_theta(int id, double theta)
{
	if (id >= 0 && id < nparams) {
		theta_bars[id] = theta;
	}
}

movable::movable(int np)
{
	init(np);
}


void movable::init(int np)
{
	index = -1;

	nparams = np;
	params = new abstract_param*[nparams];

	success = new long[nparams];
	moves = new long[nparams];
	theta_bars = new double[nparams];
	energy = UNINITIALIZED;
	prev_eng = energy;

}

double movable::propose_move()
{
    // check and perform move control here

    
    index ++;
    index %= nparams;

    if (energy == UNINITIALIZED)
    	energy = get_score();

    prev_eng = energy;
    params[index]->generate_tweak(theta_bars[index]);
    moves[index] ++;
    energy = get_score();
    return (energy - prev_eng);
}

void movable::accept_move()
{
    success[index] ++;
}

void movable::reject_move()
{
    params[index]->restore_tweak();
    energy = prev_eng;
    //cout << energy << " " << get_score() << endl;
}



movable::~movable()
{
    delete[] success;
    delete[] moves;
    delete[] theta_bars;
    delete[] params;

}

