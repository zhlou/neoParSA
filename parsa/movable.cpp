#include "movable.h"
#include <cstdlib>
using namespace std;

movable::movable()
{
    index = -1;
    theta_bars = NULL;
    success = NULL;
    moves = NULL;
    nparams = 0;
    params = NULL;
}

double movable::propose_move()
{
    // check and perform move control here

    
    index ++;
    index %= nparams;

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
}

void movable::init_stats()
{
    success = new long[nparams];
    moves = new long[nparams];
    theta_bars = new double[nparams];
    energy = get_score();
}

movable::~movable()
{
    delete[] success;
    delete[] moves;
    delete[] theta_bars;

}

