#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "annealer.h"
#include "movable.h"

using namespace std;

annealer::annealer(movable *theproblem)
{
    problem = theproblem;
    s = 0.001;
    reject_cnt = 0;
    rnd_seed = time(NULL);
    cout << "The initial energy is "<< problem->get_score() << endl;
}

double annealer::loop()
{
    double delta, crit, ran_n;
    while (!frozen()) {
        delta = problem->propose_move();
        crit = exp(-s * delta);
        ran_n = (double)rand_r(&rnd_seed)/RAND_MAX;
        if ((delta <=0.0) || crit > ran_n) {
            problem->accept_move();
            reject_cnt = 0;
        } else {
            problem->reject_move();
            reject_cnt += 1;
        }
        s *= 1.01;
    }
    cout << "Annealing stoped at s = " << s << endl;
    return problem->get_score();
}

bool annealer::frozen()
{
    const unsigned max_rej = 500;
    return (reject_cnt > max_rej);
}

annealer::~annealer()
{

}
