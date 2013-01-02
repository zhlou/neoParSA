#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <libxml/parser.h>
#include "annealer.h"
#include "movable.h"
#include "utils.h"
#include <stdexcept>

using namespace std;

annealer::annealer(movable *theproblem, xmlNode *root) :
        problem(theproblem), xmlroot(root)
{
    xmlsection = getSectionByName(root, "annealer_input");

    if (xmlsection == NULL) {
        throw runtime_error(string("Error: fail to find section annealer_input"));
    }
    s = 1.0 / getPropDouble(xmlsection, "init_T");
    lambda = getPropDouble(xmlsection, "lambda");
    init_loop = getPropInt(xmlsection, "init_loop");
    is_init = false;

    step_cnt = 0;

    rnd_seed = time(NULL);
    energy = problem->get_score();

    cout << "The initial energy is " << energy << endl;
}

double annealer::loop()
{
    if (!is_init)
        initMoves();
    while (!frozen()) {
        resetSegmentStats();
        do {
            updateStep(move());
            updateS();
            step_cnt++;
        } while (inSegment());
        updateSegment();

    }
    cout << "Annealing stopped at s = " << s << endl << "Total steps is "
            << step_cnt << endl;
    return problem->get_score();
}

double annealer::initMoves()
{
    for (int i = 0; i < init_loop; i ++) {
        updateInitStep(move());
    }
    initStats();
    is_init = true;
    return problem->get_score();
}

inline void annealer::updateInitStep(bool accept)
{
    updateStep(accept);
}

void annealer::resetSegmentStats()
{
}

bool annealer::move()
{
    double delta, crit, ran_n;
    delta = problem->propose_move();
    crit = exp(-s * delta);
    ran_n = (double) rand_r(&rnd_seed) / RAND_MAX;
    if ((delta <= 0.0) || crit > ran_n) {
        problem->accept_move();
        energy += delta;
        return true;
    } else {
        problem->reject_move();
        return false;
    }
}

annealer::~annealer()
{

}
