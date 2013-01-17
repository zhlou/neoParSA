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

}

double annealer::initMoves()
{

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

}

annealer::~annealer()
{

}
