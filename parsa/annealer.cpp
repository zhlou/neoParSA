#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <libxml/parser.h>
#include "annealer.h"
#include "movable.h"

using namespace std;

annealer::annealer(movable *theproblem, xmlNode *root) :
        xmlroot(root)
{
    xmlsection = xmlroot->children;
    xmlChar *prop;
    problem = theproblem;
    while (xmlsection != NULL) {
        if (!xmlStrcmp(xmlsection->name, (xmlChar *) "annealer_input"))
            break;
        xmlsection = xmlsection->next;
    }
    if (xmlsection == NULL) {
        throw 2;
    }
    prop = xmlGetProp(xmlsection, (xmlChar *) "init_T");
    if (prop == NULL) {
        throw 3;
    }
    s = 1.0 / atof((char *) prop);
    xmlFree(prop);
    prop = NULL;
    prop = xmlGetProp(xmlsection, (xmlChar *) "lambda");
    if (prop == NULL) {
        throw 3;
    }
    lambda = atof((char *) prop);
    xmlFree(prop);

    step_cnt = 0;

    rnd_seed = time(NULL);
    energy = problem->get_score();

    cout << "The initial energy is " << energy << endl;
}

double annealer::loop()
{
    double delta, crit, ran_n;

    while (!frozen()) {
        do {
            delta = problem->propose_move();
            crit = exp(-s * delta);
            ran_n = (double) rand_r(&rnd_seed) / RAND_MAX;
            if ((delta <= 0.0) || crit > ran_n) {
                problem->accept_move();
                energy += delta;
                updateStep(true, delta);
            } else {
                problem->reject_move();
                updateStep(false, 0.);
            }
            cool_s();
            step_cnt++;
        } while (inSegment());
        updateSegment();

    }
    cout << "Annealing stopped at s = " << s << endl << "Total steps is "
            << step_cnt << endl;
    return problem->get_score();
}

annealer::~annealer()
{

}
