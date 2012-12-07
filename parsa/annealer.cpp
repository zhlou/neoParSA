#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <libxml/parser.h>
#include "annealer.h"
#include "movable.h"

using namespace std;

annealer::annealer(movable *theproblem, xmlNode *root):
		xmlroot(root)
{
	xmlNode *section = xmlroot->children;
	xmlChar *prop;
    problem = theproblem;
    while (section != NULL) {
    	if (!xmlStrcmp(section->name, (xmlChar *)"annealer_input"))
    		break;
    	section = section->next;
    }
    if (section == NULL) {
    	throw 2;
    }
    prop = xmlGetProp(section, (xmlChar *)"init_T");
    if (prop == NULL) {
    	throw 3;
    }
    s = 1.0/atof((char *)prop);
    xmlFree(prop);
    prop = NULL;
    prop = xmlGetProp(section, (xmlChar *)"lambda");
    if (prop == NULL) {
    	throw 3;
    }
    lambda = atof((char *)prop);
    xmlFree(prop);

    reject_cnt = 0;
    rnd_seed = time(NULL);

    cout << "The initial energy is "<< problem->get_score() << endl;
}

double annealer::loop()
{
    double delta, crit, ran_n;
    long unsigned step_cnt = 0;
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
        cool_s();
        step_cnt ++;
    }
    cout << "Annealing stopped at s = " << s << endl
         << "Total steps is " << step_cnt << endl;
    return problem->get_score();
}
void annealer::cool_s()
{
	s *= (1 + lambda);
}
bool annealer::frozen()
{
    const unsigned max_rej = 100;
    return (reject_cnt > max_rej);
}

annealer::~annealer()
{

}
