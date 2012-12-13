#include "movable.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <limits>
using namespace std;
static const double UNINITIALIZED = numeric_limits<double>::max();
const double movable::theta_min = 0.;

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

movable::movable(int np, xmlNode *root)
{
    init(np);
    xmlNode *section = NULL;
    if (root != NULL) {
        section = root->children;
        while (section != NULL) {
            if (!xmlStrcmp(section->name, (xmlChar *)"move"))
                break;
            section = section->next;
        }

    }
    if (section == NULL) {
        move_gain = 0.03;
        move_interval = 100;
    } else {
        xmlChar *prop = NULL;
        if ((prop = xmlGetProp(section, (xmlChar *)"gain")) != NULL) {
            move_gain = strtod((char *)prop, NULL);
            xmlFree(prop);
            prop = NULL;
        } else
            move_gain = 0.03;
        if ((prop = xmlGetProp(section, (xmlChar *)"interval")) != NULL) {
            move_interval = atoi((char *)prop);
            xmlFree(prop);
            prop = NULL;
        } else
            move_interval = 100;

    }
}

void movable::init(int np)
{
    index = -1;

    nparams = np;
    params = new abstract_param*[nparams];

    success = new long[nparams];
    moves = new long[nparams];
    theta_bars = new double[nparams];
    for (int i = 0; i < nparams; i++) {
        success[i] = 0;
        moves[i] = 0;
        theta_bars [i] = 1.0;
    }
    energy = UNINITIALIZED;
    prev_eng = energy;
    sweep = 0;

}

double movable::propose_move()
{
    // check and perform move control here


    index++;
    index %= nparams;
    if (index == 0) {
        sweep++;
        if (sweep % move_interval == 0)
            move_control();
    }



    if (energy == UNINITIALIZED)
        energy = get_score();

    prev_eng = energy;
    params[index]->generate_tweak(theta_bars[index]);
    moves[index]++;
    energy = get_score();
    return (energy - prev_eng);
}

void movable::accept_move()
{
    success[index]++;
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

void movable::move_control()
{
    for (int i = 0; i < nparams; i++) {
        double acc_ratio = (double) success[i] / (double) moves[i];
        double x = log(theta_bars[i]);
        cout << i << "\t" << acc_ratio;
        x += move_gain * (acc_ratio - 0.44);
        theta_bars[i] = exp(x);
        if (theta_bars[i] < theta_min)
            theta_bars[i] = theta_min;
        cout << "\t" << theta_bars[i];
        success[i] = 0;
        moves[i] = 0;
    }
    cout << endl;
}
