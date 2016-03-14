/*
 * fly_pulse.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: zhlou
 */

#include <iostream>
#include <unistd.h>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

#include "annealer.h"
#include "move/feedbackMove.h"
#include "unirandom.h"
#include "oneStep.h"
#include "maxSteps.h"
#include "dynDebug.h"
#include "fly.h"

using namespace std;

int main(int argc, char **argv)
{
    if (argc <= 1) {
        cerr << "Missing input files" << endl;
        return 1;
    }
    char *docname = argv[1];
    ptree pt;
    read_xml(docname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &docroot = pt.begin()->second;

    unirand48 rnd;
    fly_params flyParams = readFlyParams(docroot);
    fly theFly(flyParams);
    oneStep::Param scheduleParam(docroot);
    maxSteps::Param frozenParam(docroot);
    annealer<fly, oneStep, maxSteps, feedbackMove>
        fly_pulse(theFly, rnd, scheduleParam, frozenParam, docroot);
    fly_pulse.setStepLog(file, (flyParams.infile_name+".steplog").c_str());
    fly_pulse.loop();
    return 0;
}
