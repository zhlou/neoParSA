/*
 * fly_expHold.cpp
 *
 *  Created on: Jun 13, 2013
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
#include "expHold.h"
#include "tempCount.h"
#include "fly.h"

using namespace std;

int main(int argc, char **argv)
{
    if (argc <= 1) {
        cerr << "Missing input files" << endl;
        return 1;
    }
    char c;
    bool isprolix = false;
    bool equil = false;
    while ( (c = getopt(argc, argv, "Ep")) != -1) {
        switch(c) {
        case 'E':
            equil = true;
            break;
        case 'p':
            isprolix = true;
            break;
        default:
            cerr << "Unrecognized option: " << c << endl;
            return 1;
        }
    }
    char *docname = argv[optind];
    ptree pt;
    read_xml(docname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &docroot = pt.begin()->second;

    unirand48 rnd;
    fly_params flyParams = readFlyParams(docroot);
    fly theFly(flyParams);

    expHold::Param scheduleParam(docroot);
    tempCount::Param tmpCntParam(docroot);
    annealer<fly, expHold, tempCount, feedbackMove>
        fly_expHold(theFly, rnd, scheduleParam, tmpCntParam, docroot);
    fly_expHold.setStepLog(file, (flyParams.infile_name+".steplog").c_str());
    if (isprolix) {
        fly_expHold.setProlix(file,(flyParams.infile_name+".prolix").c_str());
    }
    if (equil)
        fly_expHold.initMovesOnly();
    else {
        fly_expHold.loop();
        theFly.writeAnswer("eqparms");
        fly_expHold.writeResult(docroot);
        boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
        write_xml(docname, pt, std::locale(), settings);
    }
    return 0;
}
