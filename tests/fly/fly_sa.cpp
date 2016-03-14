/*
 * fly_sa.cpp
 *
 *  Created on: Mar 1, 2013
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
#include "lam.h"
#include "criCount.h"
#include "dynDebug.h"
#include "fly.h"

using namespace std;

int main(int argc, char **argv)
{
    if (argc <= 1) {
        cerr << "Missing input files" << endl;
        return 1;
    }
    char c;
    bool equil = false;
    while ( (c = getopt(argc, argv, "E")) != -1) {
    	switch(c) {
    	case 'E':
    		equil = true;
    		break;
    	}
    }
    char *docname = argv[optind];
    ptree pt;
    read_xml(docname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &docroot = pt.begin()->second;
    unirand48 rnd;
    fly_params flyParams = readFlyParams(docroot);
    fly theFly(flyParams);
    // feedbackMove<fly, debugSTD> fly_problem(theFly, rnd, docroot);
    lam::Param scheParam(docroot);
    criCount::Param criCntParam(docroot);
    annealer<fly, lam, criCount, feedbackMove>
        fly_sa(theFly,rnd, scheParam, criCntParam, docroot);
    fly_sa.setCoolLog(file, (flyParams.infile_name+".log").c_str());
    fly_sa.setProlix(file, (flyParams.infile_name+".prolix").c_str());
    if (equil)
    	fly_sa.initMoves();
    else {
        cout << "The initial energy is " << theFly.get_score() << endl;
        fly_sa.loop();
        cout << "The final energy is " << theFly.get_score() << endl;
        theFly.writeAnswer("eqparms");
        fly_sa.writeResult(docroot);
        boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
        write_xml(docname, pt, std::locale(), settings);
    }



    return 0;
}
