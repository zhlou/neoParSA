#include <iostream>
#include <fstream>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

#include "exponential.h"
#include "annealer.h"
#include "move/feedbackMove.h"
#include "tsp.h"
#include "unirandom.h"
#include "rejCount.h"
//#include "debugOut.h"

using namespace std;

int main(int argc, char **argv)
{
    const char *xmlfile = argv[1];
    ptree pt;
    read_xml(xmlfile, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &root = pt.begin()->second;
    tsp test_tsp(root);
    unirandom rnd;
    //feedbackMove<tsp, debugIGNORE> tsp_problem(test_tsp, rnd, root);
    //simpleSchedule schedule(root);
    exponential::Param scheduleParam(root);
    rejCount::Param rejCntParam(root);
    annealer<tsp, exponential, rejCount, feedbackMove>
        tsp_anneal(test_tsp, rnd, scheduleParam, rejCntParam, root);

    //tsp_problem test_problem(&test_tsp);
    //test_tsp.print_route(cout);
    //simpleSchedule tsp_anneal(&test_problem, root);
    cout.precision(16);
    cout << "The finial energy is " << tsp_anneal.loop() << endl;
    cout << "The energy cached is " << test_tsp.get_score() << endl;
    cout << "The real energy is " << test_tsp.calc_tour() << endl;
    cout << test_tsp.print_route();

    //test_tsp.print_route(cout);
    return 0;

}
