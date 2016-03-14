/*
 * printscore.cpp
 *
 *  Created on: Feb 22, 2013
 *      Author: zhlou
 */
#include <iostream>
#include <stdexcept>
#include <unistd.h>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

#include "fly.h"

using namespace std;

int main(int argc, char **argv)
{
    if (argc <= 1) {
        throw runtime_error("Missing input files");
    }
    string section;
    char c;
    while ( (c = getopt(argc, argv, "x:")) != -1 ) {
        switch (c) {
        case 'x':
            section = optarg;
            break;
        }
    }
    char *docname = argv[optind];
    ptree pt;
    read_xml(docname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &root = pt.begin()->second;
    fly_params flyParams = readFlyParams(root);
    if (!section.empty())
        flyParams.section_title = section;
    fly theFly(flyParams);
    cout.precision(flyParams.ndigits);
    cout << "chisq = " << theFly.get_score() << " rms = " << theFly.get_rms()
            << endl;
    return 0;
}
