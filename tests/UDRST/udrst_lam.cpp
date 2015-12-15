/*
 * rst_lam.cpp
 *
 *  Created on: Jul 9, 2014
 *      Author: zhlou
 */

#include <iostream>
#include <string>
#include <stdexcept>
#include <unistd.h>
#include <cstring>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

#include <libgen.h>

#include "unirandom.h"
#include "criCount.h"
#include "lam.h"
#include "move/feedbackMove.h"
#include "annealer.h"
#include "udrst.h"

using namespace std;

int main(int argc, char **argv)
{
    bool isprolix = false;
    bool issteplog = false;
    bool isequil = false;
    bool iscoollog = true;

    std::string section;
    std::string binname(basename(argv[0]));
    try {
        char c;
        while ( (c = getopt(argc, argv, "ECpl")) != -1) {
            switch(c) {
            case 'E':
                isequil = true;
                break;
            case 'C':
                iscoollog = false;
                break;
            case 'l':
                issteplog = true;
                break;
            case 'p':
                isprolix = true;
                break;
            default:
                throw std::runtime_error("Unrecognized option");
            }
        }
        if (argc <= optind) {
            throw std::runtime_error("Missing input file");
        }
    } catch (std::exception &ex) {
        std::cerr << ex.what() << std::endl;
        std::cerr << "Usage: " << binname << " [-E] [-l] [-p] input_file"
                  << std::endl;
        return -1;
    }
    std::string docname(argv[optind]);
    ptree pt;
    read_xml(docname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &xmlroot = pt.begin()->second;


    unirand48 rnd;
    udrst rst(xmlroot, rnd);
    lam::Param scheduleParam(xmlroot);
    criCount::Param frozenParam(xmlroot);
    annealer<udrst, lam, criCount, feedbackMove>
        rst_sa(rst, rnd, scheduleParam, frozenParam, xmlroot);
    string basename(docname);
    size_t sz = basename.size();
    basename.resize(sz-4);
    string outprefix = basename;
    if (isprolix) {
        rst_sa.setProlix(file, (outprefix + ".prolix").c_str());
    }

    if (iscoollog){
        rst_sa.setCoolLog(file,(basename + ".log").c_str());
    }

    if (issteplog) {
        rst_sa.setStepLog(file, (outprefix + ".steplog").c_str());
    }

    if (isequil) {
        rst_sa.initMovesOnly();
    } else {
        cout << "The initial energy is " << rst.get_score() << endl;
        rst_sa.loop();
        cout << "Final energy is " << rst.get_score() << endl;
        rst.write_section(xmlroot, "output");
        rst_sa.writeResult(xmlroot);
        boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
        write_xml(docname, pt, std::locale(), settings);
    }
    return 0;
}
