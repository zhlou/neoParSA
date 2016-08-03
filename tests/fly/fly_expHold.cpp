/*
 * fly_expHold.cpp
 *
 *  Created on: Jun 13, 2013
 *      Author: zhlou
 */

#include <iostream>
#include <unistd.h>
#include <getopt.h>

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
    bool isprolix = false;
    bool issteplog = false;
    bool equil = false;
    int iscoollog = 0;
    int saveInitState = 0;
    int readInitState = 0;
    int optIndex;
    std::string saveStatePrefix;
    std::string readStatePrefix;

    struct option long_options[] = {
        {"save-state", 1, &saveInitState, 1},
        {"read-state", 1, &readInitState, 1},
        {"cool-log", 0, &iscoollog, 1},
        {0, 0, 0, 0}
    };

    std::string binname(basename(argv[0]));
    try {
        char c;
        while ((c = getopt_long(argc, argv, "Epl", long_options, &optIndex)) != -1) {
            switch (c) {
            case 0:
                switch (optIndex) {
                case 0:
                    saveStatePrefix = optarg;
                    break;
                case 1:
                    readStatePrefix = optarg;
                    break;
                case 2:
                    break;
                default:
                    throw std::runtime_error("Unrecognized option");
                }
                break;
            case 'E':
                equil = true;
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
        std::cerr << "Usage: " << binname << " [ -x section_name ] input_file"
                << std::endl;
        return -1;
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
        fly_sa(theFly, rnd, scheduleParam, tmpCntParam, docroot);
    fly_sa.setStepLog(file, (flyParams.infile_name+".steplog").c_str());
    if (issteplog)
        fly_sa.setStepLog(file, (flyParams.infile_name+".steplog").c_str());
    if (isprolix) {
        fly_sa.setProlix(file,(flyParams.infile_name+".prolix").c_str());
    }
    if (iscoollog)
        fly_sa.setCoolLog(file,(flyParams.infile_name+".log").c_str());

    if (saveInitState) {
        fly_sa.initMoves();
        fly_sa.saveUnifiedInitState(saveStatePrefix);
    } else {
        if (readInitState) {
            fly_sa.readUnifiedInitState(readStatePrefix);
        }
        if (equil)
            fly_sa.initMovesOnly();
        else {
            fly_sa.loop();
            theFly.writeAnswer("eqparms");
            fly_sa.writeResult(docroot);
            boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
            write_xml(docname, pt, std::locale(), settings);
        }
    }
    return 0;
}
