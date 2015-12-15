/*
 * rst_expHold.cpp
 *
 *  Created on: Feb 8, 2014
 *      Author: zhlou
 */
#include <iostream>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <string>

#include <unistd.h>
#include <libgen.h>
#include <getopt.h>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

#include "annealer.h"
#include "move/feedbackMove.h"
#include "unirandom.h"
#include "criCount.h"
#include "dynDebug.h"
#include "expHold.h"
#include "udrst.h"

int main(int argc, char **argv)
{
    bool isprolix = false;
    bool issteplog = true;
    bool isequil = false;
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
        while ( (c = getopt_long(argc, argv, "EpL", long_options, &optIndex)) != -1) {
            switch(c) {
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
                isequil = true;
                break;
            case 'L':
                issteplog = false;
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
    std::string docname(argv[optind]);
    std::string filename(docname);
    ptree pt;
    read_xml(filename, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &docroot = pt.begin()->second;
    unirandom rnd;
    udrst rst(docroot, rnd);
    expHold::Param scheParam(docroot);
    criCount::Param frozenParam(docroot);
    annealer<udrst, expHold, criCount, feedbackMove>*rst_sa
            = new annealer<udrst, expHold, criCount, feedbackMove>
                  (rst, rnd, scheParam, frozenParam, docroot);
    string basename(docname);
    size_t sz = basename.size();
    basename.resize(sz-4);

    string outprefix = basename;
    if (isprolix) {
        rst_sa->setProlix(file, (outprefix + ".prolix").c_str());
    }

    if (iscoollog)
        rst_sa->setCoolLog(file,(basename + ".log").c_str());

    if (issteplog) {
        rst_sa->setStepLog(file, (outprefix + ".steplog").c_str());
    }

    if (isequil) {
        rst_sa->initMovesOnly();
    } else if (saveInitState) {
        rst_sa->initMoves();
        rst_sa->saveUnifiedInitState(saveStatePrefix);
    } else {
        if (readInitState) {
            rst_sa->readUnifiedInitState(readStatePrefix);
        }
        rst_sa->loop();
        std::cout << "The final energy is " << rst.get_score() << std::endl;
        rst.write_section(docroot, "output");
        rst_sa->writeResult(docroot);
        boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
        write_xml(filename, pt, std::locale(), settings);
    }
    delete rst_sa;

    return 0;


}
