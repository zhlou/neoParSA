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
#include <boost/optional.hpp>

#include "annealer.h"
#include "move/fixedThetaMove.h"
#include "unirandom.h"
#include "tempCount.h"
#include "dynDebug.h"
#include "expHold.h"
#include "udrst.h"

int main(int argc, char **argv)
{
    bool isprolix = false;
    bool issteplog = true;
    bool isequil = false;
    int startZero = 0;
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
        {"start-zero", 0, &startZero, 1},
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
                    saveStatePrefix = std::string(optarg);
                    break;
                case 1:
                    readStatePrefix = std::string(optarg);
                    break;
                case 2:
                case 3:
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
    std::string xmlname(argv[optind]);
    ptree pt;
    read_xml(xmlname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &pt_root = pt.begin()->second;
    unirandom rnd;
    boost::optional<unsigned int> seed = pt_root.get_optional<unsigned int>("<xmlattr>.seed");
    if (seed) {
        rnd.setSeed(*seed);
    }
    udrst rst(pt_root, rnd);
    if (startZero) {
        for (int i = 0; i < rst.get_dim(); ++i) {
            rst.set_param(i, 0.0);
        }
    }
    expHold::Param scheParam(pt_root);
    tempCount::Param frozenParam(pt_root);
    annealer<udrst, expHold, tempCount, fixedThetaMove>*rst_sa
            = new annealer<udrst, expHold, tempCount, fixedThetaMove>
                  (rst, rnd, scheParam, frozenParam, pt_root);
    string basename(xmlname);
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
        rst.write_section(pt_root, std::string("output"));
        rst_sa->writeResult(pt_root);
        boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
        write_xml(xmlname, pt, std::locale(), settings);
    }
    delete rst_sa;

    return 0;


}
