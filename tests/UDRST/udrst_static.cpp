/*
 *       Filename:  udrst_static.cpp
 *       Created :  09/21/2015 05:18:57 PM
 *       Author  :  Zhihao Lou
 */

#include <iostream>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <string>

#include <unistd.h>
#include <libgen.h>
#include <getopt.h>

#include "annealer.h"
#include "move/feedbackMove.h"
#include "unirandom.h"
#include "tempCount.h"
#include "dynDebug.h"
#include "schedule/staticLam.h"
#include "udrst.h"

int main(int argc, char **argv)
{
    bool isprolix = false;
    bool issteplog = true;
    bool isequil = false;
    int saveInitState = 0;
    int readInitState = 0;
    int optIndex;
    std::string saveStatePrefix;
    std::string readStatePrefix;

    struct option long_options[] = {
        {"save-state", 1, &saveInitState, 1},
        {"read-state", 1, &readInitState, 1},
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
    ptree pt;
    read_xml(docname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &docroot = pt.begin()->second;

    unirandom rnd;
    udrst rst(docroot, rnd);
    staticLam::Param scheParam(docroot);
    tempCount::Param frozenParam(docroot);
    annealer<udrst, staticLam, tempCount, feedbackMove>*rst_sa
            = new annealer<udrst, staticLam, tempCount, feedbackMove>
                  (rst, rnd, scheParam, frozenParam, docroot);
    string basename(docname);
    size_t sz = basename.size();
    basename.resize(sz-4);

    string outprefix = basename;
    if (isprolix) {
        rst_sa->setProlix(file, (outprefix + ".prolix").c_str());
    }

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
        write_xml(docname, pt, std::locale(), settings);
    }
    delete rst_sa;

    return 0;


}
