/*
 * tsp_lam.cpp
 *
 *  Created on: May 29, 2014
 *      Author: zhlou
 */

#include <iostream>
#include <string>
#include <stdexcept>
#include <unistd.h>
#include <cstring>
#include <libxml/parser.h>
#include <libgen.h>

#include "unirandom.h"
#include "criCount.h"
#include "lam.h"
#include "move/feedbackMove.h"
#include "annealer.h"
#include "tsp.h"

using namespace std;

int main(int argc, char **argv)
{
    bool isprolix = false;
    bool issteplog = true;
    bool isequil = false;

    const char *tourname = "tour";

    std::string section;
    std::string binname(basename(argv[0]));
    try {
        char c;
        while ( (c = getopt(argc, argv, "EpLx:")) != -1) {
            switch(c) {
            case 'E':
                isequil = true;
                break;
            case 'L':
                issteplog = false;
                break;
            case 'p':
                isprolix = true;
                break;
            case 'x':
                tourname = optarg;
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
        std::cerr << "Usage: " << binname << " [-E] [-L] [-p] input_file"
                  << std::endl;
        return -1;
    }
    char *docname = argv[optind];
    xmlDocPtr xmldoc = xmlParseFile(docname);
    xmlNodePtr xmlroot = xmlDocGetRootElement(xmldoc);

    unirand48 rnd;
    tsp theTSP(xmlroot);
    try {
        theTSP.read_tour(xmlroot, tourname);
    } catch (std::exception &ex) {
        std::cerr << ex.what() << std::endl;
        std::cerr << "Using ID rank tour instead" << std::endl;
    }
    lam::Param scheduleParam(xmlroot);
    criCount::Param frozenParam(xmlroot);
    annealer<tsp, lam, criCount, feedbackMove>
        tsp_sa(theTSP, rnd, scheduleParam, frozenParam, xmlroot);
    string basename(docname);
    size_t sz = basename.size();
    basename.resize(sz-4);
    string outprefix = basename;
    if (isprolix) {
        tsp_sa.setProlix(file, (outprefix + ".prolix").c_str());
    }

    tsp_sa.setCoolLog(file,(basename + ".log").c_str());

    if (issteplog) {
        tsp_sa.setStepLog(file, (outprefix + ".steplog").c_str());
    }

    if (isequil) {
        tsp_sa.initMovesOnly();
    } else {
        cout << "The initial energy is " << theTSP.get_score() << endl;
        tsp_sa.loop();
        cout << "Final energy is " << theTSP.get_score() << endl;
        cout << "Actual energy is " << theTSP.calc_tour() << endl;
        theTSP.write_tour(xmlroot, "tour_output");
        tsp_sa.writeResult();
        xmlSaveFormatFile(docname, xmldoc,1);
    }
    xmlFreeDoc(xmldoc);
    xmlCleanupParser();
    return 0;
}
