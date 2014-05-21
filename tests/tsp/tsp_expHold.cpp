/*
 * tsp_expHold.cpp
 *
 *  Created on: May 21, 2014
 *      Author: zhlou
 */
#include <iostream>
#include <string>
#include <stdexcept>
#include <unistd.h>
#include <cstring>
#include <libxml/parser.h>

#include "unirandom.h"
#include "expHold.h"
#include "tempCount.h"
#include "feedbackMove.h"
#include "annealer.h"
#include "tsp.h"

using namespace std;

int main(int argc, char **argv)
{
    bool isprolix = false;
    bool issteplog = true;
    bool isequil = false;

    std::string section;
    std::string binname(basename(argv[0]));
    try {
        char c;
        while ( (c = getopt(argc, argv, "EpL")) != -1) {
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
    expHold::Param scheduleParam(xmlroot);
    tempCount::Param frozenParam(xmlroot);
    annealer<tsp, expHold, tempCount, feedbackMove>
        tsp_expHold(theTSP, rnd, scheduleParam, frozenParam, xmlroot);
    string basename(docname);
    size_t sz = basename.size();
    basename.resize(sz-4);
    string outprefix = basename;
    if (isprolix) {
        tsp_expHold.setProlix(file, (outprefix + ".prolix").c_str());
    }

    tsp_expHold.setCoolLog(file,(basename + ".log").c_str());

    if (issteplog) {
        tsp_expHold.setStepLog(file, (outprefix + ".steplog").c_str());
    }

    if (isequil) {
        tsp_expHold.initMovesOnly();
    } else {
        cout << "The initial energy is " << theTSP.get_score() << endl;
        tsp_expHold.loop();
        cout << "Final energy is " << theTSP.get_score() << endl;
        theTSP.write_tour(xmlroot);
        tsp_expHold.writeResult();
        xmlSaveFormatFile(docname, xmldoc,1);
    }
    xmlFreeDoc(xmldoc);
    xmlCleanupParser();
    return 0;
}
