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
#include <libxml/parser.h>

#include "annealer.h"
#include "move/feedbackMove.h"
#include "unirandom.h"
#include "tempCount.h"
#include "dynDebug.h"
#include "expHold.h"
#include "normSchwefel.h"

int main(int argc, char **argv)
{
    bool isprolix = false;
    bool issteplog = true;
    bool isequil = false;
    int saveInitState = 0;
    int readInitState = 0;
    int optIndex;
    char *saveStatePrefix = NULL;
    char *readStatePrefix = NULL;

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
    char *docname = argv[optind];
    xmlDoc *doc = xmlParseFile(docname);
    xmlNode *docroot = xmlDocGetRootElement(doc);
    if (docroot == NULL) {
        std::cerr << "Input incorrect" << std::endl;
        return -1;
    }
    unirandom rnd;
    normSchwefel nsch(docroot, rnd);
    expHold::Param scheParam(docroot);
    tempCount::Param frozenParam(docroot);
    annealer<normSchwefel, expHold, tempCount, feedbackMove>*nsch_sa
            = new annealer<normSchwefel, expHold, tempCount, feedbackMove>
                  (nsch, rnd, scheParam, frozenParam, docroot);
    string basename(docname);
    size_t sz = basename.size();
    basename.resize(sz-4);

    string outprefix = basename;
    if (isprolix) {
        nsch_sa->setProlix(file, (outprefix + ".prolix").c_str());
    }

    nsch_sa->setCoolLog(file,(basename + ".log").c_str());

    if (issteplog) {
        nsch_sa->setStepLog(file, (outprefix + ".steplog").c_str());
    }

    if (isequil) {
        nsch_sa->initMovesOnly();
    } else if (saveInitState) {
        nsch_sa->initMoves();
        nsch_sa->saveUnifiedInitState(saveStatePrefix);
    } else {
        if (readInitState) {
            nsch_sa->readUnifiedInitState(readStatePrefix);
        }
        nsch_sa->loop();
        std::cout << "The final energy is " << nsch.get_score() << std::endl;
        nsch.writeSection((xmlChar *)"output");
        nsch_sa->writeResult(docroot);
        xmlSaveFormatFile(docname, doc, 1);
    }
    xmlFreeDoc(doc);
    xmlCleanupParser();
    delete nsch_sa;

    return 0;


}
