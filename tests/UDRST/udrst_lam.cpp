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
#include <libxml/parser.h>
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
    char *docname = argv[optind];
    xmlDocPtr xmldoc = xmlParseFile(docname);
    xmlNodePtr xmlroot = xmlDocGetRootElement(xmldoc);

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
        rst.write_section((xmlChar *)"output");
        rst_sa.writeResult();
        xmlSaveFormatFile(docname, xmldoc,1);
    }
    xmlFreeDoc(xmldoc);
    xmlCleanupParser();
    return 0;
}
