/*
 * fly_expHold.cpp
 *
 *  Created on: Jun 13, 2013
 *      Author: zhlou
 */

#include <iostream>
#include <unistd.h>
#include <libxml/parser.h>

#include "annealer.h"
#include "move/feedbackMove.h"
#include "unirandom.h"
#include "expHold.h"
#include "tempCount.h"
#include "fly.h"

using namespace std;

int main(int argc, char **argv)
{
    if (argc <= 1) {
        cerr << "Missing input files" << endl;
        return 1;
    }
    char c;
    bool isprolix = false;
    bool equil = false;
    while ( (c = getopt(argc, argv, "Ep")) != -1) {
        switch(c) {
        case 'E':
            equil = true;
            break;
        case 'p':
            isprolix = true;
            break;
        default:
            cerr << "Unrecognized option: " << c << endl;
            return 1;
        }
    }
    char *docname = argv[optind];

    xmlDoc *doc = xmlParseFile(docname);
    xmlNode *docroot = xmlDocGetRootElement(doc);
    if (docroot == NULL) {
        cerr << "Input incorrect" << endl;
        return 2;
    }

    unirand48 rnd;
    fly_params flyParams = readFlyParams(docroot);
    fly theFly(flyParams);

    expHold::Param scheduleParam(docroot);
    tempCount::Param tmpCntParam(docroot);
    annealer<fly, expHold, tempCount, feedbackMove>
        fly_expHold(theFly, rnd, scheduleParam, tmpCntParam, docroot);
    fly_expHold.setStepLog(file, (flyParams.infile_name+".steplog").c_str());
    if (isprolix) {
        fly_expHold.setProlix(file,(flyParams.infile_name+".prolix").c_str());
    }
    if (equil)
        fly_expHold.initMovesOnly();
    else {
        fly_expHold.loop();
        theFly.writeAnswer("eqparms");
        fly_expHold.writeResult();
        xmlSaveFormatFile(docname, doc, 1);
    }

    xmlFreeDoc(doc);
    xmlCleanupParser();
    return 0;
}
