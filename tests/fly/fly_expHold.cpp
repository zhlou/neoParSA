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
#include "feedbackMove.h"
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
    char *docname = argv[1];

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
        fly_expHold(theFly, &rnd, scheduleParam, tmpCntParam, docroot);
    fly_expHold.setStepLog(file, (flyParams.infile_name+".steplog").c_str());
    fly_expHold.loop();
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return 0;
}
