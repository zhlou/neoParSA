/*
 * fly_sa.cpp
 *
 *  Created on: Mar 1, 2013
 *      Author: zhlou
 */

#include <iostream>
#include <libxml/parser.h>
#include <unistd.h>

#include "annealer.h"
#include "feedbackMove.h"
#include "unirandom.h"
#include "lam.h"
#include "dynDebug.h"
#include "fly.h"

using namespace std;

int main(int argc, char **argv)
{
    if (argc <= 1) {
        cerr << "Missing input files" << endl;
        return 1;
    }
    char c;
    bool equil = false;
    while ( (c = getopt(argc, argv, "E")) != -1) {
    	switch(c) {
    	case 'E':
    		equil = true;
    		break;
    	}
    }
    char *docname = argv[optind];
    xmlDoc *doc = xmlParseFile(docname);
    xmlNode *docroot = xmlDocGetRootElement(doc);
    if (docroot == NULL) {
        cerr << "Input incorrect" << endl;
        return 2;
    }
    unirandom rnd;
    fly_params flyParams = readFlyParams(docroot);
    fly theFly(flyParams);
    // feedbackMove<fly, debugSTD> fly_problem(theFly, rnd, docroot);
    lam schedule(docroot);
    annealer<fly, lam, feedbackMove>
        fly_sa(theFly,&rnd, docroot);
    if (equil)
    	fly_sa.initMoves();
    else {
        cout << "The initial energy is " << theFly.get_score() << endl;
        fly_sa.loop();
        cout << "The final energy is " << theFly.get_score() << endl;
        theFly.writeAnswer("eqparam");
        xmlSaveFormatFile(docname, doc, 1);
    }


    return 0;
}

