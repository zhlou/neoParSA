/*
 * printscore.cpp
 *
 *  Created on: Feb 22, 2013
 *      Author: zhlou
 */
#include <iostream>
#include <stdexcept>
#include <libxml/parser.h>
#include "fly.h"

using namespace std;

int main(int argc, char **argv)
{
    if (argc <= 1) {
        throw runtime_error("Missing input files");
    }
    char *docname = argv[1];
    xmlDoc *doc = xmlParseFile(docname);
    xmlNode *root = xmlDocGetRootElement(doc);
    fly_params flyParams = readFlyParams(root);
    fly theFly(flyParams);
    cout << "chisq = " << theFly.get_score() << " rms = " << theFly.get_rms()
            << endl;
    return 0;
}


