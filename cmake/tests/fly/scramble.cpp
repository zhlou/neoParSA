/*
 * scramble.cpp
 *
 *  Created on: Mar 14, 2013
 *      Author: zhlou
 */
#include <iostream>
#include <stdexcept>
#include <unistd.h>
#include "fly.h"

int main(int argc, char **argv)
{
    if (argc <= 1) {
        throw std::runtime_error("Missing input files");
    }
    string section;
    char c;
    while ( (c = getopt(argc, argv, "x:")) != -1 ) {
        switch (c) {
        case 'x':
            section = optarg;
            break;
        }
    }
    char *docname = argv[optind];
    xmlDoc *doc = xmlParseFile(docname);
    xmlNode *root = xmlDocGetRootElement(doc);
    fly_params flyParams = readFlyParams(root, "input");
    if (!section.empty())
        flyParams.section_title = section;
    fly theFly(flyParams);
    theFly.scramble();
    theFly.writeAnswer(flyParams.section_title.c_str());
    return 0;


}



