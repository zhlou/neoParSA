/*
 * rst_scramble.cpp
 *
 *  Created on: Jul 9, 2014
 *      Author: zhlou
 */

#include <iostream>
#include <unistd.h>
#include <stdexcept>
#include <libxml/parser.h>

#include "normSchwefel.h"
#include "unirandom.h"


int main(int argc, char **argv)
{
    const char *sectionname = "Schwefel";
    const char *outname = NULL;

    try {
        char c;
        while ((c = getopt(argc, argv, "x:f:")) != -1) {
            switch(c) {
            case 'f':
                outname = optarg;
                break;
            case 'x':
                sectionname = optarg;
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
        std::cerr << "Usage: " << "rst_scramble"
                << " [-x section_name] [-f output_name] input_file" << std::endl;
        return -1;
    }
    const char *xmlname = argv[optind];
    if (NULL == outname)
        outname = xmlname;
    xmlDocPtr xmldoc = xmlParseFile(xmlname);
    xmlNodePtr xmlroot = xmlDocGetRootElement(xmldoc);

    unirandom rnd;
    normSchwefel nsch(xmlroot, rnd);
    nsch.scramble();
    std::cout << "Final score is " << nsch.get_score() << std::endl;
    nsch.writeSection(BAD_CAST sectionname);
    xmlSaveFormatFile(outname, xmldoc, 1);
    xmlFreeDoc(xmldoc);
    xmlCleanupParser();

    return 0;
}


