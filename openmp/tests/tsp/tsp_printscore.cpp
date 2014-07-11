/*
 * tsp_printscore.cpp
 *
 *  Created on: May 30, 2014
 *      Author: zhlou
 */

#include <iostream>
#include <unistd.h>
#include <stdexcept>
#include <libxml/parser.h>

#include "tsp.h"

int main(int argc, char **argv)
{
    const char *tourname = "tour";
    try {
        char c;
        while ((c = getopt(argc, argv, "x:")) != -1) {
            switch(c) {
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
        std::cerr << "Usage: " << "tsp_printscore"
                << " [-x tour_section_name] input_file" << endl;
        return -1;
    }
    const char *xmlname = argv[optind];
    xmlDocPtr xmldoc = xmlParseFile(xmlname);
    xmlNodePtr xmlroot = xmlDocGetRootElement(xmldoc);

    tsp theTSP(xmlroot);
    try {
        theTSP.read_tour(xmlroot, tourname);
    } catch (std::exception &ex) {
        std::cerr << ex.what() << std::endl;
        exit(0);
    }
    std::cout.precision(16);
    std::cout << "The tour length is " << theTSP.get_score() << endl;
    xmlFreeDoc(xmldoc);
    xmlCleanupParser();
    return 0;


}
