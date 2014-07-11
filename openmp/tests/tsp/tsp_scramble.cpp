/*
 * tsp_scramble.cpp
 *
 *  Created on: Jun 2, 2014
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
	const char *outname = NULL;
    try {
        char c;
        while ((c = getopt(argc, argv, "x:f:")) != -1) {
            switch(c) {
            case 'f':
                outname = optarg;
                break;
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
                << " [-x tour_section_name] [-f output_name] input_file" << endl;
        return -1;
    }
    const char *xmlname = argv[optind];
    if (NULL == outname)
        outname = xmlname;
    xmlDocPtr xmldoc = xmlParseFile(xmlname);
    xmlNodePtr xmlroot = xmlDocGetRootElement(xmldoc);

    tsp theTSP(xmlroot);
    try {
    	theTSP.read_tour(xmlroot, tourname);
    } catch (std::exception &ex) {
    	std::cerr << ex.what() << std::endl;
    	std::cerr << "Use rank ID tour instead" << std::endl;
    }
    std::cout << "The initial tour length is " << theTSP.get_score()
    		<< std::endl;
    size_t n = theTSP.get_ncities();
    for (size_t i = 0; i < 2*n; ++i) {
    	theTSP.generateMove(0,(double)n);
    }
    std::cout << "Final tour length is " << theTSP.get_score()
    		<< std::endl;
    theTSP.write_tour(xmlroot, tourname);
    xmlSaveFormatFile(outname, xmldoc, 1);
    xmlFreeDoc(xmldoc);
    xmlCleanupParser();

	return 0;
}



