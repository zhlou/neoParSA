/*
 * tsp_scramble.cpp
 *
 *  Created on: Jun 2, 2014
 *      Author: zhlou
 */

#include <iostream>
#include <unistd.h>
#include <stdexcept>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

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
    ptree pt;
    read_xml(xmlname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &xmlroot = pt.begin()->second;

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
    boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
    write_xml(xmlname, pt, std::locale(), settings);

	return 0;
}
