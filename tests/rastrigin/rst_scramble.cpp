/*
 * rst_scramble.cpp
 *
 *  Created on: Jul 9, 2014
 *      Author: zhlou
 */

#include <iostream>
#include <unistd.h>
#include <stdexcept>
#include <string>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

#include "rastrigin.h"
#include "unirandom.h"


int main(int argc, char **argv)
{
    const char *sectionname = "rastrigin";
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
                << " [-x tour_section_name] [-f output_name] input_file" << std::endl;
        return -1;
    }
    std::string xmlname(argv[optind]);
    if (NULL == outname)
        outname = xmlname.c_str();
    ptree pt;
    read_xml(xmlname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &xmlroot = pt.begin()->second;

    unirandom rnd;
    rastrigin rst(xmlroot, rnd);
    for (int i = 0; i < rst.getDimension(); ++i) {
        rst.generateMove(i, rastrigin::VAR_MAX * 2.0 * rnd.random() - 1.0);
    }
    std::cout << "Final score is " << rst.get_score() << std::endl;
    rst.write_section(xmlroot, sectionname);
    boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
    write_xml(xmlname, pt, std::locale(), settings);

    return 0;
}
