/*
 * rst_scramble.cpp
 *
 *  Created on: Jul 9, 2014
 *      Author: zhlou
 */

#include <iostream>
#include <unistd.h>
#include <stdexcept>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

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
    ptree pt;
    std::string xmlname(argv[optind]);
    read_xml(xmlname, pt, boost::property_tree::xml_parser::trim_whitespace);
    if (NULL == outname)
        outname = xmlname.c_str();
    ptree &xmlroot = pt.begin()->second;

    unirandom rnd;
    normSchwefel nsch(xmlroot, rnd);
    nsch.scramble();
    std::cout << "Final score is " << nsch.get_score() << std::endl;
    nsch.writeSection(xmlroot, sectionname);
    boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
    write_xml(outname, pt, std::locale(), settings);

    return 0;
}
