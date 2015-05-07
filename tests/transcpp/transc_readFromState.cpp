/* 
 * File:   transc_readFromState.cpp
 * Author: zhlou
 *
 * Created on May 4, 2015, 4:41 PM
 */
#include "pwm.h"
#include "TF.h"
#include "gene.h"
#include "datatable.h"
#include "twobit.h"
#include "organism.h"
#include "mode.h"
#include "utils.h"
#include <fstream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

using namespace std;
using boost::property_tree::ptree;
/*
 * 
 */
int main(int argc, char** argv)
{
    char *stateName = NULL;
    char c;
    while ((c = getopt(argc, argv, "r:")) != -1) {
        switch(c) {
        case 'r':
            stateName = optarg;
            break;
        }
    }
    string xmlname(argv[optind]);
    fstream infile(xmlname.c_str());
    ptree pt;
    read_xml(infile, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree& root_node = pt.get_child("Root");
    ptree& mode_node = root_node.get_child("Mode");
    ptree& input_node = root_node.get_child("Input");
    mode_ptr mode(new Mode(xmlname, mode_node));
    Organism embryo(input_node, mode);
    
    cout << "The initial energy is " << embryo.get_score() << endl;
    
    int bufSize = embryo.getStateSize();
    char * stateBuf = new char[bufSize];
    ifstream stateFile(stateName, ios::in | ios::binary);
    if (!stateFile) {
        cerr << "Fail to open file " << stateName << endl;
        exit(1);
    }
    stateFile.read(stateBuf, bufSize);
    stateFile.close();
    embryo.deserialize(stateBuf);
    cout << "The energy after reading from file is " << embryo.get_score() << endl;
    

    embryo.write("Output", root_node);
    boost::property_tree::xml_writer_settings<char> settings(' ', 2);
    write_xml(xmlname, pt, std::locale(), settings);
    

    return 0;
}

