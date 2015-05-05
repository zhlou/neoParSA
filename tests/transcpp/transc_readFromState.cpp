/* 
 * File:   transc_readFromState.cpp
 * Author: zhlou
 *
 * Created on May 4, 2015, 4:41 PM
 */
#include "flags.h"
#include "pwm.h"
#include "TF.h"
#include "gene.h"
#include "conc.h"
#include "twobit.h"
#include "organism.h"
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
    
    ptree &input = pt.get_child("Input");
    Organism embryo(input);
    
    cout << "The initial energy is " << embryo.get_score() << endl;
    
    int bufSize = embryo.getStateSize();
    char * stateBuf = new char[bufSize];
    ifstream stateFile(stateName, ios::in | ios::binary);
    stateFile.read(stateBuf, bufSize);
    stateFile.close();
    embryo.deserialize(stateBuf);
    cout << "The energy after reading from file is " << embryo.get_score() << endl;
    
    ptree output;
    embryo.write("Output", output);
    boost::property_tree::xml_writer_settings<char> settings(' ', 2);
    write_xml_element(infile, basic_string<ptree::key_type::value_type>(), output, -1,
            settings);
    

    return 0;
}

