/* 
 * File:   transcpp_expHold.cpp
 * Author: zhlou
 *
 * Created on November 18, 2014, 11:05 AM
 */

#include "flags.h"
#include "pwm.h"
#include "TF.h"
#include "gene.h"
#include "conc.h"
#include "twobit.h"
#include "organism.h"
#include <fstream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <unistd.h>
#include <libxml/parser.h>
#include "annealer.h"
#include "move/feedbackMove.h"
#include "unirandom.h"
#include "expHold.h"
#include "tempCount.h"
#include "dynDebug.h"

using namespace std;

using boost::property_tree::ptree;


int mode_verbose;

/*
 * 
 */
int main(int argc, char** argv) 
{
    bool isprolix = false;
    bool issteplog = false;
    std::string binname(basename(argv[0]));
    try {
        char c;
        while ((c = getopt(argc, argv, "pl")) != -1) {
            switch (c) {
            case 'l':
                issteplog = true;
                break;
            case 'p':
                isprolix = true;
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
        std::cerr << "Usage: " << binname << " [ -x section_name ] input_file"
                << std::endl;
        return -1;
    }
    mode_verbose = 0;
    string xmlname(argv[optind]);
    fstream infile(xmlname.c_str());

    ptree pt;
    read_xml(infile, pt, boost::property_tree::xml_parser::trim_whitespace);

    ptree & input = pt.get_child("Input");

    Organism embryo(input);
    //embryo.printParameters(cerr);

    unirand48 rnd;
    //rnd.setSeed(getpid());

    xmlDoc *doc = xmlParseFile(xmlname.c_str());
    xmlNode *docroot = xmlDocGetRootElement(doc);
    expHold::Param scheParam(docroot);
    tempCount::Param frozenParam(docroot);
    annealer<Organism, expHold, tempCount, feedbackMove>
            annealer(embryo, rnd, scheParam, frozenParam, docroot);
    if (isprolix)
        annealer.setProlix(file, (xmlname + ".prolix").c_str());
    annealer.setCoolLog(file, (xmlname + ".log").c_str());
    if (issteplog)
        annealer.setStepLog(file, (xmlname+".steplog").c_str());



    cerr << "The energy is " << embryo.get_score() << endl;
    annealer.loop();
    cerr << "The energy is " << embryo.get_score() << " after loop" << endl;
    xmlFreeDoc(doc);

    //embryo.printParameters(cerr);

    ptree output, anneal_output;
    embryo.write("Output", output);
    annealer.ptreeGetResult(output.get_child("Output"));
    //ptree& opt = output.get_child("Output");
    output.put_child("Output.anneal_output", anneal_output);
    boost::property_tree::xml_writer_settings<char> settings(' ', 2);
    write_xml_element(infile, basic_string<ptree::key_type::value_type>(), output, -1, settings);

    return 0;
}

