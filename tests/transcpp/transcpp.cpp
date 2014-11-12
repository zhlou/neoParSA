/*********************************************************************************
*                                                                                *
*     transc.cpp                                                                 *
*                                                                                *
*     Constains the main method for fitting the transcription model              *
*                                                                                *
*********************************************************************************/

#include "flags.h"
#include "pwm.h"
#include "TF.h"
#include "gene.h"
#include "nucleus.h"
#include "conc.h"
#include "twobit.h"
#include "organism.h"
#include <fstream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <unistd.h>
#include <libxml/parser.h>
#include "annealer.h"
#include "feedbackMove.h"
#include "unirandom.h"
#include "lam.h"
#include "criCount.h"
#include "dynDebug.h"

using boost::property_tree::ptree;


int mode_verbose;

int main(int argc, char* argv[])
{

  mode_verbose = 0;
  string xmlname(argv[1]);
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
  lam::Param scheParam(docroot);
  criCount::Param criCntParam(docroot);
  annealer<Organism, lam, criCount, feedbackMove>
    fly_sa(embryo,rnd, scheParam, criCntParam, docroot);
  fly_sa.setCoolLog(file, (xmlname+".log").c_str());
    fly_sa.setProlix(file, (xmlname+".prolix").c_str());

  
  cerr << "The energy is " << embryo.get_score() << endl;
  	fly_sa.initMoves();
  cerr << "The energy is " << embryo.get_score() << " after initMoves" << endl;
  fly_sa.loop();
  cerr << "The energy is " << embryo.get_score() << " after loop" << endl;
  fly_sa.writeResult();
  xmlFreeDoc(doc);
  
  //embryo.printParameters(cerr);
  
  ptree output;
  embryo.write("Output", output);
  boost::property_tree::xml_writer_settings<char> settings(' ', 2);
  write_xml_element(infile, basic_string<ptree::key_type::value_type>(), output, -1, settings);

  return 0;
}
