/*********************************************************************************
*                                                                                *
*     scramble.cpp                                                               *
*                                                                                *
*     Scrables the parameters                                                    *
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
  fstream file(xmlname.c_str());


  ptree pt;
  read_xml(file, pt, boost::property_tree::xml_parser::trim_whitespace);
  
  ptree & input = pt.get_child("Input");

  Organism embryo(input);
  embryo.scramble();

  
  ptree output;
  embryo.write("Input",  output);
  boost::property_tree::xml_writer_settings<char> settings(' ', 2);
  write_xml_element(cout, basic_string<ptree::key_type::value_type>(), output, -1, settings);

  return 0;
}
