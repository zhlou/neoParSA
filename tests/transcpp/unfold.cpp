/*********************************************************************************
*                                                                                *
*     unfold.cpp                                                                 *
*                                                                                *
*     Reads output and allows user to print various data                         *
*                                                                                *
*********************************************************************************/
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
#include <getopt.h>

#include <libxml/parser.h>
#include "annealer.h"
#include "feedbackMove.h"
#include "unirandom.h"
#include "lam.h"
#include "criCount.h"
#include "dynDebug.h"

using boost::property_tree::ptree;


int mode_verbose;

static const char *optString = "hi:s:";

static const struct option longOpts[] = {
    { "help",        no_argument,       NULL, 'h' },
    { "input-file",  required_argument, NULL, 'i' },
    { "section",     required_argument, NULL, 's' },
    { "sites",       no_argument,       NULL,  0  },
    { "occupancy",   no_argument,       NULL,  0  },
    { "modeocc",     no_argument,       NULL,  0  },
    { "effocc",      no_argument,       NULL,  0  },
    { "subgroups",   no_argument,       NULL,  0  },
    { "scores",      required_argument, NULL, 's' },
    { "section",     required_argument, NULL, 's' },
    { "gene",        required_argument, NULL,  0 },
    { "tf",          required_argument, NULL,  0 },
    { "rate",        no_argument,       NULL,  0 },
    { "data",        no_argument,       NULL,  0 },
    { "params",      no_argument,       NULL,  0 },
    { "check-scale", no_argument,       NULL,  0 },
    { "score",       no_argument,       NULL,  0 }
};

void display_usage()
{
  cerr << endl << "\t Usage" << endl << endl
       << "\t unfold [options] -i [infile]" << endl << endl
       << "\t Options" << endl
       << "\t --help    [-h]   print this message" << endl
       << "\t --section [-s]   use section of input file (default Output)" << endl
       << "\t --score          prints the score of the entire fit" << endl
       << "\t --sites          prints binding sites" << endl
       << "\t --occupancy      prints binding site fractional occupancy" << endl
       << "\t --modeocc        prints the occupancy split into tf modes" << endl
       << "\t --effocc         prints the effective (activating) occupancy" << endl
       << "\t --subgroups      prints subgroups" << endl
       << "\t --scores         prints pwm scores" << endl
       << "\t --rate           prints rate for each gene" << endl
       << "\t --data           prints rate data for each gene" << endl
       << "\t --params         prints the parameter table" << endl
       << "\t --check-scale    checks that the scale function works with the scoring function" << endl
       << "\t --gene    [name] prints only for gene with name" << endl
       << "\t --tf name [name] prints only for tf with name" << endl << endl;
  exit(1);
}
       

int main(int argc, char* argv[])
{
  int opt = 0;
  int longIndex = 0;
  string section_name("Output");
  
  string gene_name("");
  string tf_name("");
  
  bool rate       = false;
  bool sites      = false;
  bool data       = false;
  bool score      = false;
  bool occupancy  = false;
  bool modeocc    = false;
  bool effocc     = false;
  bool subgroups  = false;
  bool params     = false;
  bool checkscale = false;
  
  string infile_name;
  
  opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  while(opt != -1)
  {
    switch (opt)
    {
      case 'h':
        display_usage();
        break;
      case 'i':
        infile_name = optarg;
        break;
      case 's':
        section_name = optarg;
        break;
      case 0: 
        if (longOpts[longIndex].name == "sites")
          sites = true;
        if (longOpts[longIndex].name == "score")
          score = true;
        if (longOpts[longIndex].name == "rate")
          rate = true;
        if (longOpts[longIndex].name == "data")
          data = true;
        if (longOpts[longIndex].name == "gene")
          gene_name = optarg;
        if (longOpts[longIndex].name == "tf")
          tf_name = optarg;
        if (longOpts[longIndex].name == "occupancy")
          occupancy = true;
        if (longOpts[longIndex].name == "modeocc")
          modeocc = true;
        if (longOpts[longIndex].name == "effocc")
          effocc = true;
        if (longOpts[longIndex].name == "subgroups")
          subgroups = true;
        if (longOpts[longIndex].name == "section")
          section_name = optarg;
        if (longOpts[longIndex].name == "params")
          params = true;
        if (longOpts[longIndex].name == "check-scale")
          checkscale = true;
        
        break;
    }
    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  }
  
  if (infile_name == "")
    display_usage();
  
  ifstream infile(infile_name.c_str());

  ptree pt;
  read_xml(infile, pt, boost::property_tree::xml_parser::trim_whitespace);
  ptree& section = pt.get_child(section_name);

  Organism embryo(section);

  
  tfs_ptr   tfs   = embryo.getTFs();
  genes_ptr genes = embryo.getGenes();
  
  if (sites)
  {
    if      (gene_name != "" && tf_name != "")
    {
      Gene& gene = genes->getGene(gene_name);
      TF&   tf   = tfs->getTF(tf_name);
      embryo.printSites(gene, tf, cout);
    } 
    else if (gene_name != "" && tf_name == "")
    {
      Gene& gene = genes->getGene(gene_name);
      embryo.printSites(gene, cout);
    }
    else if (gene_name == "" && tf_name != "")
    {
      TF&   tf   = tfs->getTF(tf_name);
      embryo.printSites(tf, cout);
    } else
      embryo.printSites(cout);
  }
  
  if (rate)
  {
    embryo.printRate(cout);
  }
  if (checkscale)
  {
    embryo.checkScale(cout);
  }
  if (data)
    embryo.printRateData(cout);
  if (score)
  {
    embryo.printScore(cout);
  }
  if (params)
    embryo.printParameters(cout);

  if (occupancy)
  {
    if (gene_name == "")
    {
      cerr << "ERROR: must specify a gene to print occupancy!" << endl;
      exit(1);
    }
    Gene& gene = genes->getGene(gene_name);
    embryo.printOccupancy(gene, cout);
  }
  if (modeocc)
  {
    if (gene_name == "")
    {
      cerr << "ERROR: must specify a gene to print occupancy!" << endl;
      exit(1);
    }
    Gene& gene = genes->getGene(gene_name);
    embryo.printModeOccupancy(gene, cout);
  }
  if (effocc)
  {
    if (gene_name == "")
    {
      cerr << "ERROR: must specify a gene to print occupancy!" << endl;
      exit(1);
    }
    Gene& gene = genes->getGene(gene_name);
    embryo.printEffectiveOccupancy(gene, cout);
  }
  
  
  if (subgroups)
  {
    if (gene_name == "")
    {
      cerr << "ERROR: must specify a gene to print subgroups!" << endl;
      exit(1);
    }
    Gene& gene = genes->getGene(gene_name);
    embryo.printSubgroups(gene, cout);
  }
  
  return 0;
}
