/*********************************************************************************
*                                                                                *
*     gene.cpp                                                                   *
*                                                                                *
*     Contains class definition for genes                                        *
*                                                                                *
*********************************************************************************/

#include "gene.h"

#include <boost/foreach.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/range/adaptor/map.hpp>

# define foreach_ BOOST_FOREACH

using namespace boost::adaptors;

/********************************   Gene    *************************************/


/*    Constructors    */


Gene::Gene() {}

Gene::Gene(string name, string head, int left, int right, int t, promoter_ptr p, string& seq)
{
  gname       = name;
  header      = head;
  sequence    = string2int(seq);
  right_bound = right;
  left_bound  = left;
  promoter    = p;
  tss         = t;
}


/*    Getters   */
  
const vector<int> & Gene::getSequence() { return(sequence); }
  
const string Gene::getSequenceString() const { return(int2string(sequence)); }

const string & Gene::getName() const { return(gname); }

int Gene::getRightBound() const { return(right_bound); }
    
int Gene::getLeftBound() const {return(left_bound); }

int Gene::length() const { return(sequence.size()); }


/*    Setters   */

void Gene::setSequence(vector<int> s)
{
  sequence = s;
  right_bound = 0; // default bound is 0
  left_bound = - length();
}

void Gene::setSequence(string s)
{
  sequence = string2int(s);
  right_bound = 0; // default bound is 0
  left_bound = -length();
}

void Gene::setName(string s) { gname = s; }

void Gene::setRightBound(int r)
{
  right_bound = r;
  left_bound = r-length();
}

void Gene::setLeftBound(int l)
{
  left_bound = l;
  right_bound = l+length();
}


/*    Methods   */


/*    Print   */

void Gene::print(ostream& os) const
{
  os << ">" << gname 
     << "|" << right_bound
     << "|" << length() << endl;
  os << int2string(sequence) << endl;
}

void Gene::write(ostream& os) const
{
  ptree pt;
  write(pt);
  boost::property_tree::xml_writer_settings<char> settings(' ', 2);
  write_xml_element(os, basic_string<ptree::key_type::value_type>(), pt, -1, settings);
}

void Gene::write(ptree& pt) const
{
  ptree& gene = pt.add("Gene","");
  gene.put("<xmlattr>.name",       gname);
  gene.put("<xmlattr>.header",     header);
  gene.put("<xmlattr>.left_bound", left_bound);
  gene.put("<xmlattr>.right_bound",right_bound);
  gene.put("<xmlattr>.TSS",        tss);
  gene.put("<xmlattr>.promoter",   promoter->getName());
}
  

/******************************   GeneContainer    ***********************************/
      
/*    Constructor   */

GeneContainer::GeneContainer() {}

GeneContainer::GeneContainer(ptree& pt, promoters_ptr p)
{
  promoters = p;
  read(pt);
}


/*    Setters   */



void GeneContainer::add(gene_ptr gene) { genes.push_back(gene); }


/*    Getters   */

Gene& GeneContainer::getGene(int index) {return *genes[index]; }

Gene& GeneContainer::getGene(const string& n) 
{
  int ngenes = genes.size();
  for (int i=0; i<ngenes; i++)
  {
    if (genes[i]->getName() == n)
      return(*genes[i]);
  }
  cerr << "ERROR: getGene() could not find gene with name " << n << endl;
  exit(1);
  return(*genes[0]);
}

gene_ptr GeneContainer::getGeneptr(int index) {return genes[index]; }

gene_ptr GeneContainer::getGeneptr(const string& n) 
{
  int ngenes = genes.size();
  for (int i=0; i<ngenes; i++)
  {
    if (genes[i]->getName() == n)
      return(genes[i]);
  }
  cerr << "ERROR: getGene() could not find gene with name " << n << endl;
  exit(1);
  return(genes[0]);
}


int GeneContainer::size() const { return genes.size(); }


/*    I/O    */


void GeneContainer::read(ptree & pt)
{
  ptree& genes_node = pt.get_child("Genes");
  foreach_(ptree::value_type const& source, genes_node)
  {
    if (source.first != "Source") continue;

    string name = source.second.get<string>("<xmlattr>.name");
    string file = source.second.get<string>("<xmlattr>.file");
    string type = source.second.get<string>("<xmlattr>.type");
    if (type == "twobit")
    {
      twobit_ptr tmp(new TwoBit(file));
      twobits[name] = tmp;
      readTwoBitGenes( (ptree&) source.second, name);
    }
    else if (type == "fasta")
    {
      fasta_ptr tmp(new Fasta(file));
      fastas[name] = tmp;
      readFastaGenes( (ptree&) source.second, name);
    }
    else
    {
      cerr << "ERROR: unrecognized sequence source file type " << type << endl;
      exit(1);
    }
  }
}

void GeneContainer::readTwoBitGenes(ptree& gene_nodes, string& twobit_name)
{
  foreach_(ptree::value_type const& gene_node, gene_nodes)
  {
    if (gene_node.first != "Gene") continue;
    
    bool include = gene_node.second.get<bool>("<xmlattr>.include", true);
    if (!include) continue;
    
    int error = numeric_limits<int>::max();
    
    string gname  = gene_node.second.get<string>("<xmlattr>.name");
    string header = gene_node.second.get<string>("<xmlattr>.header");
    
    int right_bound = gene_node.second.get<int>("<xmlattr>.right_bound", error); 
    int left_bound  = gene_node.second.get<int>("<xmlattr>.left_bound",  error);
    int tss         = gene_node.second.get<int>("<xmlattr>.TSS",         error);        
    
    if (right_bound == error || left_bound == error || tss == error)
    {
      cerr << "ERROR: right_bound, left_bound, and TSS must be set to positive integers when reading from twobit files" << endl;
      exit(1);
    }
    
    string promoter_name  = gene_node.second.get<string>("<xmlattr>.promoter");
    promoter_ptr promoter = promoters->getPromoter(promoter_name);
    
    string sequence = twobits[twobit_name]->getSequence(header, left_bound, right_bound);
    gene_ptr tmp(new Gene(gname, header, left_bound, right_bound, tss, promoter, sequence));
    genes.push_back(tmp);
    string type("twobit");
    gene_map[type][twobit_name].push_back(tmp);
  }
}

void GeneContainer::readFastaGenes(ptree& gene_nodes, string& fasta_name)
{
  foreach_(ptree::value_type const& gene_node, gene_nodes)
  {
    if (gene_node.first != "Gene") continue;
    
    bool include = gene_node.second.get<bool>("<xmlattr>.include", true);
    if (!include) continue;
    
    int error = numeric_limits<int>::max();
    
    string gname  = gene_node.second.get<string>("<xmlattr>.name");
    string header = gene_node.second.get<string>("<xmlattr>.header");
    
    string& sequence = fastas[fasta_name]->getSeq(header);
    
    int right_bound = gene_node.second.get<int>("<xmlattr>.right_bound", error); 
    int left_bound  = gene_node.second.get<int>("<xmlattr>.left_bound",  error);
    // unless otherwise specified, we assume TSS is at position -1
    int tss         = gene_node.second.get<int>("<xmlattr>.TSS",         -1); 
    
    // if no left or right bound is specified, make seq immediately left of tss
    if (right_bound == error && left_bound == error)
    {
      right_bound = -1; 
      left_bound  = right_bound - sequence.length();
    } 
    else if (right_bound == error && left_bound != error)
    {
      right_bound = left_bound + sequence.length();
    }
    else // we ignore left_bound if right bound is set 
    {
      left_bound  = right_bound - sequence.length();
    }
    
    
    string promoter_name  = gene_node.second.get<string>("<xmlattr>.promoter");
    promoter_ptr promoter = promoters->getPromoter(promoter_name);
    
    gene_ptr tmp(new Gene(gname, header, left_bound, right_bound, tss, promoter, sequence));
    genes.push_back(tmp);
    string type("fasta");
    gene_map[type][fasta_name].push_back(tmp);
  }
}
        
      
void GeneContainer::write(ptree& pt) 
{
 
  ptree& genes_node = pt.add("Genes","");
  foreach_(string type, gene_map | map_keys)
  {
    map<string, gene_ptr_vector>& submap = gene_map[type];
    foreach_(string source_name, submap | map_keys)
    {
      ptree& source_node = genes_node.add("Source","");
  
      string file_name;
      if (type == "twobit")
        file_name = twobits[source_name]->getFileName();
      else if (type == "fasta")
        file_name = fastas[source_name]->getFileName();
      
      source_node.put("<xmlattr>.name", source_name);
      source_node.put("<xmlattr>.file", file_name);
      source_node.put("<xmlattr>.type", type);
      
      int ngenes = gene_map[type][source_name].size();
      for (int i=0; i<ngenes; i++)
        gene_map[type][source_name][i]->write(source_node);
    }
  }   
}

void GeneContainer::write(ostream& os) 
{
  ptree pt;
  write(pt);
  boost::property_tree::xml_writer_settings<char> settings(' ', 2);
  write_xml_element(os, basic_string<ptree::key_type::value_type>(), pt, -1, settings);
}

