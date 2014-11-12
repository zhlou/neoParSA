/*********************************************************************************
*                                                                                *
*     fasta.cpp                                                                  *
*                                                                                *
*     contains methods for parsing a fasta file into a map of sequences          *
*                                                                                *
*********************************************************************************/

#include "fasta.h"

#include <boost/foreach.hpp>
#include <boost/range/adaptor/map.hpp>
#include <limits>
#include <map>

#define foreach_ BOOST_FOREACH

using namespace boost::adaptors;

/*    Constructors    */

Fasta::Fasta() {}

Fasta::Fasta(string& fname) 
{
  file_name = fname;
  read(fname); 
}


/*    Getters   */

string& Fasta::getSeq(string& name)
{
  
  if (seqs.find(name) == seqs.end())
  {
    cerr << "ERROR: could not find seq with name " << name << endl;
    exit(1);
  }
  else
    return seqs[name];
  
}

/*    I/O   */

void Fasta::read(string& fname)
{
  string line;
  ifstream file(fname.c_str());
  string seq_name;
  
  if (file.is_open())
  {
    while(getline(file, line))
    {
      if (!line.empty())
      {
        if (line.at(0) == '>') // 
        {
          nseqs++;
          seq_name = line.substr(1,line.length()-1);
        } 
        else
        {
          seqs[seq_name] = seqs[seq_name] + line;
        }
      }
    }
  }
}

void Fasta::print(ostream& os)
{
  foreach_(string seq_name, seqs | map_keys)
  {
    os << '>' << seq_name << endl;
    os << seqs[seq_name] << endl;
  }
}
          
          
          
