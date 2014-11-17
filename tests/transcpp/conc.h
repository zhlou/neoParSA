/*********************************************************************************
*                                                                                *
*     conc.h                                                                     *
*                                                                                *
*     Contains class which holds tf concentration at each nucleus                *
*                                                                                *
*********************************************************************************/

#ifndef CONC_H
#define CONC_H

#include <vector>
#include <string>
#include <cstdlib>
#include <map>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <limits>
#include <iostream>

#include "scalefactor.h"

using namespace std;
using boost::property_tree::ptree;



class ConcContainer
{
private: 
  // vectorized data; must check to make sure all rows are the same length!
  map<string, vector<double> > conc;
  map<string, vector<double> > scaled_conc;
  
  vector<string> names;
  vector<int> ids;
  
  scale_factors_ptr scales;
  map<string, scale_factor_ptr> scales_map;
  
  bool hasName(const string & n);
  bool hasID(int id);
public:
  // Constructors
  ConcContainer();
  ConcContainer(ptree& pt, scale_factors_ptr p, string name);
  
  // Getters
  // Missing the data is given the value NA=numeric_limits<double>max();
  double& getConcByID(int id, const string &, bool scaled);
  double& getConcByID(const string &, int id, bool scaled);
  double& getConcByIndex(int idx, const string &);
  double& getConcByIndex(const string &, int idx);
  
  vector<double>& getConc(const string& n) {return scaled_conc[n];}
  
  scale_factor_ptr getScale(const string& n) {return scales_map[n];}
  
  int index2id(int index) { return ids[index]; }
  
  int size();
  
  // Methods
  void scale_data();
  
  // Output
  void read(ptree& pt, scale_factors_ptr p, string name);
  void write(ostream & os, string node_name);
  void write(ptree & pt, string node_name);
};

typedef boost::shared_ptr<ConcContainer> conc_ptr;


























#endif
