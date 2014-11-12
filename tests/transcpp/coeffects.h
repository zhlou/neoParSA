/*********************************************************************************
*                                                                                *
*     coeffects.h                                                                *
*                                                                                *
*     The structures for coactivation and corepression.                          *
*                                                                                *
*********************************************************************************/

#ifndef COEFFECTS_H
#define COEFFECTS_H

#include "parameter.h"
#include "distance.h"

#include <map>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>


/***********************    Coeffect Class   *******************************/

class Coeffect
{
private:
  // sources of information
  distances_ptr distances;
  
  // tfs involved
  string actor;
  string target;
  int    coef_idx; // what coef of the target this causes
  
  // parameters
  param_ptr    efficiency;
  distance_ptr dist;
  
  // Orientations
  bool HH;
  bool HT;
  bool TT;
  
public:
  // Constructors
  Coeffect();
  Coeffect(distances_ptr, ptree&);
  
  // Getters
  void   getParameters(param_ptr_vector& p);
  
  string& getActor()  { return actor; }
  string& getTarget() { return target; }
  double  getEfficiency() { return efficiency->getValue(); }
  int     getIdx()        { return coef_idx; }
  
  double distFunc(double d) { return dist->getDistFunc(d); }
  
  // Setters
  void setDist(distance_ptr d)    { dist = d;}
    
  // I/O
  void read(ptree& pt);
  void write(ptree& pt);
  void print(ostream& os);
};
  
typedef boost::shared_ptr<Coeffect> coeffect_ptr;

typedef vector< pair<string,coeffect_ptr> > coeffect_pairs;

/* A container which stores multiple Promoter objects and interfaces with
genes */
class CoeffectContainer
{
private:
  // sources of information
  distances_ptr distances;
  
  vector<coeffect_ptr> coeffects;
  
public:
  //Constructors
  CoeffectContainer();
  CoeffectContainer(distances_ptr, ptree&);
  
  // Getters
  coeffect_pairs getTargets(string tfname);
  void           getParameters(param_ptr_vector& p);
  
  // Setters
  void add(coeffect_ptr c) {coeffects.push_back(c);}
  void setDistances(distances_ptr d) {distances = d;}
  
  // I/O
  void read(ptree& pt, distances_ptr);
  void write(ptree& pt);
  void print(ostream& os);
};
  
typedef boost::shared_ptr<CoeffectContainer> coeffects_ptr;
  
  
  
#endif
