/*********************************************************************************
*                                                                                *
*     scalefactor.h                                                             *
*                                                                                *
*     Contains a simple structure for scale factors                              *
*                                                                                *
*********************************************************************************/

#ifndef SCALEFACTOR_H
#define SCALEFACTOR_H

#include "parameter.h"

#include <map>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>

using namespace std;

class ScaleFactor
{
private:
  /* I use an affine ( ax + b ) transformation for scaling data. I generally
  assume that b=0, but some users set offset lower or are more or less 
  conservative in background removal, so realistically b could be in the range
  -10 to 10. */
  param_ptr A;
  param_ptr B;
  
  string name;
  
public:
  // Constructors
  ScaleFactor();
  ScaleFactor(ptree& pt);
  
  // Getters
  void getParameters(param_ptr_vector& p);
  string& getName() { return name; }
  
  // Method
  double scale(double x);
  double unscale(double x);
  
  // I/O
  void read(ptree& pt);
  void write(ptree& pt);
};

typedef boost::shared_ptr<ScaleFactor> scale_factor_ptr;

class ScaleFactorContainer
{
private:
  map<string, scale_factor_ptr> scales;
  
public:
  // Constructors
  ScaleFactorContainer();
  ScaleFactorContainer(ptree& pt);
  
  // Getters
  void             getParameters(param_ptr_vector& p);
  scale_factor_ptr getScaleFactor(string name) { return scales[name]; }
  
  // I/O
  void read(ptree& pt);
  void write(ptree& pt);
};

typedef boost::shared_ptr<ScaleFactorContainer> scale_factors_ptr;

#endif
