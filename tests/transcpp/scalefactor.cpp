/*********************************************************************************
*                                                                                *
*     scalefactor.cpp                                                            *
*                                                                                *
*     Contains a simple structure for scale factors                              *
*                                                                                *
*********************************************************************************/

#include "scalefactor.h"


#include <boost/foreach.hpp>
#include <boost/range/adaptor/map.hpp>

#define foreach_ BOOST_FOREACH

using namespace boost::adaptors;

/************************    ScaleFactor Class   ********************************/

/*    Constructors    */

ScaleFactor::ScaleFactor() 
{
  name = "default";
  param_ptr Aparam(new Parameter());
  param_ptr Bparam(new Parameter());
  A = Aparam;
  B = Bparam;
  A->set(1);
  B->set(0);
  A->setAnnealed(false);
  B->setAnnealed(false);
  A->setLimits(1,1);
  B->setLimits(0,0);
  A->setParamName("ScaleFactor");
  B->setParamName("ScaleFactor");
}

ScaleFactor::ScaleFactor(ptree& pt) 
{
  read(pt);
}


/*    Getters   */

void ScaleFactor::getParameters(param_ptr_vector& p)
{
  if (A->isAnnealed())
    p.push_back(A);
  if (B->isAnnealed())
    p.push_back(B);
}


/*    Methods   */

double ScaleFactor::scale(double x)
{
  return x * A->getValue() + B->getValue();
}

double ScaleFactor::unscale(double x)
{
  return (x - B->getValue())/A->getValue();
}


/*    I/O   */

void ScaleFactor::read(ptree& pt)
{
  name = pt.get<string>("<xmlattr>.name");
  param_ptr Aparam(new Parameter(pt.get_child("A")));
  param_ptr Bparam(new Parameter(pt.get_child("B")));
  A = Aparam;
  B = Bparam;
  A->setParamName("ScaleFactor");
  B->setParamName("ScaleFactor");
}

void ScaleFactor::write(ptree& pt)
{
  ptree& factor_node = pt.add("ScaleFactor", "");
  factor_node.put("<xmlattr>.name", name);
  ptree& A_node = factor_node.add("A", "");
  A->write(A_node);
  ptree& B_node = factor_node.add("B", "");
  B->write(B_node);
}


/*******************    ScaleFactor Container Class   ***************************/


/*    Constructors    */

ScaleFactorContainer::ScaleFactorContainer() {}

ScaleFactorContainer::ScaleFactorContainer(ptree& pt) 
{
  read(pt);
}


/*    Getters   */

void ScaleFactorContainer::getParameters(param_ptr_vector& p)
{
  foreach_(scale_factor_ptr scale, scales | map_values)
    scale->getParameters(p);
}


/*    I/O   */

void ScaleFactorContainer::read(ptree& pt)
{
  ptree& scale_factors_node = pt.get_child("ScaleFactors");
  
  foreach_(ptree::value_type& node, scale_factors_node)
  {
    if (node.first != "ScaleFactor") continue;
    
    scale_factor_ptr scale(new ScaleFactor((ptree&) node.second));
    string& name = scale->getName();
    scales[name] = scale;
  }
  
  scale_factor_ptr scale(new ScaleFactor());
  scales["default"] = scale;
}

void ScaleFactorContainer::write(ptree& pt)
{
  ptree& scales_node = pt.add("ScaleFactors", "");
  
  foreach_(scale_factor_ptr scale, scales | map_values)
  {
    scale->write(scales_node);
  }
}
