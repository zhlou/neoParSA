/*********************************************************************************
*                                                                                *
*     promoter.cpp                                                               *
*                                                                                *
*     Contains parameters that are shared among many nuclei (Q,theta)            *
*                                                                                *
*********************************************************************************/

#include "promoter.h"

#include <math.h>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/ref.hpp>
#include <limits>

using namespace boost::adaptors;

# define foreach_ BOOST_FOREACH


/************************    Functions   ****************************************/

double Arrhenius(double M, double max, double theta, double Q)
{
  double exponent = exp(Q*M - theta);
  return max * (exponent) / (1 + exponent);
}


/************************    Promoter Class   ***********************************/

/*    Constructors    */

Promoter::Promoter() {}

Promoter::Promoter(ptree& pt) { read(pt); }


/*    Getters   */

void Promoter::getParameters(param_ptr_vector& p)
{
  foreach_(param_ptr param, params | map_values)
  {
    if (param->isAnnealed())
      p.push_back(param);
  }
}

/*    I/O   */

void Promoter::read(ptree& pt)
{
  name      = pt.get<string>("<xmlattr>.name");
  func_name = pt.get<string>("<xmlattr>.function");
  
  if (func_name == "Arrhenius")
  {
    param_ptr qparam(    new Parameter(pt.get_child("Q")));
    param_ptr maxparam(  new Parameter(pt.get_child("Rmax")));
    param_ptr thetaparam(new Parameter(pt.get_child("Theta")));
    
    qparam->setParamName("Q");
    maxparam->setParamName("Rmax");
    thetaparam->setParamName("Theta");
    
    params["Q"]     = qparam;
    params["Rmax"]  = maxparam;
    params["Theta"] = thetaparam;
    
    rateFunc = bind(&Arrhenius, _1, 
      boost::ref(params["Rmax"]->getValue()), 
      boost::ref(params["Theta"]->getValue()), 
      boost::ref(params["Q"]->getValue()));
  }
  else
  {
    cerr << "ERROR: read promoter could not find function with name " << func_name << endl;
    exit(1);
  }
}
    
void Promoter::write(ptree& pt)
{
  ptree& promoter_node = pt.add("Promoter","");
  
  promoter_node.put("<xmlattr>.name",     name);
  promoter_node.put("<xmlattr>.function", func_name);
  
  if (func_name == "Arrhenius")
  {
    ptree& q_node     = promoter_node.add("Q    ","");
    ptree& max_node   = promoter_node.add("Rmax ","");
    ptree& theta_node = promoter_node.add("Theta","");
    
    params["Q"]->write(q_node);
    params["Rmax"]->write(max_node);
    params["Theta"]->write(theta_node);
  }
}


/*******************    Promoter Container Class   ******************************/


/*    Constructors    */

PromoterContainer::PromoterContainer() {}

PromoterContainer::PromoterContainer(ptree& pt) { read(pt); }


/*    Getters   */

void PromoterContainer::getParameters(param_ptr_vector& p)
{
  foreach_(promoter_ptr promoter, promoters | map_values)
    promoter->getParameters(p);
}


/*    I/O   */

void PromoterContainer::read(ptree& pt)
{
  ptree& promoters_node = pt.get_child("Promoters");
  
  foreach_(ptree::value_type& node, promoters_node)
  {
    if (node.first != "Promoter") continue;
    
    promoter_ptr promoter(new Promoter((ptree&) node.second));
    string& name = promoter->getName();
    promoters[name] = promoter;
  }
}

void PromoterContainer::write(ptree& pt)
{
  ptree& promoters_node = pt.add("Promoters", "");
  
  foreach_(promoter_ptr promoter, promoters | map_values)
  {
    promoter->write(promoters_node);
  }
}
    
