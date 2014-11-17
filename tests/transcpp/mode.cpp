/*********************************************************************************
*                                                                                *
*     mode.cpp                                                                   *
*                                                                                *
*     Information for how the program is run. If it is something that will       *                     
*     change the output of the program in any way, it should be here so that     *                       
*     a change in output is directly reflected in the input file. Anything       *                     
*     that does not change output (i.e. verbose level) should be handled by      *                      
*     flags, with the possible exception of modes that change how efficiently    *                        
*     the code runs.                                                             *
*                                                                                *
*********************************************************************************/

#include "mode.h"

#include <boost/foreach.hpp>
#include <boost/range/adaptor/map.hpp>
using namespace boost::adaptors;
# define foreach_ BOOST_FOREACH



/*    Constructors    */

Mode::Mode() {}

Mode::Mode(ptree& pt) { read(pt); }



bool Mode::getWeightFunction(string s)
{
  return weight_functions[s];
}

/*    I/O   */

template<typename T>
void Mode::readNode(ptree& pt, string node_name, T* value, T def)
{
  if (pt.count(node_name))
  {
    ptree& node = pt.get_child(node_name);
    *value = node.get<T>("<xmlattr>.value", def);
  }
  else
    *value = def;
}

void Mode::readWeightFunc(ptree& pt)
{
  if (pt.count("WeightFunction"))
  {
    ptree& weight_func_node = pt.get_child("WeightFunction");
    min_weight = weight_func_node.get<int>("<xmlattr>.min", 0);
    foreach_(ptree::value_type const& v, weight_func_node)
    {
      if (v.first != "<xmlattr>")
      {
        weight_functions[v.first] = v.second.get<bool>("<xmlattr>.value", false);
      }
    }
  }
  else
    min_weight = 0.0;
}
    
    
  
void Mode::read(ptree& pt)
{
  ptree& mode_node = pt.get_child("Mode");
  
  readNode<bool>(  mode_node, string("MultipleSubgroups"), &multiple_subgroups, false    );
  readNode<string>(mode_node, string("OccupancyMethod"),   &occupancy_method,   "dynamic");
  readNode<string>(mode_node, string("ScoreFunction"),     &score_function,     "sse"    );
  
  readWeightFunc(mode_node);
}

void Mode::writeWeightFunc(ptree& pt)
{
  ptree& weight_func_node = pt.add("WeightFunction","");
  weight_func_node.put("<xmlattr>.min",min_weight);
  foreach_(string s, weight_functions | map_keys)
  {
    ptree& func_node = weight_func_node.add(s,"");
    func_node.put("<xmlattr>.value", weight_functions[s]);
  }
}

void Mode::write(ptree& pt)
{
  ptree& mode_node = pt.add("Mode", "");
  
  ptree& multiple_subgroups_node = mode_node.add("MultipleSubgroups", "");
  ptree& occupancy_method_node   = mode_node.add("OccupancyMethod  ", "");
  ptree& score_function_node     = mode_node.add("ScoreFunction    ", "");

  
  multiple_subgroups_node.put("<xmlattr>.value", multiple_subgroups);
  occupancy_method_node.put("<xmlattr>.value", occupancy_method);
  score_function_node.put("<xmlattr>.value", score_function);

  writeWeightFunc(mode_node);
}
