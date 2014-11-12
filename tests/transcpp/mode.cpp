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

/*    Constructors    */

Mode::Mode() {}

Mode::Mode(ptree& pt) { read(pt); }










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
      
void Mode::read(ptree& pt)
{
  ptree& mode_node = pt.get_child("Mode");
  
  readNode<bool>(  mode_node, string("MultipleSubgroups"), &multiple_subgroups, false    );
  readNode<string>(mode_node, string("OccupancyMethod"),   &occupancy_method,   "dynamic");
  readNode<string>(mode_node, string("ScoreFunction"),     &score_function,     "sse"    );
  readNode<string>(mode_node, string("WeightFunction"),    &weight_function,    "area"   );
  readNode<bool>(  mode_node, string("Unscale"),           &unscale,            true     );
}

void Mode::write(ptree& pt)
{
  ptree& mode_node = pt.add("Mode", "");
  
  ptree& multiple_subgroups_node = mode_node.add("MultipleSubgroups", "");
  ptree& occupancy_method_node   = mode_node.add("OccupancyMethod  ", "");
  ptree& score_function_node     = mode_node.add("ScoreFunction    ", "");
  ptree& weight_function_node    = mode_node.add("WeightFunction   ", "");
  ptree& unscale_node            = mode_node.add("Unscale          ", "");
  
  multiple_subgroups_node.put("<xmlattr>.value", multiple_subgroups);
  occupancy_method_node.put("<xmlattr>.value", occupancy_method);
  score_function_node.put("<xmlattr>.value", score_function);
  weight_function_node.put("<xmlattr>.value", weight_function);
  unscale_node.put("<xmlattr>.value", unscale);
}
