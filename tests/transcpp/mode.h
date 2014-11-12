/*********************************************************************************
*                                                                                *
*     mode.h                                                                     *
*                                                                                *
*     Information for how the program is run. If it is something that will       *                     
*     change the output of the program in any way, it should be here so that     *                       
*     a change in output is directly reflected in the input file. Anything       *                     
*     that does not change output (i.e. verbose level) should be handled by      *                      
*     flags, with the possible exception of modes that change how efficiently    *                        
*     the code runs.                                                             *
*                                                                                *
*********************************************************************************/

#ifndef MODE_H
#define MODE_H

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/shared_ptr.hpp>

using namespace std;
using boost::property_tree::ptree;

class Mode
{
private:
  bool   multiple_subgroups;
  bool   unscale;
  string occupancy_method; // dynamic, combinations, sampling
  string score_function;   // see, cc
  string weight_function;  
  
  template< typename T> 
  void readNode(ptree&, string, T*, T);
  
public:
  // Constructors
  Mode();
  Mode(ptree& pt);
  
  // Getters
  bool   multipleSubgroups() { return multiple_subgroups; }
  string occupancyMethod()   { return occupancy_method;   }
  string scoreFunction()     { return score_function;     }
  string weightFunction()    { return weight_function;    }
  bool   getUnscale()        { return unscale;            }
  
  // Setters
  void setMultipleSubgroups(bool v) { multiple_subgroups = v; }
  void setOccupancyMethod(string v) { occupancy_method   = v; }
  void setScoreFunction(string v)   { score_function     = v; }
  void setWeightFunction(string v)  { weight_function    = v; }
  void setUnscale(bool v)           { unscale            = v; }
  
  // I/O
  void read(ptree& pt);
  void write(ptree& pt);
  
};


typedef boost::shared_ptr<Mode> mode_ptr;














#endif
