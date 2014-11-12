/*********************************************************************************
*                                                                                *
*     parameter.h                                                                *
*                                                                                *
*     Contains class description for parameters                                  *
*                                                                                *
*********************************************************************************/


#ifndef PARAMETER_H
#define PARAMETER_H

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/shared_ptr.hpp>

using namespace std;
using boost::property_tree::ptree;

/* for now paramters will be of type double, but I will consider templating this
in the future */


class Parameter
{
private:
  string param_name;
  string tf_name;
  
  double value;           // the value of the parameter
  double previous_value;  // the last value the parameter had
  
  double lim_low;
  double lim_high;
  
  bool out_of_bounds;     // is it outside of the limits?
  bool anneal;            // are we ever tweaking it?
  
  bool checkLimits();
  
  bool tf_name_set;
  
public:
  // Constructor 
  Parameter();
  Parameter(ptree& pt);
  
  void read(ptree& pt);
  
  // Getters
  double& getValue();
  double getPrevious() const;
  double getLimHigh() const;
  double getLimLow() const;
  
  bool isOutOfBounds();  
  bool isAnnealed();
  
  const string& getParamName();
  const string& getTFName();
  
  bool is_tf_param();
  
  // Setters
  void set(double);
  void setParamName(const string&);
  void setTFName(const string&);
  void setLimits(double low, double high);
  void setAnnealed(bool b) { anneal = b; }
  
  
  // Methods
  void scramble(double);
  void tweak(double delta);
  void restore();
  
  // Output
  void write(ptree& pt) const;
  void print(ostream& os);
};

typedef boost::shared_ptr<Parameter> param_ptr;
typedef vector<param_ptr> param_ptr_vector;






#endif
