/*********************************************************************************
*                                                                                *
*     distance.h                                                                 *
*                                                                                *
*     Contains functions and interfaces for distance functions                   *
*                                                                                *
*********************************************************************************/

# ifndef DISTANCE_H
# define DISTANCE_H

#include "parameter.h"

#include <map>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>



/************************    Distance Functions   *******************************/

/*  All distance functions take a distance and return a value between 0 and 1 as
a function of that distance. They accept a vector of doubles as input, that way 
they can be pointed to by any of the interactions */


// returns 1 if less than or equal to max, else return 0
double Uniform(double distance, double max);

// returns 1 up to a, then linear decay to b
double Trapezoid(double distance, double a, double b);

// linear decay
double Linear(double distance, double max);

// A linearly decaying sine wave
double Sin(double distance, double period, double offset, double max);
  

/************************    Distance Class   ***********************************/

// interface and contain for distance functions and their parameters
class Distance
{
private:
  string name;
  string func_name;
  
  map<string, param_ptr> params;
  
  double max_distance;
  
  boost::function<double (double)> distFunc;
  
  void read(ptree& pt);
  
public:
  Distance();
  Distance(ptree& pt);
  
  double getMaxDistance();
  double getDistFunc(double);
  string getName() { return(name); }
  
  void setParam(string, double);
  void setDistFunc(string funcname);
  
  void print(ostream& os);
  void write(ptree&);
};

typedef boost::shared_ptr<Distance> distance_ptr;

class DistanceContainer
{
private:
  map<string, distance_ptr> distances;
  
  
public:
  DistanceContainer();
  DistanceContainer(ptree& pt);
  
  void add(string n, distance_ptr d) {distances[n] = d;} 
  void add(ptree& pt);
  
  void print(ostream& os);
  void write(ptree&);
  
  distance_ptr getDistance(string name);
  
  double distFunc(string name, double distance);
};
  
  


typedef boost::shared_ptr<DistanceContainer> distances_ptr;





#endif
  

