/*********************************************************************************
*                                                                                *
*     interactions.h                                                             *
*                                                                                *
*     Contains classes or structures for tf-tf interactions                      *
*                                                                                *
*********************************************************************************/

# include parameter.h

/* Interactions all have an actor, a set of targets, a distance function and 
parameters, and a parameter governing their effect. For now I will keep separate 
structures for quenching, cooporativity, coactivation, etc. There is no actor 
function because each actor will own the structure */


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
  
// interface and contain for distance functions and their parameters
class Distance
{
private:
  string func_name;
  vector<Parameter> params;
  
  void read(ptree&pt);
  
public:
  Distance();
  Distance(ptree& pt);
  
  
  
  

/***********************    Interaction Classes   *******************************/


  
  
class Cooporativity
{
private:
  string target;
  vector<double> distparams;
  void (*distfunc)(double, vector<double>);
  double value;
public:
  // constructors
  Cooporativity(string t, 
    double v, 
    string dfunc,
    vector<double> dparams);
  
  //Cooporativity(ptree& pt);
  
  // methods
  double distFunc(double dist);
};
