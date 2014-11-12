/*********************************************************************************
*                                                                                *
*     interactions.cpp                                                           *
*                                                                                *
*     Contains classes or structures for tf-tf interactions                      *
*                                                                                *
*********************************************************************************/

/* Interactions all have an actor, a set of targets, a distance function and 
parameters, and a parameter governing their effect. For now I will keep separate 
structures for quenching, cooporativity, coactivation, etc. There is no actor 
function because each actor will own the structure */


/************************    Distance Functions   *******************************/

/*  All distance functions take a distance and return a value between 0 and 1 as
a function of that distance. They accept a vector of doubles as input, that way 
they can be pointed to by any of the interactions */

typedef int (distptr*)(double, vector<double>);

distptr GetDistFunc(string func_name)
{
  if (func_name == "Uniform")
    return Uniform;
  else
  {
    string err;
    err << "GetDistFunc could not find function with name " << func_name;
    error(err);
    return void;
  }
}
    

double Uniform(double distance, double max)
{
  if (abs(distance) <= max)
    return 1;
  else
    return 0;
}

double Trapezoid(double distance, double a, double b)
{
  if (distance < 0)
    distance = -distance; 
  
  if (distance <= a)
    return 1;
  else if (distance < b)
  {
    double x   = distance      - a;
    double max = b - a;
    return (1 - x/max);
  } else
    return 0;
}

double Linear(double distance, double max)
{
  if (distance < 0)
    distance = -distance; 
  
  if (distance < max)
    return (1 - distance/max);
  else
    return 0;
}

// a decaying sin wave
double Sin(double distance, period, offset, max)
{
  const double pi = 3.1415926535897;
  double a      = 2*pi/period;
  
  if (distance < 0)
    distance = -distance; 
  
  if (distance < max)
    return sin(a*distance+offset)+1)*(1-distance/max)/2;
  else
    return 0;
}
    
  

/***************************    Cooporativity   *********************************/


Cooporativity::Cooporativity(
  string t, 
  double v, 
  string dfunc,
  vector<double> dparams)
{
  target     = t;
  distparams = dparams;
  distfunc   = GetDistFunc(dfunc);
  distparams = dparams;
}

Cooporativity::distFunc(double dist)
{
  return distfunc(double, dparams);
}


