/*********************************************************************************
*                                                                                *
*     parameter.cpp                                                              *
*                                                                                *
*     Contains class description for parameters                                  *
*                                                                                *
*********************************************************************************/


#include "parameter.h"
#include "utils.h"

#include <cmath>


/*    Constructors    */

template< typename T> 
Parameter<T>::Parameter() 
{
  anneal      = false;
  tf_name_set = false;
  seed        = 1;
  param_name  = "not set";
  move_func   = "ResetAll";
  setTypeName();
}

template< typename T> 
Parameter<T>::Parameter(ptree& pt)
{
  seed = 1;
  node = &pt;
  setTypeName();
  read(pt);
}

template< typename T> 
Parameter<T>::Parameter(string name, ptree& pt)
{
  seed       = 1;
  node       = &pt;
  param_name = name;
  setTypeName();
  read(pt);
}

template< typename T >
void Parameter<T>::setTypeName() { warning("typenames must be set to export parameters to R or Matlab"); }

template< >
void Parameter<int>::setTypeName()      { type = "int";      }

template< >
void Parameter<double>::setTypeName()   { type = "double";   }

template< >
void Parameter<Sequence>::setTypeName() { type = "Sequence"; }

template< >
void Parameter<PWM>::setTypeName()      { type = "PWM";      }


/*    Getters   */

template< typename T> 
T& Parameter<T>::getValue() { return value; }

template< typename T> 
T Parameter<T>::getPrevious() const  { return previous_value; }

template< typename T> 
double Parameter<T>::getLimHigh() const { return lim_high; }

template< typename T> 
double Parameter<T>::getLimLow() const { return lim_low; }


/*    Setters   */

template< typename T> 
void Parameter<T>::set(T v)
{
  value = v;
}

template< typename T> 
void Parameter<T>::setLimits(double low, double high)
{
  lim_high = high;
  lim_low  = low;
}

/*    Methods   */

template< typename T> 
bool Parameter<T>::checkLimits()
{
  if (value > lim_high || value < lim_low)
  {
    out_of_bounds = true;
    return true;
  }
  else
  {
    out_of_bounds = false;
    return false;
  }
}

template<> 
bool Parameter<Sequence>::checkLimits()
{
  out_of_bounds = false;
  return false;
}

template<> 
bool Parameter<PWM>::checkLimits()
{
  out_of_bounds = false;
  return false;
}



template< typename T> 
void Parameter<T>::tweak(double delta)
{
  previous_value = value;
  value += delta;
  checkLimits();
}

/*  Naively, to tweak sequence I will simply pick some number of bases, based
on an integer rounding of delta, then mutate them to a random base. I will
use delta as a seed for random  */

template<> 
void Parameter<Sequence>::tweak(double delta)
{
  previous_value = value;
  
  int nedits = ceil(delta); // the number of bases to mutate
  vector<int>& seq = value.getSequence(); // the seq to tweak
  int length = seq.size(); // the last index in the sequence
  
  for (int i=0; i<nedits; i++)
  {
    // output = min + (rand_r(seed) % (int)(max - min + 1))
    int pos  = rand_r(&seed) % length;
    int base = rand_r(&seed) % (int) 4;
    seq[pos] = base;
  }
}

template<> 
void Parameter<PWM>::tweak(double delta)
{
  previous_value = value;
  
  vector<vector<double> > & pwm = value.getPWM();
  int length = pwm.size();
  
  int pos1 = rand_r(&seed) % length;
  int pos2 = rand_r(&seed) % (int) 4;
  
  pwm[pos1][pos2] += delta;
  
  value.setNscore();
}

template< typename T> 
void Parameter<T>::scramble(double rand_uniform)
{
  value = (lim_high - lim_low)*rand_uniform + lim_low;
  stringstream tmp;
  tmp << setprecision(5);
  tmp.str("");
  tmp << value;
  node->put("<xmlattr>.value", tmp.str());
  if (checkLimits())
  {
    stringstream err;
    err << "ERROR: scrambled variable out of bounds" << endl;
    error(err.str());
  }
}

template<> 
void Parameter<Sequence>::scramble(double rand_uniform)
{
  vector<int>& seq = value.getSequence();
  int length = seq.size();
  
  seed = (unsigned int) rand_uniform * 100000;
  
  for (int i=0; i<length; i++)
    seq[i] = rand_r(&seed) % (int) 4;
}

// I am making scrambled pwms according to 
template<> 
void Parameter<PWM>::scramble(double rand_uniform)
{
  vector<vector<double> > & pwm = value.getPWM();
  int length = pwm.size();
  
  seed = (unsigned int) rand_uniform * 100000;
  
  for (int i=0; i<length; i++)
  {
    for (int j=0; j<4; j++)
      pwm[i][j] = -10 + (rand_r(&seed) % (int)(21));
  }   
}
  


/*  Serialize and Deserialize functions  */

template< typename T> 
void Parameter<T>::serialize(void *buf) const
{
  T * dest = static_cast<T *>(buf);
  *dest = value;
}

template< typename T> 
void Parameter<T>::deserialize(void const *buf) 
{
  T const * from = static_cast<T const *>(buf);
  value = *from;
}


template<> 
void Parameter<Sequence>::serialize(void *buf) const
{
  error("serialize not implemented for parameter of type Sequence");
}

template<> 
void Parameter<Sequence>::deserialize(void const *buf)
{
  error("deserialize not implemented for parameter of type Sequence");
}


template<> 
void Parameter<PWM>::serialize(void *buf) const
{
  error("serialize not implemented for parameter of type PWM");
}

template<> 
void Parameter<PWM>::deserialize(void const *buf)
{
  error("deserialize not implemented for parameter of type PWM");
}


template< typename T> 
size_t Parameter<T>::getSize()
{
  return sizeof(T);
}

template<> 
size_t Parameter<Sequence>::getSize()
{
  error("getSize not implemented for parameter of type Sequence");
  return 0;
}

template<> 
size_t Parameter<PWM>::getSize()
{
  error("getSize not implemented for parameter of type PWM");
  return 0;
}













template< typename T> 
void Parameter<T>::restore()
{
  value   = previous_value;
}


/*    I/O   */


template< typename T> 
void Parameter<T>::read(ptree& pt)
{
  node        = &pt;
  tf_name_set = false;
  
  // get the name of the move function to be used, if not found, call reset all.
  move_func    = pt.get<string>("<xmlattr>.move",    string("ResetAll"));
  
  /* // need to let parameter get mode pointer 
  if (mode->getVerbose() >= 1 && (move_func == string("ResetAll") || restore_func == string("ResetAll")))
    cerr << "WARNING: A move or restore function was set to ResetAll. This can be very slow" << endl;
  */
   
  value        = pt.get<T>("<xmlattr>.value");
  anneal       = pt.get<bool>(  "<xmlattr>.anneal");
               
  lim_low      = pt.get<double>("<xmlattr>.lim_low");
  lim_high     = pt.get<double>("<xmlattr>.lim_high");
  
  previous_value = value;
  
  out_of_bounds = checkLimits();
}

template<> 
void Parameter<Sequence>::read(ptree& pt)
{
  stringstream err;
  err << "Reading of sequence is handled by gene.cpp!" << endl;
  error(err.str());
}

template<> 
void Parameter<PWM>::read(ptree& pt)
{
  stringstream err;
  err << "Reading of pwm is handled by TF.cpp!" << endl;
  error(err.str());
}


template< typename T> 
void Parameter<T>::write(ptree& pt, int precision) const
{
  stringstream tmp;
  tmp << setprecision(precision);
  tmp.str("");
  tmp << value;
  pt.put("<xmlattr>.value", tmp.str());
  tmp.str("");
  tmp << lim_low;
  pt.put("<xmlattr>.lim_low",  tmp.str());
  tmp.str("");
  tmp << lim_high;
  pt.put("<xmlattr>.lim_high", tmp.str());
  pt.put("<xmlattr>.anneal", anneal);
  pt.put("<xmlattr>.move", move_func);
}

template<> 
void Parameter<Sequence>::write(ptree& pt, int precision) const
{
  stringstream err;
  err << "ERROR: writing sequence not implemented yet" << endl;
  error(err.str());
}

template<> 
void Parameter<PWM>::write(ptree& pt, int precision) const
{
  stringstream err;
  err << "ERROR: writing PWM not implemented yet" << endl;
  error(err.str());
}

template< typename T> 
void Parameter<T>::printHeader(ostream& os)
{
  os << setprecision(4)
     << setw(10) << "TF"
     << setw(18) << "type"
     << setw(12) << "value"
     << setw(6)  << "tweak"
     << setw(10) << "lim_low"
     << setw(10) << "lim_high"
     << endl;
}

template< typename T> 
void Parameter<T>::print(ostream& os)
{
  os << setprecision(4)
     << setw(10) << tf_name
     << setw(18) << param_name
     << setw(12) << value
     << setw(6)  << anneal
     << setw(12) << lim_low
     << setw(12) << lim_high
     << endl;
}

template<> 
void Parameter<Sequence>::print(ostream& os)
{
  stringstream err;
  err << "ERROR: printing sequence not implemented yet" << endl;
  error(err.str());
}

template<> 
void Parameter<PWM>::print(ostream& os)
{
  stringstream err;
  err << "ERROR: printing PWM not implemented yet" << endl;
  error(err.str());
}

template class Parameter<double>;
template class Parameter<int>;
template class Parameter<Sequence>;
template class Parameter<PWM>;


  
