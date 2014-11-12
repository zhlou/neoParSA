/*********************************************************************************
*                                                                                *
*     parameter.cpp                                                              *
*                                                                                *
*     Contains class description for parameters                                  *
*                                                                                *
*********************************************************************************/


#include "parameter.h"

/* for now paramters will be of type double, but I will consider templating this
in the future */

Parameter::Parameter() {}

Parameter::Parameter(ptree& pt)
{
  read(pt);
}

void Parameter::read(ptree& pt)
{
  tf_name_set = false;
  
  value  = pt.get<double>("<xmlattr>.value");
  anneal = pt.get<bool>("<xmlattr>.anneal"); 
  
  lim_low  = pt.get<double>("<xmlattr>.lim_low");
  lim_high = pt.get<double>("<xmlattr>.lim_high");
  
  previous_value = value;
  
  out_of_bounds = checkLimits();
}

bool Parameter::checkLimits()
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

double& Parameter::getValue() { return value; }

double Parameter::getPrevious() const  { return previous_value; }

double Parameter::getLimHigh() const { return lim_high; }
double Parameter::getLimLow() const { return lim_low; }

bool Parameter::isOutOfBounds() { return out_of_bounds; }

bool Parameter::isAnnealed() { return anneal; }

void Parameter::set(double v)
{
  value = v;
}

void Parameter::setParamName(const string& pname) { param_name = pname; }

void Parameter::setTFName(const string& tfname) 
{ 
  tf_name = tfname; 
  tf_name_set = true;
}

void Parameter::setLimits(double low, double high)
{
  lim_high = high;
  lim_low  = low;
}

bool Parameter::is_tf_param() { return tf_name_set; }

const string& Parameter::getParamName() { return param_name; }
const string& Parameter::getTFName() { return tf_name; }

void Parameter::tweak(double delta)
{
  previous_value = value;
  value += delta;
  checkLimits();
  /*if (out_of_bounds)
  {
    value = previous_value;
  } */
}

void Parameter::scramble(double rand_uniform)
{
  value = (lim_high - lim_low)*rand_uniform + lim_low;
  if (checkLimits())
  {
    cerr << "ERROR: scrambled variable out of bounds" << endl;
    exit(1);
  }
}
  

void Parameter::restore()
{
  value   = previous_value;
}

void Parameter::write(ptree& pt) const
{
  stringstream tmp;
  tmp << setprecision(5);
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
}


void Parameter::print(ostream& os)
{
  os << setprecision(4);
  os << "name:     " << setw(20) << param_name << "  ";
  os << "value:    " << setw(6)  << value      << "  ";
  os << "anneal:   " << setw(6)  << anneal     << "  ";
  os << "lim_low:  " << setw(6)  << lim_low    << "  ";
  os << "lim_high: " << setw(6)  << lim_high   << "  ";
  os << endl;
}
  
