/*********************************************************************************
*                                                                                *
*     utils.cpp                                                                  *
*                                                                                *
*     Contains some general functions useful for a number of different classes   *
*                                                                                *
*********************************************************************************/

# include "utils.h"
# include <cstdlib>
# include <cmath>

vector<int> string2int(const string& s)
{
  char t;
  int i, len;
  len = s.length();
  
  vector<int> out;
  
  for (i=0; i<len; i++)
  {
    t=toupper(s[i]);
    switch (t)
    {
    case 'A':
      out.push_back(0);
      break;
    case 'C':
      out.push_back(1);
      break;
    case 'G':
      out.push_back(2);
      break;
    case 'T':
      out.push_back(3);
      break;
    case 'N':
      out.push_back(4);
      break;
    }
  }
  return(out);
}

string int2string(const vector<int>& s)
{
  int i,t, len;
  len = s.size();
  
  string out;
  
  for (i=0; i<len; i++)
  {
    t = s[i];
    switch (t)
    {
    case 0:
      out.push_back('A');
      break;
    case 1:
      out.push_back('C');
      break;
    case 2:
      out.push_back('G');
      break;
    case 3:
      out.push_back('T');
      break;
    case 4:
      out.push_back('N');
      break;
    }
  }
  return(out);
}

void error(const string & error_message)
{
  cerr << "ERROR: " << error_message << endl;
  exit(1);
}

void warning(const string & warn_message)
{
  cerr << "WARNING: " << warn_message << endl;
}

double correlation(vector<double>& x, vector<double>& y)
{
  double NA = numeric_limits<double>::infinity();
  
  unsigned int length = x.size();
  if (length != y.size())
    error("correlation() x and y of different lengths!");
  
  double sx, sxx, sy, syy, sxy;
  sx = sxx = sy = syy = sxy = 0.0;
  double tx, ty;
  
  for (unsigned int i=0; i<length; i++)
  {
    tx = x[i];
    ty = y[i];
    
    if (tx == NA || ty == NA) continue;
    
    sx  += x[i];
    sy  += y[i];
  }
  double mx = sx/length;
  double my = sy/length;
  for (unsigned int i=0; i<length; i++)
  {    
    tx = x[i];
    ty = y[i];
    
    if (tx == NA || ty == NA) continue;
    
    sxx  += (tx-mx)*(tx-mx);
    syy  += (ty-my)*(ty-my);
    sxy  += (ty-my)*(tx-mx);
  }
  double r = sxy/sqrt(sxx*syy);
  
  return(r);
}


