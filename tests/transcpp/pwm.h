/*********************************************************************************
*                                                                                *
*     pwm.h                                                                      *
*                                                                                *
*     Contains objects and functions for use in creating pwms, log odd pwms      *
*     and using them to score a sequence                                         *
*                                                                                *
*********************************************************************************/

using namespace std;

# include <iostream>
# include <fstream>
# include <vector>
# include <iomanip>
# include <string>
# include <sstream>
# include <cstdlib>
# include <math.h>

class Pwm 
{
private:
  string pname;
  int width;
  // counts and frequencies as x[position][base]
  vector< vector<int> > counts;
  vector< vector<double> > frequencies;
public:
  Pwm();
  //Pwms are constructed from count data
  Pwm(vector< vector<int> > c, string n);
  // a default pwm can be specified with w, giving each prob of 0.25
  Pwm(int w, string n);
  //log odd pwm can be derived with pseudo count and gc content  
  vector< vector<double> > getLogOddPwm(int pseudo_count, double gc) const;
  string getName(void) const;
  void printCounts(ostream& os) const;
  void printFrequencies(ostream& os) const;
};

class PwmContainer // bit of a misnomer since it is a vector... 
{
private:
  vector<Pwm> pwms;
  int npwms;
public:
  // Constructor
  PwmContainer();
  
  // Getters
  int length() const;
  const Pwm& getPwm(string& n);
  
  // Setters
  void addPwm(Pwm p);
  void addPwm();
  void readPwms(string& fname);
};
