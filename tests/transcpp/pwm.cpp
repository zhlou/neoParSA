/*********************************************************************************
*                                                                                *
*     pwm.cpp                                                                    *
*                                                                                *
*     Contains objects and functions for use in creating pwms, log odd pwms      *
*                                                                                *
*********************************************************************************/

# include "pwm.h"


/*********************************   PwmContainer    ********************************/


/*    Constructor   */

PwmContainer::PwmContainer() {npwms = 0;}


/*    Getters   */

const Pwm& PwmContainer::getPwm(string& n)
{
  int i;
  for (i=0; i<npwms; i++)
  {
    if (!pwms[i].getName().compare(n))
      return(pwms[i]);
  }
  cerr << "ERROR: getPWM() PWM with name " << n << " not found" << endl;
  exit(1);
}

int PwmContainer::length() const { return(npwms);}


/*    Setters   */

void PwmContainer::addPwm(Pwm p) { pwms.push_back(p); }

void PwmContainer::addPwm() 
{ 
  Pwm p;
  pwms.push_back(p); 
}

void PwmContainer::readPwms(string& fname)
{
  string line, pname;
  vector< vector<int> > counts;
  vector<int> cline;
  string sub;
  int first = 1;
  int temp;
  
  ifstream file (fname.c_str());
  if (file.is_open())
  {
    while (getline(file,line))
    {
      if (!line.empty())
      {
        if (line.at(0)=='>') // new pwm
        {
          npwms++;
          if (!first) 
            pwms.push_back(Pwm(counts,pname));
          pname = line.substr(1,line.length()-1);
          counts.clear();
          first = 0;
        } else {
          cline.clear();
          istringstream iss(line);
          while (iss) // iss splits into whitespace fields
          {
            iss >> sub;
            temp = atoi(sub.c_str());
            cline.push_back(temp);
          }
          counts.push_back(cline);
        }
      }   
    }
    pwms.push_back(Pwm(counts,pname));
  } else {
    cerr << "ERROR: readPWMs() could not find file with name " << fname << endl;
    exit(1);
  }
  file.close();
}
  

/*********************************   Pwm    ************************************/


/*    Constructors    */

Pwm::Pwm() {}

Pwm::Pwm(int w, string n)
{
  int i,j;
  
  counts.resize(w);
  frequencies.resize(w);
  width = w;
  pname.assign(n);
  
  for (i=0; i<width; i++)
  {
    for (j=0; j<4; j++)
    {
      counts[i].resize(4);
      frequencies[i].resize(4);
      
      counts[i][j]      =   25;
      frequencies[i][j] = 0.25;
    }
  }
}

Pwm::Pwm(vector<vector<int> > c, string n)
{
  int i, j;
  double sum;
  
  counts = c;
  width = counts.size();
  frequencies.resize(width);
  pname.assign(n);
  
  for (i = 0; i<width; i++)
  {
    frequencies[i].resize(4);
    sum = 0;
    for (j=0; j<4; j++)
      sum += counts[i][j];
    for (j=0; j<4; j++)
      frequencies[i][j] = counts[i][j]/sum;
  }
}


/*    Getters   */

string Pwm::getName() const { return(pname); }


/*    Methods   */

vector< vector<double> > Pwm::getLogOddPwm(int pseudo_count, double gc) const
{
  int i,j;
  double a,c,g,t, f; 
  vector< vector<double> > out(width);
  vector<int> nalign(width);
  c = g = gc/2;
  a = t = (1-gc)/2;
  
  for (i=0; i<width; i++)
  {
    out[i].resize(4);
    nalign[i] = 0; 
    for (j=0; j<4; j++)
      nalign[i] += counts[i][j];
    nalign[i] += pseudo_count;
    
    f = (counts[i][0]+a*pseudo_count)/nalign[i];
    out[i][0] = log(f/a);
    
    f = (counts[i][1]+c*pseudo_count)/nalign[i];
    out[i][1] = log(f/c);
    
    f = (counts[i][2]+g*pseudo_count)/nalign[i];
    out[i][2] = log(f/g);
    
    f = (counts[i][3]+t*pseudo_count)/nalign[i];
    out[i][3] = log(f/t);
  }
  return(out);
}
  

/*    Print   */

void Pwm::printCounts(ostream& os) const
{
  int i,j;
  int w = 8;
  os << ">" << pname << "_pcm" << endl;
  os << setw(w) << "A"
     << setw(w) << "C"
     << setw(w) << "G"
     << setw(w) << "T"
     << endl;
     
  for (i=0; i<width; i++)
  {
    for (j=0; j<4; j++)
    {
      os << setw(w) << counts[i][j];
    }
    os << endl;
  }
}

void Pwm::printFrequencies(ostream& os) const
{
  int i,j;
  int w = 8;
  os << ">" << pname << "_pfm" << endl;
  os << setw(w) << setprecision(2) << "A"
     << setw(w) << setprecision(2) << "C"
     << setw(w) << setprecision(2) << "G"
     << setw(w) << setprecision(2) << "T"
     << endl;
     
  for (i=0; i<width; i++)
  {
    for (j=0; j<4; j++)
    {
      os << setw(w) << frequencies[i][j];
    }
    os << endl;
  }
}

  
