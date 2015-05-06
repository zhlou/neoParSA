
#include "pwm.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <limits.h>
#include <boost/lexical_cast.hpp>
#define to_ boost::lexical_cast
#define to_string_ boost::lexical_cast<string>

PWM::PWM() 
{
  // default parameters
  pseudo     = 1.0;
  gc         = 0.5;
  maxscore   = 0;
  input_type = -1;
  
}

PWM::PWM(vector<vector<double> >& t, int type)
{
  // default parameters
  source     = string("");
  pseudo     = 1.0;
  gc         = 0.5;
  maxscore   = 0;
  input_type = type;

  setPWM(t, type);
}

void PWM::setPWM(vector<vector<double> >& t, int type, double g, double p)
{
  source     = string("");
  pseudo     = p;
  gc         = g;
  input_type = type;
  maxscore   = 0;
  
  setPWM(t, type);
}

void PWM::setPWM(vector<vector<double> >& e, int type)
{
  source     = string("");
  mat        = e;
  input_type = type;
  maxscore   = 0;
  
  // we need to verify the integrity of this pwm
  int pwmlen=mat.size();
  if (pwmlen <= 0)
    error("ERROR: setPWM() with matrix size less than 1");
  for (int i=0; i<pwmlen; i++)
  { 
    int n = mat[i].size();
    if (n == 4)
      mat[i].push_back(0);
    else if (n < 4 || n >5)
      error("ERROR: setPWM() does not have width of 4!");
  }
  switch (type)
  {
    case PCM:
      PCM2PFM();
      PFM2PSSM();
      break;
    case PFM:
      PFM2PSSM();
      break;
    case PSSM:
      break;
    case BEM:
      break;
    default:
      error("setPWM() unrecognized pwm type");
      break;
  }
  setNscore();
}

void PWM::setGC(double gc)
{
  if (input_type == -1)
    this->gc = gc;
  else
  {
    vector<vector<double> > replacement = getPWM(input_type);
    this->gc = gc;
    setPWM(replacement, input_type);
  }
}

void PWM::setPseudo(double pseudo)
{
  if (input_type == -1)
    this->pseudo = pseudo;
  else
  {
    vector<vector<double> > replacement = getPWM(input_type);
    this->pseudo = pseudo;
    setPWM(replacement, input_type);
  }
}

  
vector<vector<double> >& PWM::getPWM()
{
  return mat;
}

vector<vector<double> > PWM::getPWM(int type)
{
  // work with a copy of the matrix to start
  vector<vector<double> > out = mat;
  
  switch (type)
  {
    case PCM:
      PSSM2PFM(out);
      PFM2PCM(out);
      break;
    case PFM:
      PSSM2PFM(out);
      break;
    case PSSM:
      calc_max_score();
      break;
    case BEM:
      calc_max_score();
      break;
    default:
      error("getPWM() unrecognized pwm type " + to_string_(type));
      break;
  }
  return out;
}

void PWM::PCM2PFM()
{
  // we have counts. We need to add pseudo count and divide rows
  int pwmlen = mat.size();
  position_counts.resize(pwmlen);
  double gc_adjusted_count;
  
  for (int i=0; i<pwmlen; i++)
  {
    double rowsum = 0;
    for (int j=0; j<4; j++)
    {
      if (j==0 || j==3) // we have an A or T
        gc_adjusted_count = pseudo*(1-gc)/2;
      else             // we have a G or C
        gc_adjusted_count = pseudo*(gc)/2;

      rowsum       += mat[i][j];
      mat[i][j] += gc_adjusted_count;
    }
    position_counts[i] = rowsum;
    rowsum += pseudo;

    for (int j=0; j<4; j++)
      mat[i][j] /= rowsum;
  }
}

void PWM::PFM2PCM(vector<vector<double> >& t)
{
  unsigned int pwmlen = t.size();
  
  if (position_counts.size() != pwmlen)
  {
    warning("Position counts not set. Using default of 100.");
    position_counts.resize(pwmlen);
    for (unsigned int i=0; i<pwmlen; i++)
      position_counts[i] = 100;
  }
  
  double gc_adjusted_count;
  for (unsigned int i=0; i<pwmlen; i++)
  {
    double rowsum = position_counts[i] + pseudo;
    for (int j=0; j<4; j++)
    {
      t[i][j] *= rowsum;
      
      if (j==0 || j==3) // we have an A or T
        gc_adjusted_count = pseudo*(1-gc)/2;
      else             // we have a G or C
        gc_adjusted_count = pseudo*(gc)/2;

      t[i][j] -= gc_adjusted_count;
    }
  }
}

void PWM::PFM2PSSM()
{
  int pwmlen = mat.size();
  double bkgd;
  
  maxscore = 0;
  
  for (int i=0; i<pwmlen; i++)
  {
    double position_max = 0;
    for (int j=0; j<4; j++)
    {
      if (j==0 || j==3)   // we have an A or T
        bkgd = (1-gc)/2;
      else                // we have a G or C
        bkgd = (gc)/2;
      
      mat[i][j] = log( mat[i][j]/bkgd );
      position_max = max(position_max, mat[i][j]);
    }
    maxscore += position_max;
  }
}

void PWM::calc_max_score()
{
  int pwmlen = mat.size();
  maxscore = 0;
  for (int i=0; i<pwmlen; i++)
  {
    double position_max = 0;
    for (int j=0; j<4; j++)
      position_max = max(position_max, mat[i][j]);
    maxscore += position_max;
  }
}

/* this will check to make sure frequencies add to 1 */
void PWM::PSSM2PFM(vector<vector<double> >& t)
{
  int pwmlen = t.size();
  double bkgd;
  
  for (int i=0; i<pwmlen; i++)
  {
    double sum = 0;    
    for (int j=0; j<4; j++)
    {
      if (j==0 || j==3)   // we have an A or T
        bkgd = (1-gc)/2;
      else                // we have a G or C
        bkgd = (gc)/2;
      
      t[i][j] = bkgd * exp(t[i][j]);
      sum += t[i][j];
    }
    for (int j=0; j<4; j++)
      t[i][j] /= sum;
  }
}

void PWM::setNscore()
{
  int pwmlen = mat.size();
  for (int i=0; i<pwmlen; i++)
  {
    double minscore = mat[i][0];
    for (int j=1; j<4; j++)
      minscore = min(minscore, mat[i][j]);
    mat[i][4] = minscore;
  }
}

/* returns the threshold that will give this p value. It works... but I have 
no idea how. It was taken from MOODs in BioPerl */
double PWM::pval2score(double p)
{
  int n = mat.size();

  double dp_multi = 100.0;
  
  vector<double> bg;
  bg.push_back( (1-gc)/2 );
  bg.push_back( (  gc)/2 );
  bg.push_back( (  gc)/2 );
  bg.push_back( (1-gc)/2 );

  vector<vector<int> > tmat(4,vector<int>(n,0));

  int maxT = 0;
  int minV = INT_MAX;

  for (int i = 0; i < n; ++i)
  {
      for (int j = 0; j < 4; ++j)
      {
          if (mat[i][j] > 0.0){
              tmat[j][i] = (int) ( dp_multi * mat[i][j] + 0.5 );
          }
          else {
              tmat[j][i] = (int) ( dp_multi * mat[i][j] - 0.5 );
          }
      }
  }

  for (int i = 0; i < n; ++i)
  {
      int max = tmat[0][i];
      int min = max;
      for (int j = 1; j < 4; ++j)
      {
          int v = tmat[j][i];
          if (max < v)
              max = v;
          else if (min > v)
              min = v;
      }
      maxT += max;
      if (minV > min)
          minV = min;
  }

  int R = maxT - n * minV;
  
  vector<double> table0(R + 1, 0.0);
  vector<double> table1(R + 1, 0.0);

  for (int j = 0; j < 4; ++j)
      table0[tmat[j][0] - minV] += bg[j];

  for (int i = 1; i < n; ++i)
  {
      for (int j = 0; j < 4; ++j)
      {
          int s = tmat[j][i] - minV;
          for (int r = s; r <= R; ++r)
              table1[r] += bg[j] * table0[r - s];
      }
      for (int r = 0; r <= R; ++r)
      {
          table0[r] = table1[r];
          table1[r] = 0.0;
      }
  }

  double sum = 0.0;
  
  for (int r = R; r >= 0; --r)
  {
      sum += table0[r];
      if (sum > p)
      {
          return (double) ((r + n * minV + 1) / dp_multi);
      }
  }

	return (double) ((n * minV) / dp_multi);
}

/* a rather inefficient way of converting a score to a pvalue */
double PWM::score2pval(double s)
{
  if (s > maxscore)
    error("score greater than maxscore!");
  
  //calculates log10 p vals up to 10
  int precision = 1000; // number of steps to take
  if (s2p.size() == 0)
  {
    s2p.resize(precision);
    for (int i=1; i<precision; i++)
    {
      double mlog10p = (double) i/(precision/10.0); 
      double pval = pow(10.0,-mlog10p);
      s2p[i] = pval2score(pval);
    }
  }
      
  for (int i=(precision-1); i>0; i--)
  {
    if (s2p[i] < s)
    {
      double mlog10p = (double) i/(precision/10.0); 
      double pval = pow(10.0,-mlog10p);
      return(pval);
    }
  }
  return 1/std::pow(4.0, (int) mat.size()); // the max pvalue any matrix can have
}   

void PWM::subscore(const vector<int> & s, double * out)
{
  int i,j,k;
  out[0]=0.0;
  out[1]=0.0;
  int pwmlen = mat.size();
  for (i=0, j=(pwmlen-1); i<pwmlen; i++, j--)
  {
    vector<double>& row = mat[i]; 
    // forward sequence, i iterates forward
    k = s[i];
    out[0] += row[k];
    
    // reverse sequence, j iterates back over subseq, 3-n give complement
    if (s[j] !=4) 
      k = 3-s[j];
    else 
      k = 4;
    out[1] += row[k];
  }
  out[2] = max(out[0],out[1]);
}

void PWM::score(const vector<int>& s, TFscore &t)
{
  /* I do something a little unorthodox here. I want to control for boundary
  effects, so I pad the beginning and end of the sequence with Ns and score
  using the boundary behavior. Unlike the transc code, I report the score for the
  middle of the binding site rather than the m position. This means my output is
  of the same length as the sequence and it is the same length for every factor */
  
  int i, j, k;
  int mdist, ndist;
  int start, end;
  int len;
  int pwmlen = mat.size();
  vector<int> sub(pwmlen);
  
  ndist = pwmlen/2; // returns floor for middle position
  mdist = pwmlen-ndist;
  
  double * out = new double[3];
  len = s.size();

  t.fscore.resize(len);
  t.rscore.resize(len);
  t.mscore.resize(len);
  //t.tfname = tfname;
  
  for (i=0; i<len; i++)
  {
    start = i - mdist;
    end   = i + ndist;
    // j = index in s
    // j = index in pwm
    for (j=start, k=0; j<end; j++, k++)
    {
      if (j<0 || j>=len)
      {
        sub[k]=4;
      } else 
      {
        sub[k]=s[j];
      }
    }
    subscore(sub,out);
    t.fscore[i] = out[0];
    t.rscore[i] = out[1];
    t.mscore[i] = out[2];
  } 
  delete[] out;
}

  


void PWM::print(ostream& os, int precision)
{
  int p = precision;
  int w = p + 7;
  int pwmlen = mat.size();
  vector<char> bases;
  bases.push_back('A');
  bases.push_back('C');
  bases.push_back('G');
  bases.push_back('T');

  for (int i=0; i<4; i++)
  {
    os << setw(w) << setprecision(p) << bases[i];
    for (int j=0; j<pwmlen; j++)
      os << setw(w) << setprecision(p) << mat[j][i];
    os << endl;
  }
}
