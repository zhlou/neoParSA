#ifndef PWM_H
#define PWM_H

#include "utils.h"
#include <vector>

using namespace std;

/* there are several ways we specify pwm information:

PCM, or Position Count Matrix, contains the counts of 
observed bases at each position

PFM, or Position Frequency Matrix, contains the
probablity of a base at each position

PSSM, the Position Specific Scoring Matrix, contains
the log-likelihood of each base at each position,
taking into account the background model

BEM, the Binding Energy Model, contains what can be
interpreted as the ddG for a substitution of each
base at each position, in the case where lambda=kT
THIS MODE HAS BEEN REMOVED! The conversion of PSSM to BEM
is lossy and I decided there was no real benefit to using
BEM over PSSM. BEMs can still be entered, but they are scored
as PSSMs with a maxscore of 0

*/

enum pwm_types { PCM = 0, PFM = 1, PSSM = 2, BEM = 3};

/* struct holding the score of a TF over sequence */
struct TFscore 
{
  //TF* tf;
  vector<double> fscore;
  vector<double> rscore;
  vector<double> mscore;
  double maxscore;
};

class PWM
{
private:
  string source;
  
  int input_type;
  // the underlying matrix
  vector<vector<double> > mat;
  
  /* the observed total counts at each position,
  necessary for backwards converstion */
  vector<int>    position_counts;
  /* scores the pvalue at a given score */
  vector<double> s2p;

  double pseudo;   // pseudo count used to convert from PCM
  double gc;       // the GC content of sequences used to generate PWM
                   // Note: this is an experimental parameter for generating a PWM, not a general property of drosophila
  double maxscore; // the maximum score a PSSM can give

  // private methods
  void subscore(const vector<int> & s, double * out);
  
public:
  // constructors
  PWM();
  PWM(vector<vector<double> >& t, int type);
  
  // setters
  void setSource(string source) { this->source = source; }
  void setGC(double gc);         
  void setPseudo(double pseudo); 
  
  void setPWM(vector<vector<double> >& t, int type); // use default gc=0.25, pseudo=1
  void setPWM(vector<vector<double> >& t, int type, double gc, double pseudo);

  // getters
  const string& getSource()    { return source; }
  double getGC()        { return gc; }
  double getPseudo()    { return pseudo; }
  double getMaxScore()  { return maxscore; }
  int    length()       { return mat.size(); }
  int    getInputType() { return input_type; }
  vector<vector<double> >& getPWM(); // efficient, return reference to pwm
  vector<vector<double> >  getPWM(int type); // return a copy, for conveneince, not inner loop

  // forward converstions
  void PCM2PFM();
  void PFM2PSSM();

  // reverse conversions (dont touch the actual matrix since PSSMs are always used internally)
  void PFM2PCM(vector<vector<double> >& t);
  void PSSM2PFM(vector<vector<double> >& t);

  // methods
  void   setNscore();
  void   calc_max_score();
  double pval2score(double pval);  // returns the threshold that would yeild a given p-value
  double score2pval(double score); // returns the pvalue of a given score
  void   score(const vector<int>& s, TFscore &t);
  
  void print(ostream& os, int precision);
};

#endif

