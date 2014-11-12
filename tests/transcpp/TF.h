/*********************************************************************************
*                                                                                *
*     TF.h                                                                       *
*                                                                                *
*     Contains class definition for Transcription Factors                        *
*                                                                                *
*********************************************************************************/

#ifndef TF_H
#define TF_H

#include "parameter.h"
#include "cooperativity.h"
#include "coeffects.h"
#include <vector>
#include <string>
#include <cstdlib>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>



using namespace std;
using boost::property_tree::ptree;


/* struct holding the score of a TF over sequence */
struct TFscore 
{
  //TF* tf;
  vector<double> fscore;
  vector<double> rscore;
  vector<double> mscore;
  double maxscore;
};

class TF : boost::noncopyable
{
private: 
  vector< vector<double> > energy; // binding energy model
  string tfname;
  int       pwmlen;
  int       bsize;
  double    gc;        // if the pwm was obtained from drosophila sequence it has a bias
                       // against gc. This corrects for said bias.  
                       
  param_ptr kmax;      // true binding affinity confounded with relative protein levels
  param_ptr threshold;
  param_ptr lambda;
  
  vector<param_ptr> coefs;
  
  double    maxscore;  // maximum score a TF can give
  int       Nbehavior; // how to treat Ns in sequence 1: conservative (min) 2: median 3: aggressive (max)

  vector< pair<TF*, coop_ptr> >     coops;
  vector< pair<TF*, coeffect_ptr> > coeffects;
  
  void subscore(const vector<int> & s, double* out); // for scoring sequence of pwmlen
  void readPWM(ptree& pt);
  void PCM2PFM(double pseudo);
  void PFM2PSSM();
  void PSSM2BEM();
  void setNscore();
  
public:
  // constructors
  TF();           // empty constructor
  TF(ptree& pt);  // create from property tree
  
  // getters
  const string&  getName() const;
  double         getThreshold();
  int            getBindingSize() const;
  double         getKmax();
  double&        getCoef();
  double         getModifiedCoef();
  double         getLambda();
  double         getMaxScore() const;
  bool           hasChanged() const;
  void           getParameters(param_ptr_vector& p); // pushes parameters onto argument
  bool           checkCoops(TF*, char, char); // check if this cooperates with tf
  coop_ptr       getCoop(TF*);
  bool           neverActivates();
  bool           neverQuenches();
  vector<double> getCoefs();
  int            getNumModes() { return coefs.size(); }
  
  vector< pair<TF*, coeffect_ptr> >& getTargets() { return coeffects; }   
  
  // setters
  void set(ptree& pt); // set entire TF
  void setName(string n);
  void setPwm(vector< vector<double> > p);
  void setKmax(double k);
  void setLambda(double l);
  void setThreshold(double t);
  void setBindingSize(int b);
  void setCoops(vector< pair<TF*,coop_ptr> > p) { coops = p; }
  void setCoeffects(vector< pair<TF*, coeffect_ptr> > p) { coeffects = p; }
  
  // methods
  TFscore score(const string & s);     // score a string with tf
  TFscore score(const vector<int>& s); // score an int vector (faster) 
  void    score(const vector<int>& s, TFscore &t); // pass by reference (fastest)
  
  // I/O
  void print(ostream& os);
  void write(ostream& os); // write XML to output
  void write(ptree& tfsnode); // write XML to property tree

};

typedef boost::shared_ptr<TF> tf_ptr;
typedef vector<tf_ptr> tf_ptr_vector;

class TFContainer
{
private:
  tf_ptr_vector tfs;

public:
  // Constructor
  TFContainer();
  TFContainer(ptree& pt);
  
  // Getters
  TF& getTF(const string& n);
  TF& getTF(int index);
  tf_ptr getTFptr(int index);
  tf_ptr getTFptr(const string& n);
  int size() const;
  void getParameters(param_ptr_vector& p);
  
  // Setters
  void add(ptree& pt);
  void add(tf_ptr t);
  void setCoops(coops_ptr);
  void setCoeffects(coeffects_ptr);
  
  // Output
  void write(ostream& os) const;
  void write(ptree& pt) const;
};

typedef boost::shared_ptr<TFContainer> tfs_ptr;

#endif
