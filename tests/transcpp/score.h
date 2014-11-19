/*********************************************************************************
*                                                                                *
*     score.h                                                                    *
*                                                                                *
*     Presumably we could use any number of functions to score the model.        *
*     We have some set of predictions and outputs, expressed here as pointer     *
*     vectors, as well as some rule for weighting the scores of each             *
*                                                                                *
*********************************************************************************/

#ifndef SCORE_H
#define SCORE_H

#include "organism.h"
#include "gene.h"
#include "mode.h"
#include "scalefactor.h"

// forward declare organism so that we can get info from it
class Organism;

/*****************  Scoring Functions  ******************************************/

/* Each scoring function takes two vectors of double pointers and outputs a single
double representing the score of that function */

//double sse(vector<double*>& data, vector<double*>& prediction, scale_factor_ptr);

class Score
{
private:
  mode_ptr  mode;
  genes_ptr genes;
  
  vector<int> ids;
  
  vector<vector <double*> > data;
  vector<vector <double*> > prediction;

  vector<scale_factor_ptr> scale;
  
  vector<double> scores;
  vector<double> weights;
  double min_weight;
  
  double score;
  
  // weighting functions
  typedef void (Score::*WFP)();
  vector<WFP> weight_functs;
  void area();
  void height();
  void height2();
  void gene_number();
  void nuc_number();
  
  // scoring functions
  typedef void (Score::*SFP)();
  SFP scoreFunc;
  void sse();
  void chisq();
  void pdiff();
  void rms();
  void arkim();
  void sum_slope_squares();
  
  
public:
  
  // Constructors
  Score();
  Score(Organism* parent);
  
  // Getters
  double getScore();
  
  // Setters
  void set(Organism* parent);
  void setWeights();
  
  // Method
  void checkScale(ostream& os);
  
  // I/O
  void print(ostream& os);
};
  
 
  




#endif
  
  
