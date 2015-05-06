/*********************************************************************************
*                                                                                *
*     organism.h                                                                 *
*                                                                                *
*     An organism is defined here as a collection of nuclei. This is the master  *
*     class which holds most of the data. It creates and array of nuclei         *
*     dynamicaly based on the data present to be scored against.                 *
*     Its main purpose is to be an interface for the model and the annealer.     *
*                                                                                *
*********************************************************************************/

#ifndef ORGANISM_H
#define ORGANISM_H

#include "mode.h"
#include "distance.h"
#include "gene.h"
#include "TF.h"
#include "datatable.h"
#include "score.h"

#include "nuclei.h"
#include "promoter.h"
#include "parameter.h"
#include "quenching.h"
#include "bindings.h"
#include "scalefactor.h"
#include "coeffects.h"
#include "competition.h"

using namespace std;

class   Nuclei;
class   Score;
typedef boost::shared_ptr<Score>  score_ptr;
typedef boost::shared_ptr<Nuclei> nuclei_ptr;

class Organism
{
private:
 
  mode_ptr          mode;
  distances_ptr     distances;
  tfs_ptr           master_tfs;      // may vary between nuclei
  genes_ptr         master_genes;    // may vary between nuclei
  table_ptr         tfdata;
  table_ptr         ratedata;
  promoters_ptr     promoters;
  scale_factors_ptr scale_factors;
  coeffects_ptr     coeffects;
  score_ptr         score_class;
  coops_ptr         coops;
  competition_ptr   competition;
  
  param_ptr_vector params; 
  param_ptr_vector all_params;
  
  nuclei_ptr nuclei;
  
  // we need to translate what is in nuclei objects to an array that corresponds
  // to the embryo itself. When we read in we store the order we read 
  vector<string> ids;
  
  int move_count;
  double score_out;
  double previous_score_out;
  
  /*  Move Functions  */
  void setMoves();
  typedef void (Organism::*MFP)(int);
  // just map input string to a move function, used for initialization
  map<string, MFP> move_map;
  map<string, MFP> restore_map;
  
  // the move function that gets called during annealing
  vector<MFP> moves;
  vector<MFP> restores;
  
  // the move functions
  // void ResetAll(int); 
  void moveScores(int);
  void movePWM(int);
  void moveSites(int);
  void moveLambda(int);
  void moveKmax(int);
  void moveCoopD(int);
  void moveKcoop(int);
  void moveCoef(int);
  void moveQuenching(int);
  void moveQuenchingCoef(int);
  void moveCoeffect(int);
  void moveCoeffectEff(int);
  void movePromoter(int);
  void null_function(int);  
  
  // the restore functions
  void restoreScores(int);
  void restorePWM(int);
  void restoreSites(int);
  void restoreLambda(int);
  void restoreKmax(int);
  void restoreCoopD(int);
  void restoreKcoop(int);
  void restoreCoef(int);
  void restoreQuenching(int);
  void restoreQuenchingCoef(int);
  void restoreCoeffect(int); 
  void restoreCoeffectEff(int);   

public:
  // Constructors
  Organism();
  Organism(ptree &pt, mode_ptr);
  Organism(string fname, string section);
  
  void initialize(ptree& pt);
  
  // Getters
  table_ptr         getTFData()         {return tfdata;         }
  table_ptr         getRateData()       {return ratedata;       }
  promoters_ptr     getPromoters()      {return promoters;      } 
  distances_ptr     getDistances()      {return distances;      }
  genes_ptr         getGenes()          {return master_genes;   }
  tfs_ptr           getTFs()            {return master_tfs;     }
  mode_ptr          getMode()           {return mode;           }
  scale_factors_ptr getScales()         {return scale_factors;  }
  nuclei_ptr        getNuclei()         {return nuclei;         }
  competition_ptr   getCompetition()    {return competition;    }
  vector<string>    getIDs()            {return ids;            }
  double*           getPrediction(Gene&,string&);
  double*           getData(Gene&,string&);
  int               getNNuc()           {return ratedata->getNames("ID").size();}
  int               getNGenes()         {return master_genes->size();}
  double            getTotalScore()     {return get_score();}
  vector<double>&   getN(int gidx);
  bindings_ptr      getBindings();
  
  // Setters
  void populate_nuclei();
  void setMode(mode_ptr mode) {this->mode = mode;}
  void addTF(tf_ptr tf) { master_tfs->add(tf); Recalculate(); }
  
  // Methods
  void Recalculate(); 
  void ResetAll(int); 
  void score();
  void calc_f();
  void scramble();
  void permute(string& table, string& by);
  void checkScale(ostream&);
  void move(int idx);

  iparam_ptr getParam(int idx) {return params[idx];}
  string getParamName(int idx);

  template < typename T >
  int  getParamValue(int idx);
  
  // Matlab
  int test_int;
  int test() {test_int++; return test_int;}
  
  // I/O
  void write(string, ptree& pt);
  void printParameters(ostream& os);
  
  void printSites(ostream& os);                     
  void printSites(Gene& gene, ostream& os);         
  void printSites(TF& tf, ostream& os);             
  void printSites(Gene& gene, TF& tf, ostream& os); 
  
  void printScores(Gene& gene, ostream& os);         
  
  void printSubgroups(Gene& gene, ostream& os);
  void printOccupancy(Gene& gene, ostream& os, bool invert);
  void printModeOccupancy(Gene& gene, ostream& os, bool invert);
  void printEffectiveOccupancy(Gene& gene, ostream& os, bool invert);
  
  void printRate(ostream& os, bool invert);
  void printRateData(ostream& os, bool invert);
  void printScore(ostream& os);
  
  void printR2D(ostream& os);
  void printN2D(ostream& os);
  
  
  // Annealing
  int    getDimension() const; 
  double get_score();
  void   generateMove(int idx, double theta);
  void   restoreMove(int idx);
  void   serialize(void *buf) const;
  void   deserialize(void const *buf);
  int    getStateSize();
};
  

typedef boost::shared_ptr<Organism> organism_ptr;





























#endif
