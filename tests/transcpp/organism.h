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
#include "conc.h"
#include "score.h"

#include "nuclei.h"
#include "promoter.h"
#include "parameter.h"
#include "quenching.h"
#include "bindings.h"
#include "scalefactor.h"
#include "coeffects.h"

using namespace std;

class   Nuclei;
class   Score;
typedef boost::shared_ptr<Score>  score_ptr;
typedef boost::shared_ptr<Nuclei> nuclei_ptr;

class Organism
{
private:
  
  mode_ptr       mode;
  distances_ptr  distances;
  tfs_ptr        master_tfs;      // may vary between nuclei
  genes_ptr      master_genes;    // may vary between nuclei
  conc_ptr       tfdata;
  conc_ptr       ratedata;
  promoters_ptr  promoters;
  coops_ptr      coops;
  coeffects_ptr  coeffects;
  score_ptr      score_class;
  
  scale_factors_ptr scale_factors;
  
  param_ptr_vector params; // typedef in parameter.h
     
  vector<nuclei_ptr> nuclei;
  
  ptree annealer_input;
  ptree move;
  ptree count_criterion;
  ptree mix;
  ptree lam;
  void readAnnealing(ptree& pt);
  
  // we need to translate what is in nuclei objects to an array that corresponds
  // to the embryo itself. When we read in we store the order we read 
  vector<int> ids;
  
  double score_out;
  
  // move functions
  typedef void (Organism::*MFP)(int);
  map<string, MFP> move_map;
  map<string, MFP> restore_map;
  
  void setMoves();
  void checkParameters();
  
  void moveThreshold(int);
  void restoreThreshold(int);
  
  void moveKmax(int);
  void restoreKmax(int);

  void moveLambda(int);
  void restoreLambda(int);
  
  void moveCoef(int);
  void restoreCoef(int);
  
  void moveQuenchingCoef(int);
  void restoreQuenchingCoef(int);
  
  void moveQuenching(int);
  void restoreQuenching(int);
  
  void moveCoeffect(int);
  void restoreCoeffect(int);
  
  void moveScaleFactor(int);
  void moveQ(int);
    
public:
  // Constructors
  Organism();
  Organism(ptree &pt);
  
  void resetAll();
  
  // Getters
  conc_ptr       getTFData()         {return tfdata;         }
  conc_ptr       getRateData()       {return ratedata;       }
  promoters_ptr  getPromoters()      {return promoters;      } 
  distances_ptr  getDistances()      {return distances;      }
  genes_ptr      getGenes()          {return master_genes;   }
  tfs_ptr        getTFs()            {return master_tfs;     }
  mode_ptr       getMode()           {return mode;           }
  scale_factors_ptr getScales()      {return scale_factors;  }
  vector<int>    getIDs()            {return ids;            }
  double*        getPrediction(Gene&,int);
  double*        getData(Gene&,int,bool);
  
  // Setters
  void populate_nuclei();
  void populate_nuclei_multiple();
  void populate_nuclei_single();
  
  // Methods
  void score();
  void calc_f();
  void scramble();
  void checkScale(ostream&);

  
  void moveKcoop(int);
  void restoreKcoop(int);
  void setParam(int idx, double val) {params[idx]->set(val);}
  
  // I/O
  void write(string, ptree& pt);
  void printParameters(ostream& os);
  
  void printSites(ostream& os);                     
  void printSites(Gene& gene, ostream& os);         
  void printSites(TF& tf, ostream& os);             
  void printSites(Gene& gene, TF& tf, ostream& os); 
  
  void printSubgroups(Gene& gene, ostream& os);
  void printOccupancy(Gene& gene, ostream& os);
  void printModeOccupancy(Gene& gene, ostream& os);
  void printEffectiveOccupancy(Gene& gene, ostream& os);
  
  void printRate(ostream& os);
  void printRateData(ostream& os);
  void printScore(ostream& os);
  
  
  // Annealing
  int    getDimension() const; 
  double get_score();
	void   generateMove(int idx, double theta);
	void   restoreMove(int idx);
	void   serialize(void *buf) const;
	void   deserialize(void const *buf);
	int    getStateSize();
};
  
  


































#endif
