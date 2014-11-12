/*********************************************************************************
*                                                                                *
*     bindings.h                                                                 *
*                                                                                *
*     Constains a map of all site on all genes                                   *
*                                                                                *
*********************************************************************************/

#ifndef BINDINGS_H
#define BINDINGS_H

#include "TF.h"
#include "gene.h"
#include "bindingsite.h"
#include "conc.h"

#include <boost/shared_ptr.hpp>

using namespace std;

/* Binding sites were initialized as poitners because it makes it simple to
ensure that subgroups and quenching interactions never point to null. This
may or may not be what we actually want and should be considered in the
future */


typedef map<Gene*, map<TF*, site_ptr_vector> > site_map;
typedef map<Gene*, map<TF*, TFscore> >         scores_map;


class Bindings
{
private:
  int nnuc;
  
  genes_ptr genes;
  tfs_ptr   tfs;
  conc_ptr  tfdata;
  
  map<int, int> id_2_idx;
  map<int, int> idx_2_id;
  
  
  //scores[gene][tf] = TFscore
  scores_map scores;
  scores_map saved_scores;
  
  // sites[gene][tf] -> BindingSite
  site_map sites;                    
  site_map saved_sites;
  
  map<TF*, vector<double> > conc;
  
  bool hasScores(Gene&, TF&);
  bool hasSites(Gene&, TF&);
  void createSite(site_ptr_vector& tmp_sites, Gene& gene, TF& tf,
                  int pos, double bsize, double score, char orientation, 
                  double lambda, double kmax, double maxscore,
                  vector<double>& v, int nmodes);
  
public:
  
  // Constructors
  Bindings();
  
  // Setters
  void setGenes(genes_ptr g) {genes = g;  }  
  void setTFs(tfs_ptr t)     {tfs = t;    }
  void setTFData(conc_ptr c) {tfdata = c; }
  
  void create();
  void clear();
  
  void createScores();
  void createSites();
  
  void createScores(Gene&, TF&);  // score a gene for a tf 
  void createSites(Gene&, TF&);          
  
  void addSite(Gene* g, TF* t, site_ptr b) { sites[g][t].push_back(b); }
  
  // Functions for moving
  void saveScores(TF&);
  void updateScores(TF&);
  void updateScores();
  void restoreScores(TF&);
  
  void saveSites(TF&);
  void updateSites(TF&);
  void updateSites();
  void restoreSites(TF&);
  
  void updateK(TF&);
  void updateKandLambda(TF&);
  
  void saveAllOccupancy();
  void restoreAllOccupancy();
  
  void saveOccupancy();
  void restoreOccupancy();
  
  void saveModeOccupancy();
  void restoreModeOccupancy();
  
  void saveEffectiveOccupancy();
  void restoreEffectiveOccupancy();
  
  site_ptr_vector& getSites(Gene&, TF&);
  map<TF*, site_ptr_vector>& getSites(Gene&);
  
  void addNuc(int nuc_id);
  int  getNnuc() {return nnuc; }
  
  void printSites(ostream& os);
  void printSites(Gene& gene, ostream& os);
  void printSites(TF& tf,     ostream& os);
  void printSites(Gene& gene, TF& tf, ostream& os);
  
  void printScores(ostream& os);
  
  void printTotalOccupancy(Gene& gene, ostream& os);
  void printModeOccupancy(Gene& gene, ostream& os);
  void printEffectiveOccupancy(Gene& gene, ostream& os);
  
  
};




typedef boost::shared_ptr<Bindings>            bindings_ptr;








































#endif


