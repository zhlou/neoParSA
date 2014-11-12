/*********************************************************************************
*                                                                                *
*     occupancy.h                                                                *
*                                                                                *
*     It has become increasingly clear over the course of writing this code      *
*     that most of the intense calculations are much more efficient if they      *
*     can simply traverse an array of occupancies. This class is designed        *
*     to be a matrix of occupancies and effective occupancies and all subgroup   *
*     and quenching calculations will take a pointer to this object to do        *
*     their calculations                                                         *
*                                                                                *
*********************************************************************************/

#ifndef OCCUPANCY_H
#define OCCUPANCY_H

#include "gene.h"
#include "TF.h"
#include "bindings.h"
#include "conc.h"

#include <boost/shared_ptr.hpp>

using namespace std;

//class Bindings;
//typedef boost::shared_ptr<Bindings> bindings_ptr;

struct OccData
{
  vector< vector<double> > kv;
  vector< vector<double> > occupancy;
  vector< vector<double> > effective_occupancy;
  
  vector< vector<double> > saved_kv;
  vector< vector<double> > saved_occupancy;
  vector< vector<double> > saved_effective_occupancy;
};



// occupancy will hold only for a single gene
class Occupancy
{
private:
  tfs_ptr       tfs;
  genes_ptr     genes;
  bindings_ptr  bindings;
  conc_ptr      tf_data;
  
  int nnuc;
  
  // this is a one-to-one match to the bindings array
  // occ_data[gene][tf].x[site][nuc]
  map<Gene*, map<TF*, OccData> > occ_data;
  
  
public:
  // Constructors
  Occupancy();
  
  // Setters
  void setGenes(genes_ptr x)       { genes    = x; }
  void setTFs(tfs_ptr x)           { tfs      = x; }
  void setTFData(conc_ptr x)       { tf_data  = x; }
  void setBindings(bindings_ptr x) { bindings = x; }
  void setNnuc(int x)              { nnuc     = x; }
  
  void create();
  
  // Getters
  int getNnuc() { return nnuc; }
  
  // Methods
  void update(Gene& gene, TF& tf);
  void calcKV(Gene& gene, TF& tf);
  
  void update(TF& tf);
  void updateKV(TF& tf);
  
  void saveOccupancy();
  void restoreOccupancy();
  
  void saveEffectiveOccupancy();
  void restoreEffectiveOccupancy();
  
  void saveAllOccupancy();
  void restoreAllOccupancy();
  
  // I/O
  void printOccupancy(Gene& gene, ostream& os);
  void printEffectiveOccupancy(Gene& gene, ostream& os);
};

typedef boost::shared_ptr<Occupancy> occupancy_ptr;



#endif

