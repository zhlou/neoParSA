/*********************************************************************************
*                                                                                *
*     nucleus.h                                                                  *
*                                                                                *
*     Constains the nucleus class, which holds references to each member         *
*                                                                                *
*********************************************************************************/

#ifndef NUCLEUS_H
#define NUCLEUS_H 

#include "nuclei.h"
#include "gene.h"
#include "distance.h"
#include "TF.h"
#include "conc.h"
#include "bindingsite.h"
#include "subgroup.h"
#include "promoter.h"
#include "bindings.h"
#include "quenching.h"
#include "flags.h"

#include <map>
#include <boost/utility.hpp>

using namespace std;



class Nuclei;

class Nucleus : boost::noncopyable
{
private:
  int id;
  int idx;
  genes_ptr     genes;            // references to the genes in this nucleus
  tfs_ptr       tfs;              // references to the TFs in this nucleus
  bindings_ptr  bindings;         // all the binding site data
  conc_ptr      tfdata;           // the concentration of each TF in this nucleus
  //conc_ptr      ratedata;
  promoters_ptr promoters;
  subgroups_ptr subgroups;
  quenching_ptr quenching;
  
  map<TF*, double> conc;
  map<Gene*, double> Ns; // total activation input
  map<Gene*, double> Rs; // dmRNA/dt
 
  void calcN();
  void calcN(Gene& gene);
  void calcR();
  void calcR(Gene& gene);
  
public:
  // Constructors
  Nucleus();
  Nucleus(int id, int idx, Nuclei* parent);
  
  // Getters
  double getRate(Gene&);
  double getN(Gene&);
  int    getID() {return id;}
  
  // Methods
  bool hasSameTFs(tfs_ptr t);
  
  void updateKV(TF&);
  void updateKV();
     
  void calcQuenching();
  void calc_f();
   
  void updateN();
  void updateR();
  
  // Print
  void printSites(ostream& os);
};
  
  

typedef boost::shared_ptr<Nucleus> nucleus_ptr;

  
  
  
  
  
  
  
#endif
  
  
  

