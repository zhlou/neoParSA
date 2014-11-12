/*********************************************************************************
*                                                                                *
*     nuclei.h                                                                   *
*                                                                                *
*     A nuclei object contains all nuclei which share the same TFs. This means   *
*     they share subgroups and quenching interactions.                           *
*                                                                                *
*     In the original code, many of the calculations on occupancy and quenching  *
*     were multiplications by 0. This class was created to prevent such needless *
*     calculations.                                                              *
*                                                                                *
*********************************************************************************/

#ifndef NUCLEI_H
#define NUCLEI_H 

#include "organism.h"

#include "gene.h"
#include "conc.h"
#include "promoter.h"
#include "subgroup.h"
#include "quenching.h"
#include "TF.h"
#include "mode.h"

/* For the most part, nuclei does not own the private data inside it. It simply
points to the data from it's parent class (Organism). The notable exceptions
are the QuenchingInteractions and Subgroups */

class Organism;

class Nuclei
{
private:
  int           n;
  
  // owned by parent object (organism)
  genes_ptr     genes;
  tfs_ptr       tfs;
  conc_ptr      tfdata;
  distances_ptr distances;
  promoters_ptr promoters;
  mode_ptr      mode;
  
  // owned by Nuclei
  bindings_ptr  bindings;
  subgroups_ptr subgroups;
  quenching_ptr quenching;
  modifying_ptr coeffects;
  
  vector<int> IDs;

  map<Gene*, vector<double> > Ns;
  map<Gene*, vector<double> > Rs;
  
public:
  // Constructors
  Nuclei(Organism* parent, tfs_ptr);
  
  //  Setters
  void setTFs(tfs_ptr);
  void setGenes(genes_ptr);
  void setBindings(bindings_ptr);
  void setTFData(conc_ptr);
  void setDistances(distances_ptr);
  void setPromoters(promoters_ptr);
  
  void addNuc(int id);
  
  void create();
  void createBindings();
  void createOccupancy();
  void createSubgroups();
  void createCoeffects();
  void createQuenching();
  
  //  Getters
  tfs_ptr       getTFs();
  genes_ptr     getGenes();
  bindings_ptr  getBindings();
  conc_ptr      getTFData();
  distances_ptr getDistances();
  promoters_ptr getPromoters();
  subgroups_ptr getSubgroups();
  quenching_ptr getQuenching();
  
  vector<double>& getRate(Gene& gene) {return Rs[&gene];}
  double&         getRate(Gene& gene, int id);
  vector<int>&    getIDs() {return IDs;}

  int size() {return n;}
  
  // Methods
  bool compareTFs(tfs_ptr); // are the tfs the same as in this object?

  void saveSubgroups()    {subgroups->save();}
  void saveQuenching()    {quenching->save();}
  void saveCoeffects()    {coeffects->save();}
  
  
  void updateSubgroups()  {subgroups->update();}
  void updateQuenching()  {quenching->update();}
  void updateCoeffects()  {coeffects->update();}
  
  void restoreSubgroups() {subgroups->restore();}
  void restoreQuenching() {quenching->restore();}
  void restoreCoeffects() {coeffects->restore();}
  
  void calcQuenching()    { quenching->calc();}
  void calcCoeffects()    { coeffects->calc();}
  void calc_f()           { subgroups->calc_f();}
  
  void calcN();
  void calcR();
  void calcN(Gene&);
  void calcR(Gene&);
  
  void saveSites(TF& tf)           { bindings->saveSites(tf);          }
  void saveScores(TF& tf)          { bindings->saveScores(tf);         }
  void updateScores(TF& tf)        { bindings->updateScores(tf);       }
  void updateScores()              { bindings->updateScores();         }
  void updateSites(TF& tf)         { bindings->updateSites(tf);        }
  void updateSites()               { bindings->updateSites();          }
  void updateK(TF& tf)             { bindings->updateK(tf);            }
  void updateKandLambda(TF& tf)    { bindings->updateKandLambda(tf);   }
                                   
  void restoreScores(TF& tf)       { bindings->restoreScores(tf);      }
  void restoreSites(TF& tf)        { bindings->restoreSites(tf);       }
                                
  void saveAllOccupancy()          { bindings->saveAllOccupancy();          }
  void restoreAllOccupancy()       { bindings->restoreAllOccupancy();       }
  
  void saveModeOccupancy()    { bindings->saveModeOccupancy();    }
  void restoreModeOccupancy() { bindings->restoreModeOccupancy(); }
  
  void saveEffectiveOccupancy()    { bindings->saveEffectiveOccupancy();    }
  void restoreEffectiveOccupancy() { bindings->restoreEffectiveOccupancy(); }
  
  void updateN() {calcN();}
  void updateR() {calcR();}
  
  
  // I/O
  void printSubgroups(Gene& g, ostream& os) {subgroups->print(g, os);}
  
  void printOccupancy(Gene& g, ostream& os)          {bindings->printTotalOccupancy(g,os);    }
  void printModeOccupancy(Gene& g, ostream& os)      {bindings->printModeOccupancy(g,os);     }
  void printEffectiveOccupancy(Gene& g, ostream& os) {bindings->printEffectiveOccupancy(g,os);}
  
  void printSites(ostream& os)                     { bindings->printSites(os);}                    
  void printSites(Gene& gene, ostream& os)         { bindings->printSites(gene, os);}             
  void printSites(TF& tf, ostream& os)             { bindings->printSites(tf, os);}            
  void printSites(Gene& gene, TF& tf, ostream& os) { bindings->printSites(gene, tf, os);}
};



typedef boost::shared_ptr<Nuclei>  nuclei_ptr;



#endif
