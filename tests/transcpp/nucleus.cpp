/*********************************************************************************
*                                                                                *
*     nucleus.cpp                                                                *
*                                                                                *
*     Constains the nucleus class, which holds references to each member         *
*                                                                                *
*********************************************************************************/

#include "distance.h"
#include "nucleus.h"
#include "nuclei.h"
#include "gene.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <boost/make_shared.hpp>
#include <boost/ref.hpp>



/*    Constructors    */

Nucleus::Nucleus() {}

Nucleus::Nucleus(int i, int x, Nuclei* parent)
{
  id        = i;
  idx       = x;
  genes     = parent->getGenes();
  tfs       = parent->getTFs();
  tfdata    = parent->getTFData();
  promoters = parent->getPromoters();
  subgroups = parent->getSubgroups();
  bindings  = parent->getBindings();
  quenching = parent->getQuenching();
  
}


/*    Getters   */

double Nucleus::getN(Gene& gene)
{
  return Ns[&gene];
}

/*    Methods   */

bool Nucleus::hasSameTFs(tfs_ptr t)
{
  if (t->size() != tfs->size())
    return false;
  else
  {
    for (int i=0; i<t->size(); i++)
    {
      if (t->getTFptr(i).get() != tfs->getTFptr(i).get())
        return false;
    }
  }
  return true;
}

void Nucleus::calc_f()
{
  int ngenes = genes->size();
  subgroups->calc_f(idx);
}

void Nucleus::calcQuenching()
{
  quenching->calc(idx);
}
  


    

/*    Update    */

/*  All the update functions rescore, and generally have complementary
    restore functions   */
    

/*    Print   */

void Nucleus::printSites(ostream & os)
{
  int i, j, k;
  int w = 10;

  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    const string& gname = gene.getName();
    os << gname << endl;
    os << setprecision(3)
       << setw(w) << "name"
       << setw(w) << "start"
       << setw(w) << "end"
       << setw(w) << "score"
       << setw(w) << "K" 
       << setw(w) << "K*Kmax"
       << setw(w) << "kv"
       << setw(w) << "f" 
       << setw(w) << "ef" << endl;
    for (j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      //const string& tfname = tf.getName();
      site_ptr_vector& sites = bindings->getSites(gene, tf);
      int nsites = sites.size();
        for (k=0; k<nsites; k++)
      {
        site_ptr b = sites[k];
        os << setw(w) << b->tf->getName()
           << setw(w) << b->m 
           << setw(w) << b->n 
           << setw(w) << b->score 
           << setw(w) << b->K_exp_part
           << setw(w) << b->K_exp_part_times_kmax
           << setw(w) << b->kv[id]
           << setw(w) << b->occupancy[id]
           << setw(w) << b->effective_occupancy[id] << endl;
      }
    }
  }
}

