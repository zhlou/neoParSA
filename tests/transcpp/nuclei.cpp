/*********************************************************************************
*                                                                                *
*     nuclei.cpp                                                                 *
*                                                                                *
*     A nuclei object contains all nuclei which share the same TFs. This means   *
*     they share subgroups and quenching interactions.                           *
*                                                                                *
*     In the original code, many of the calculations on occupancy and quenching  *
*     were multiplications by 0. This class was created to prevent such needless *
*     calculations.                                                              *
*                                                                                *
*********************************************************************************/

#include "nuclei.h"

#include <boost/foreach.hpp>
#include <boost/range/adaptor/map.hpp>
#include <limits>

using namespace boost::adaptors;
# define foreach_ BOOST_FOREACH

/*    Constructors    */

Nuclei::Nuclei(Organism* parent, tfs_ptr t) :
  bindings(bindings_ptr(new Bindings)),
  subgroups(subgroups_ptr(new Subgroups)),
  quenching(quenching_ptr(new QuenchingInteractions)),
  coeffects(modifying_ptr(new ModifyingInteractions))
{
  n         = 0;
  tfs       = t;
  genes     = parent->getGenes();
  tfdata    = parent->getTFData();
  distances = parent->getDistances();
  promoters = parent->getPromoters();
  mode      = parent->getMode();
}
  

/*    Setters    */

void Nuclei::setTFs(tfs_ptr x)             {tfs       = x;}
void Nuclei::setGenes(genes_ptr x)         {genes     = x;}
void Nuclei::setBindings(bindings_ptr x)   {bindings  = x;}
void Nuclei::setTFData(conc_ptr x)         {tfdata    = x;}
void Nuclei::setDistances(distances_ptr x) {distances = x;}
void Nuclei::setPromoters(promoters_ptr x) {promoters = x;}

void Nuclei::create()
{
  createBindings();
  createSubgroups();
  createCoeffects();
  createQuenching();

  calc_f();
  calcCoeffects();
  calcQuenching();
  calcN();
  calcR();
}

void Nuclei::createBindings()
{
  bindings->setGenes(genes);
  bindings->setTFs(tfs);
  bindings->setTFData(tfdata); 
  
  for (int i=0; i<IDs.size(); i++)
    bindings->addNuc(IDs[i]);
  
  bindings->createScores();
  bindings->createSites();
}
   
void Nuclei::createSubgroups()
{
  subgroups->create(genes, tfs, bindings, mode);
}

void Nuclei::createCoeffects()
{
  coeffects->create(genes, tfs, bindings, distances);
}

void Nuclei::createQuenching()
{
  quenching->create(genes, tfs, bindings, distances);
}
  
void Nuclei::addNuc(int id)
{
  IDs.push_back(id);
  n=IDs.size();
  
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    Ns[&gene].resize(n);
    Rs[&gene].resize(n);
  }
    
}

/*    Getters   */

tfs_ptr       Nuclei::getTFs()          {return tfs;      }
genes_ptr     Nuclei::getGenes()        {return genes;    }
bindings_ptr  Nuclei::getBindings()     {return bindings; }
conc_ptr      Nuclei::getTFData()       {return tfdata;   }
distances_ptr Nuclei::getDistances()    {return distances;}
promoters_ptr Nuclei::getPromoters()    {return promoters;}
quenching_ptr Nuclei::getQuenching()    {return quenching;}
subgroups_ptr Nuclei::getSubgroups()    {return subgroups;}
  
double& Nuclei::getRate(Gene& gene, int id)
{
  int nids = IDs.size();
  for (int i=0; i<nids; i++)
  {
    if (id == IDs[i])
    {
      return Rs[&gene][i];
    }
  }
  cerr << "ERROR: could not find id in this set of nuclei!" << endl;
  exit(1);
}

/*    Methods   */

bool Nuclei::compareTFs(tfs_ptr t)
{
  if (tfs->size() != t->size())
    return false;
  else
  {
    for (int i=0; i<tfs->size(); i++)
    {
      if (tfs->getTFptr(i) != t->getTFptr(i))
        return false;
    }
  }
  return true;
}

void Nuclei::calcN()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    calcN(gene);
  }
}


void Nuclei::calcN(Gene& gene)
{
  vector<double>& tNs = Ns[&gene];
  tNs.resize(n);
  for (int i=0; i<n;i++)
    tNs[i]=0;
  
  int ntfs = tfs->size();
  
  for (int i = 0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    
    if (tf.neverActivates()) continue;
    
    vector<double> coefs    = tf.getCoefs();
    site_ptr_vector& tsites = bindings->getSites(gene,tf);
    int ntsites = tsites.size();
    int nmodes  = coefs.size();
    for (int j=0; j<nmodes; j++)
    {
      double efficiency = coefs[j];
      if (efficiency <= 0) continue;

      for (int k=0; k<ntsites; k++)
      {
        vector<double>& eff_occ = tsites[k]->effective_occupancy[j];
        for (int l=0; l<n; l++)
        {
          tNs[l] += eff_occ[l]*efficiency;
        }
      }
    }
  }
}

void Nuclei::calcR()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    calcR(gene);
  }
}

void Nuclei::calcR(Gene& gene)
{
  vector<double>& tNs = Ns[&gene];
  vector<double>& tRs = Rs[&gene];
  
  for (int i=0; i<n; i++)
    tRs[i] = gene.getRate(tNs[i]);
}

