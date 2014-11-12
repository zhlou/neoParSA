/*********************************************************************************
*                                                                                *
*     bindings.cpp                                                               *
*                                                                                *
*     Constains a map of all site on all genes                                   *
*                                                                                *
*********************************************************************************/

#include "bindings.h"
#include "flags.h"

#include <boost/foreach.hpp>
#include <boost/range/adaptor/map.hpp>
#include <limits>

#define foreach_ BOOST_FOREACH

using namespace boost::adaptors;

Bindings::Bindings() {nnuc=0;}

void Bindings::clear()
{
  scores.clear();
  saved_scores.clear();
  sites.clear();
  saved_sites.clear();
}

map<TF*, site_ptr_vector>& Bindings::getSites(Gene& gene)
{
  return sites[&gene];
}

site_ptr_vector& Bindings::getSites(Gene& gene, TF& tf) 
{
  return sites[&gene][&tf];
}

bool Bindings::hasScores(Gene& gene, TF& tf)
{
  if (scores.find(&gene) == scores.end())
    return false;
  else
  {
    if (scores[&gene].find(&tf) == scores[&gene].end())
      return false;
  }
  return true;
}

bool Bindings::hasSites(Gene& gene, TF& tf)
{
  if (sites.find(&gene) == sites.end())
    return false;
  else
  {
    if (sites[&gene].find(&tf) == sites[&gene].end())
      return false;
  }
  return true;
}

void Bindings::createScores()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      createScores(gene, tf);
    }
  }
}

void Bindings::createSites()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      createSites(gene, tf);
    }
  }
}

void Bindings::create()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      createScores(gene, tf);
      createSites(gene, tf);
    }
  }
}

// find sites for a particular gene and tf
void Bindings::createScores(Gene& gene, TF& tf) 
{
  scores[&gene][&tf] = tf.score(gene.getSequence());
}


void Bindings::createSites(Gene& gene, TF& tf)
{
  site_ptr_vector& tmp_sites = sites[&gene][&tf];
  
  TFscore& t         = scores[&gene][&tf];
  int      len       = t.mscore.size();
  int      nmodes    = tf.getNumModes();
  double   threshold = tf.getThreshold();
  double   bsize     = tf.getBindingSize();
  double   kmax      = tf.getKmax();
  double   maxscore  = tf.getMaxScore();
  double   lambda    = tf.getLambda();
  
  vector<double>& v = conc[&tf];
  
  for (int k=0; k<len; k++)
  {
    
    if (t.mscore[k] < threshold) continue;
  
    double fscore = t.fscore[k];
    double rscore = t.rscore[k];
    
    if (fscore >= threshold)
      createSite(tmp_sites, gene, tf, k, bsize, fscore, 'F', lambda, kmax,maxscore,v,nmodes);
    if (rscore >= threshold)
      createSite(tmp_sites, gene, tf, k, bsize, rscore, 'R', lambda, kmax,maxscore,v,nmodes);

  }
  updateK(tf);
}

void Bindings::createSite(site_ptr_vector& tmp_sites, Gene& gene, TF& tf,
                          int pos, double bsize, double score, char orientation, 
                          double lambda, double kmax, double maxscore,vector<double>& v, int nmodes)
{
  site_ptr b(new BindingSite);
  b->tf                    = &tf;
  b->orientation           = orientation;
  b->m                     = pos    - bsize/2;
  b->n                     = b->m + bsize-1;
  b->score                 = score;
  b->K_exp_part            = exp((b->score - maxscore)/lambda);
  double K_exp_part_times_kmax = kmax * b->K_exp_part;
  b->K_exp_part_times_kmax = K_exp_part_times_kmax;
  
  vector<double>& kv = b->kv;
  kv.resize(nnuc);
  for (int i=0; i<nnuc; i++)
    kv[i] = K_exp_part_times_kmax * v[i];
  
  b->total_occupancy.resize(nnuc);

  vector< vector<double> >& mode_occupancy      = b->mode_occupancy;
  vector< vector<double> >& effective_occupancy = b->effective_occupancy;
  
  mode_occupancy.resize(nmodes);
  effective_occupancy.resize(nmodes);
  for (int i=0; i<nmodes; i++)
  {
    mode_occupancy[i].resize(nnuc);
    effective_occupancy[i].resize(nnuc);
  }
  tmp_sites.push_back(b);
}
      

void Bindings::updateK(TF& tf)
{
  vector<double>&  v    = conc[&tf];
  double           kmax = tf.getKmax();
  
  typedef map<TF*, site_ptr_vector> submap;
  foreach_(submap i , sites | map_values) 
  {
    site_ptr_vector& tmp_sites = i[&tf];
    int nsites = tmp_sites.size();
    for (int k=0; k<nsites; k++)
    {
      BindingSite* b = tmp_sites[k].get();
      double K_exp_part_times_kmax = kmax * b->K_exp_part;
      b->K_exp_part_times_kmax = K_exp_part_times_kmax;
      vector<double>& kv = b->kv;
      for (int j=0; j<nnuc; j++)
        kv[j] = K_exp_part_times_kmax * v[j];
    }
  }
}

void Bindings::updateKandLambda(TF& tf)
{
  vector<double>&     v = conc[&tf];
  
  double   kmax     = tf.getKmax();
  double   lambda   = tf.getLambda();
  double   maxscore = tf.getMaxScore();
  
  typedef map<TF*, site_ptr_vector> submap;
  foreach_(submap i , sites | map_values) 
  {
    site_ptr_vector& tmp_sites = i[&tf];
    int nsites = tmp_sites.size();
    for (int k=0; k<nsites; k++)
    {
      BindingSite* b = tmp_sites[k].get();
      b->K_exp_part            = exp((b->score - maxscore)/lambda);
      
      double K_exp_part_times_kmax = kmax * b->K_exp_part;
      
      b->K_exp_part_times_kmax = K_exp_part_times_kmax;
      vector<double>& kv = b->kv;
      for (int j=0; j<nnuc; j++)
        kv[j] = K_exp_part_times_kmax * v[j];
    }
  }
}

void Bindings::addNuc(int nuc_id)
{
  
  id_2_idx[nuc_id] = nnuc; 
  idx_2_id[nnuc]   = nuc_id;
  nnuc++;
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    const string& tfname = tf.getName();
    conc[&tf].push_back(tfdata->getConcByID(tfname, nuc_id, true));
  }
}



void Bindings::saveScores(TF& tf)
{
  foreach_(Gene* gene, scores | map_keys) 
    saved_scores[gene][&tf] = scores[gene][&tf];
}

void Bindings::restoreScores(TF& tf)
{
  foreach_(Gene* gene, scores | map_keys) 
    scores[gene][&tf] = saved_scores[gene][&tf];
}

void Bindings::updateScores(TF& tf)
{
  foreach_(Gene* gene, scores | map_keys)
    scores[gene][&tf] = tf.score(gene->getSequence());
}

void Bindings::updateScores()
{
  foreach_(Gene* gene, scores | map_keys)
  {
    foreach_(TF* tf, scores[gene] | map_keys)
      scores[gene][tf] = tf->score(gene->getSequence());
  }
}

void Bindings::updateSites()
{
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    updateSites(tf);
  }
}

void Bindings::updateSites(TF& tf)
{

  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    site_ptr_vector& gene_sites = sites[&gene][&tf];
    gene_sites.clear();
    createSites(gene,tf);
  }
}

void::Bindings::saveSites(TF& tf)
{
  foreach_(Gene* gene, scores | map_keys)
    saved_sites[gene][&tf] = sites[gene][&tf];
}

void::Bindings::restoreSites(TF& tf)
{
  foreach_(Gene* gene, scores | map_keys)
    sites[gene][&tf] = saved_sites[gene][&tf];
}


void Bindings::saveOccupancy()
{
  typedef map<TF*, site_ptr_vector>& submap;
  foreach_(submap i , sites | map_values) 
  {
    foreach_(site_ptr_vector& j, i | map_values)
    {
      int nsites = j.size();
      for (int k=0; k<nsites; k++)
        j[k]->saved_total_occupancy = j[k]->total_occupancy;
    }
  }
}
        

void Bindings::restoreOccupancy()
{
  typedef map<TF*, site_ptr_vector>& submap;
  foreach_(submap i , sites | map_values) 
  {
    foreach_(site_ptr_vector& j, i | map_values)
    {
      int nsites = j.size();
      for (int k=0; k<nsites; k++)
        j[k]->total_occupancy = j[k]->saved_total_occupancy;
    }
  }
}


/* this function saves effective occupancy and resets occupancy to 
simple fractional occupancy */
void Bindings::saveEffectiveOccupancy()
{
  typedef map<TF*, site_ptr_vector>& submap;
  foreach_(submap i , sites | map_values) 
  {
    foreach_(site_ptr_vector& j, i | map_values)
    {
      int nsites = j.size();
      for (int k=0; k<nsites; k++)
      {
        j[k]->saved_effective_occupancy = j[k]->effective_occupancy;
      }
    }
  }
}


void Bindings::restoreEffectiveOccupancy()
{
  typedef map<TF*, site_ptr_vector>& submap;
  foreach_(submap i , sites | map_values) 
  {
    foreach_(site_ptr_vector& j, i | map_values)
    {
      int nsites = j.size();
      for (int k=0; k<nsites; k++)
        j[k]->effective_occupancy = j[k]->saved_effective_occupancy;
    }
  }
}

void Bindings::saveModeOccupancy()
{
  typedef map<TF*, site_ptr_vector>& submap;
  foreach_(submap i , sites | map_values) 
  {
    foreach_(site_ptr_vector& j, i | map_values)
    {
      int nsites = j.size();
      for (int k=0; k<nsites; k++)
      {
        j[k]->saved_effective_occupancy = j[k]->effective_occupancy;
        j[k]->saved_mode_occupancy      = j[k]->mode_occupancy;
      }
    }
  }
}


void Bindings::restoreModeOccupancy()
{
  typedef map<TF*, site_ptr_vector>& submap;
  foreach_(submap i , sites | map_values) 
  {
    foreach_(site_ptr_vector& j, i | map_values)
    {
      int nsites = j.size();
      for (int k=0; k<nsites; k++)
      {
        j[k]->effective_occupancy = j[k]->saved_effective_occupancy;
        j[k]->mode_occupancy      = j[k]->saved_mode_occupancy;
      }
    }
  }
}


void Bindings::saveAllOccupancy()
{
  typedef map<TF*, site_ptr_vector>& submap;
  foreach_(submap i , sites | map_values) 
  {
    foreach_(site_ptr_vector& j, i | map_values)
    {
      int nsites = j.size();
      for (int k=0; k<nsites; k++)
      {
        site_ptr b = j[k];
        b->saved_kv                  = b->kv;
        b->saved_effective_occupancy = b->effective_occupancy;
        b->saved_mode_occupancy      = b->mode_occupancy;
        b->saved_total_occupancy     = b->total_occupancy;
      }
    }
  }
}

void Bindings::restoreAllOccupancy()
{
  typedef map<TF*, site_ptr_vector>& submap;
  foreach_(submap i , sites | map_values) 
  {
    foreach_(site_ptr_vector& j, i | map_values)
    {
      int nsites = j.size();
      for (int k=0; k<nsites; k++)
      {
        site_ptr b = j[k];
        b->kv                  = b->saved_kv;
        b->mode_occupancy      = b->saved_mode_occupancy;
        b->effective_occupancy = b->saved_effective_occupancy;
        b->total_occupancy     = b->saved_total_occupancy;
      }
    }
  }
}


void Bindings::printSites(ostream& os)
{
  foreach_(Gene* gene, sites | map_keys)
    printSites(*gene, os);
}

void Bindings::printSites(Gene& gene, ostream& os)
{
  os << gene.getName() << endl;
  printSiteHeader(os);
  foreach_(TF* tf, sites[&gene] | map_keys)
  {
    site_ptr_vector& tmp_sites = sites[&gene][tf];
    int nsites = tmp_sites.size();
    for (int i=0; i<nsites; i++)
    {
      printSite(tmp_sites[i].get(), os);
      os << endl;
    }
  }
}

void Bindings::printSites(TF& tf, ostream& os)
{
  foreach_(Gene* gene, sites | map_keys)
  {
    os << gene->getName() << endl;
    printSiteHeader(os);

    site_ptr_vector& tmp_sites = sites[gene][&tf];
    int nsites = tmp_sites.size();
    for (int i=0; i<nsites; i++)
    {
      printSite(tmp_sites[i].get(), os);
      os << endl;
    }
  }
}

void Bindings::printSites(Gene& gene, TF& tf, ostream& os)
{
  os << gene.getName() << endl;
  printSiteHeader(os);
  site_ptr_vector& tmp_sites = sites[&gene][&tf];
  int nsites = tmp_sites.size();
  for (int i=0; i<nsites; i++)
  {
    printSite(tmp_sites[i].get(), os);
    os << endl;
  }
}

void Bindings::printTotalOccupancy(Gene& gene, ostream& os)
{
  int w = 10;
  os << setprecision(3);
  
  // loop through everything and figure out the header
  os << setw(w) << "id";
  
  // find all the ids, using map to prevent duplicates
  // simultaneously print binding site headers
  foreach_(TF* tf, sites[&gene] | map_keys)
  {
    site_ptr_vector& tmp_sites = sites[&gene][tf];
    int nsites = tmp_sites.size();
    for (int i=0; i<nsites; i++)
    {
      string header = tf->getName();
      stringstream convert;
      convert << i;
      string num = convert.str();
      header += num;
      //header << i;
      os << setw(w) <<  header;
    }
  }
  
  os << endl;
  
  for (int nuc=0; nuc<nnuc; nuc++)
  {
    int id = idx_2_id[nuc];
    os << setw(w) << id;
    
    foreach_(TF* tf, sites[&gene] | map_keys)
    {
      site_ptr_vector& tmp_sites = sites[&gene][tf];
      int nsites = tmp_sites.size();
      for (int i=0; i<nsites; i++)
        os << setw(w) << tmp_sites[i]->total_occupancy[nuc];
    }
    os << endl;
  }
}

void Bindings::printEffectiveOccupancy(Gene& gene, ostream& os)
{
    int w = 10;
  os << setprecision(3);
  
  // loop through everything and figure out the header
  os << setw(w) << "id";
  
  // find all the ids, using map to prevent duplicates
  // simultaneously print binding site headers
  foreach_(TF* tf, sites[&gene] | map_keys)
  {
    site_ptr_vector& tmp_sites = sites[&gene][tf];
    int nsites = tmp_sites.size();
    int nmodes = tf->getNumModes();
    for (int i=0; i<nsites; i++)
    {      
      for (int j=0; j<nmodes; j++)
      {
        string header = tf->getName();
        stringstream convert;
        convert << i;
        if (j>0) convert << "_" << j+1;
        string num = convert.str();
        header += num;
        //header << i;
        os << setw(w) <<  header;
      }
    }
  }
  
  os << endl;
  
  for (int nuc=0; nuc<nnuc; nuc++)
  {
    int id = idx_2_id[nuc];
    os << setw(w) << id;
   
    foreach_(TF* tf, sites[&gene] | map_keys)
    {
      site_ptr_vector& tmp_sites = sites[&gene][tf];
      int nsites = tmp_sites.size();
      int nmodes = tf->getNumModes();
      for (int i=0; i<nsites; i++)
      {
        for (int j=0; j<nmodes; j++)
          os << setw(w) << tmp_sites[i]->effective_occupancy[j][nuc];
      }
    }
    os << endl;
  }
}

void Bindings::printModeOccupancy(Gene& gene, ostream& os)
{
    int w = 10;
  os << setprecision(3);
  
  // loop through everything and figure out the header
  os << setw(w) << "id";
  
  // find all the ids, using map to prevent duplicates
  // simultaneously print binding site headers
  foreach_(TF* tf, sites[&gene] | map_keys)
  {
    site_ptr_vector& tmp_sites = sites[&gene][tf];
    int nsites = tmp_sites.size();
    int nmodes = tf->getNumModes();
    for (int i=0; i<nsites; i++)
    {
      for (int j=0; j<nmodes; j++)
      {
        string header = tf->getName();
        stringstream convert;
        convert << i;
        if (j>0) convert << "_" << j+1;
        string num = convert.str();
        header += num;
        //header << i;
        os << setw(w) <<  header;
      }
    }
  }
  
  os << endl;
  
  for (int nuc=0; nuc<nnuc; nuc++)
  {
    int id = idx_2_id[nuc];
    os << setw(w) << id;
   
    foreach_(TF* tf, sites[&gene] | map_keys)
    {
      site_ptr_vector& tmp_sites = sites[&gene][tf];
      int nsites = tmp_sites.size();
      int nmodes = tf->getNumModes();
      for (int i=0; i<nsites; i++)
      {
        for (int j=0; j<nmodes; j++)
          os << setw(w) << tmp_sites[i]->mode_occupancy[j][nuc];
      }
    }
    os << endl;
  }
}
