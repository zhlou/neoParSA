/*********************************************************************************
*                                                                                *
*     bindings.cpp                                                               *
*                                                                                *
*     Constains a map of all site on all genes                                   *
*                                                                                *
*********************************************************************************/

#include "bindings.h"

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <limits>

#define foreach_ BOOST_FOREACH
#define to_string_ boost::lexical_cast<string>

Bindings::Bindings() {nnuc=0;}

void Bindings::clear()
{
  scores.clear();
  saved_scores.clear();
  sites.clear();
  saved_sites.clear();
}

boost::unordered_map<TF*, site_ptr_vector>& Bindings::getSites(Gene& gene)
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

void Bindings::createScores(Gene& gene)
{
  int ntfs   = tfs->size();
  
  for (int j=0; j<ntfs; j++)
  {
    TF& tf = tfs->getTF(j);
    createScores(gene, tf);
  }
}

void Bindings::createSites()
{
  int ngenes = genes->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    createSites(gene);
  }
}

void Bindings::createSites(Gene& gene)
{
  int ntfs   = tfs->size();

  for (int j=0; j<ntfs; j++)
  {
    TF& tf = tfs->getTF(j);
    createSites(gene, tf);
  }
  order_sites(gene);
  //for (int j=0; j<ntfs; j++)
  //{
  //  TF& tf = tfs->getTF(j);
  //  //cerr << "verify after creating 1 for " << tf.getName() << endl;
  //  // verify_order(gene, tf);
  //}
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
    order_sites(gene);
    //for (int j=0; j<ntfs; j++)
    //{
    //  TF& tf = tfs->getTF(j);
    //  //cerr << "verify after creating 2 for " << tf.getName() << endl;
    //  // verify_order(gene, tf);
    //}
  }
}

void Bindings::create(Gene& gene)
{
  int ntfs   = tfs->size();

  for (int j=0; j<ntfs; j++)
  {
    TF& tf = tfs->getTF(j);
    createScores(gene, tf);
    createSites(gene, tf);
  }
  order_sites(gene);
  //for (int j=0; j<ntfs; j++)
  //  {
  //    TF& tf = tfs->getTF(j);
  //    //cerr << "verify after creating 3 for " << tf.getName() << endl;
  //    verify_order(gene, tf);
  //  }
}

// find sites for a particular gene and tf
void Bindings::createScores(Gene& gene, TF& tf) 
{
  scores[&gene][&tf] = tf.score(gene.getSequence());
}


void Bindings::createSites(Gene& gene, TF& tf)
{
  site_ptr_vector& tmp_sites = sites[&gene][&tf];
  if (tmp_sites.size() != 0) error("you didnt clear tmp_sites!");
  
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
  if (!mode->getSelfCompetition())
    trimOverlaps(gene,tf);
  
  updateK(gene, tf);
}

void Bindings::add_to_ordered(Gene& gene, TF& tf)
{
  site_ptr_vector& tf_sites    = sites[&gene][&tf];
  vector<BindingSite*>& fsites = ordered_sites_f[&gene];
  
  int ntfsites = tf_sites.size();
  int nfsites  = fsites.size();
  int total    = nfsites + ntfsites;
  fsites.resize(total);
  
  int i, j;
  for (i=nfsites, j=0; j<ntfsites; i++, j++)
    fsites[i] = tf_sites[j].get();
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
  
  b->index_in_site_map = tmp_sites.size();
  tmp_sites.push_back(b);
  //ordered_sites_f[&gene].push_back(b.get());
}


void Bindings::trimOverlaps(Gene& gene, TF& tf)
{
  site_ptr_vector& tf_sites = sites[&gene][&tf];
  
  int nsites = tf_sites.size();
  int i = 0;
  
  site_ptr_vector::iterator sites_begin = tf_sites.begin();
  while (i<(nsites-1))
  {
    BindingSite& site1 = *tf_sites[i];
    BindingSite& site2 = *tf_sites[i+1];
    
    if (bad_overlap_function(site1, site2))
    {
      if (site1.K_exp_part > site2.K_exp_part)
      {
        tf_sites.erase(sites_begin + i + 1);
        nsites--;
      } 
      else
      {
        tf_sites.erase(sites_begin + i);
        nsites--;
      }
    } 
    else
      i++;
  }
  
  // now we need to set the index in bindings back
  nsites = tf_sites.size();
  for (int i=0; i<nsites; i++)
  {
    BindingSite& b = *tf_sites[i];
    b.index_in_site_map = i;
  }
}
    
  
  
  
void Bindings::order_sites(Gene& gene)
{
  vector<BindingSite*>& fsites = ordered_sites_f[&gene];
  vector<BindingSite*>& rsites = ordered_sites_r[&gene];
  
  fsites.clear();
  int ntfs = tfs->size();
  boost::unordered_map<TF*, site_ptr_vector>& gene_sites = sites[&gene];
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tf_sites = gene_sites[&tf];
    int nsites = tf_sites.size();
    int cur_total = fsites.size();
    fsites.resize(cur_total+nsites);
    for (int j=0; j<nsites; j++)
      fsites[cur_total + j] = tf_sites[j].get();
  }
  
  // sort forward sites
  
  std::sort(fsites.begin(), fsites.end(), compareBindingSiteRight);

  // sort reverse sites
  int nsites = fsites.size();
  
  rsites.resize(nsites);
  for (int i=0; i<nsites; i++)
    rsites[i] = fsites[nsites - i - 1];
  
  std::sort(rsites.begin(), rsites.end(), compareBindingSiteLeft);
 
  // in general, we want to know where sites are in each structure, so we will keep an index of each
  for (int i=0; i<nsites; i++)
  {
    fsites[i]->index_in_ordered_f = i;
    rsites[i]->index_in_ordered_r = i;
  }
}

void Bindings::updateK(Gene& gene, TF& tf)
{
  vector<double>&  v    = conc[&tf];
  double           kmax = tf.getKmax();
  
  site_ptr_vector& tmp_sites = sites[&gene][&tf];
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

void Bindings::updateK(TF& tf)
{
  vector<double>&  v    = conc[&tf];
  double           kmax = tf.getKmax();
  
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    site_ptr_vector& tmp_sites = sites[&gene][&tf];
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

void Bindings::updateKandLambda(Gene& gene, TF& tf)
{
  vector<double>&     v = conc[&tf];
  
  double   kmax     = tf.getKmax();
  double   lambda   = tf.getLambda();
  double   maxscore = tf.getMaxScore();
  
  site_ptr_vector& tmp_sites = sites[&gene][&tf];
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

void Bindings::updateKandLambda(TF& tf)
{
  vector<double>&     v = conc[&tf];
  
  double   kmax     = tf.getKmax();
  double   lambda   = tf.getLambda();
  double   maxscore = tf.getMaxScore();
  
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    site_ptr_vector& tmp_sites = sites[&gene][&tf];
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

void Bindings::addNuc(string& nuc_id)
{
  
  //id_2_idx[nuc_id] = nnuc; 
  //idx_2_id[nnuc]   = nuc_id;
  
  IDs.push_back(nuc_id);
  nnuc++;
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    const string& tfname = tf.getName();
    conc[&tf].push_back(tfdata->getDataPoint("TF",tfname, "ID",nuc_id));
  }
}



void Bindings::saveScores(TF& tf)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    saved_scores[&gene][&tf] = scores[&gene][&tf];
  }
}

void Bindings::saveScores(Gene& gene, TF& tf)
{
  saved_scores[&gene][&tf] = scores[&gene][&tf];
}

void Bindings::restoreScores(TF& tf)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    scores[&gene][&tf] = saved_scores[&gene][&tf];
  }
}

void Bindings::restoreScores(Gene& gene, TF& tf)
{
  scores[&gene][&tf] = saved_scores[&gene][&tf];
}

void Bindings::updateScores(Gene& gene, TF& tf)
{
  scores[&gene][&tf] = tf.score(gene.getSequence());
}

void Bindings::updateScores(TF& tf)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    scores[&gene][&tf] = tf.score(gene.getSequence());
  }
}

void Bindings::updateScores()
{
  int ntfs   = tfs->size();
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      scores[&gene][&tf] = tf.score(gene.getSequence());
    }
  }
}

void Bindings::updateScores(Gene& gene)
{
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    scores[&gene][&tf] = tf.score(gene.getSequence());
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

void Bindings::updateSites(Gene& gene)
{
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    updateSites(gene, tf);
  }
}

void Bindings::updateSites(TF& tf)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    updateSites(gene,tf);
  }
}

void Bindings::eraseTF(Gene& gene, TF& tf)
{
  /*
  vector<BindingSite*>& fsites = ordered_sites_f[&gene];
  vector<BindingSite*>& rsites = ordered_sites_r[&gene];
  
  // these are populated by base, so they are in order!
  site_ptr_vector& tf_sites = sites[&gene][&tf];
  
  vector<BindingSite*>::iterator fbegin = fsites.begin();
  vector<BindingSite*>::iterator rbegin = rsites.begin();
  
  int ntfsites = tf_sites.size();
  int i, j;
  for (i=0, j=(ntfsites-1); i<ntfsites; i++, j--)
  {
    BindingSite& fsite = *tf_sites[j];
    BindingSite& rsite = *tf_sites[i];

    fsites.erase(fbegin + fsite.index_in_ordered_f);
    rsites.erase(rbegin + rsite.index_in_ordered_r);
  }
  */
  site_ptr_vector& tf_sites = sites[&gene][&tf];
  tf_sites.clear();
} 
  

  
void Bindings::updateSites(Gene& gene, TF& tf)
{
  //cerr << "verify before update for " << tf.getName() << endl;
  //verify_order(gene, tf);
  eraseTF(gene,tf);
  createSites(gene,tf);
  order_sites(gene);
  //cerr << "verify after update for " << tf.getName() << endl;
  //verify_order(gene, tf);
}

/* 
a function to verify that the ordered sites and tfs all point to the right 
places. We need to do do comparisons to do wo. First, we need to verify that 
everything in the site_ptr_vector is in the ordered sites at the right place,
second, we need to verify that all the ordered sites point to valid objects that
are in the site_ptr_vector
*/

void Bindings::verify_order(Gene& gene, TF& tf)
{
  vector<BindingSite*>& fsites = ordered_sites_f[&gene];                       
  vector<BindingSite*>& rsites = ordered_sites_r[&gene];                         
  site_ptr_vector& tf_sites    = sites[&gene][&tf];
  
  // verify that ever site points to something in the ordered list
  int nsites = tf_sites.size();
  //cerr << nsites << " for tf " << tf.getName() << endl;
  for (int j=0; j<nsites; j++)
  {
    BindingSite& b = *tf_sites[j];
    // check that this index exists in the sites
    if (b.index_in_ordered_f >= (int) fsites.size())
      error("Tried to access index " + to_string_(b.index_in_ordered_f) + " of " + to_string_(fsites.size()) + " from forward ordered sites");
    if (b.index_in_ordered_r >= (int) rsites.size())
      error("Tried to access index " + to_string_(b.index_in_ordered_r) + " of " + to_string_(rsites.size()) + " from reverse ordered sites");
    // check that this points to itself
    if (b.index_in_site_map != j)
      error(" index in site map malformed for site " + to_string_(j) + " of tf " + tf.getName());
    
    if (fsites[b.index_in_ordered_f] != &b)
      error(" index in ordered forward malformed for site " + to_string_(j) + " of tf " + tf.getName());
    
    if (rsites[b.index_in_ordered_r] != &b)
      error(" index in ordered reverse malformed for site " + to_string_(j) + " of tf " + tf.getName());
  }
  
  // verify that the ordered sites are formed correctly
  nsites = fsites.size();
  if (nsites != (int) rsites.size())
    error("different numbers of sites in forward and reverse ordered lists");
  
  for (int i=0; i<nsites; i++)
  {
    BindingSite& sitef = *fsites[i];
    BindingSite& siter = *rsites[i];
    
    TF& tff = *sitef.tf;
    TF& tfr = *siter.tf;
    
    if (&tff == &tf)
    {
      if (sitef.index_in_site_map >= (int) tf_sites.size())
        error("Tried to access index " + to_string_(sitef.index_in_site_map) + " of " + to_string_(tf_sites.size()) + " from sites for tf " + tff.getName());
      if (tf_sites[sitef.index_in_site_map].get() != &sitef)
        error(" ordered f site " + to_string_(i) + " does not point to a valid site" + to_string_(sitef.index_in_ordered_f));
    }
    if (&tfr == &tf)
    {
      if (siter.index_in_site_map >= (int) tf_sites.size())
        error("Tried to access index " + to_string_(siter.index_in_site_map) + " of " + to_string_(tf_sites.size()) + " from sites for tf " + tfr.getName());
      if (tf_sites[siter.index_in_site_map].get() != &siter)
        error(" ordered r site " + to_string_(i) + " does not point to a valid site " + to_string_(siter.index_in_ordered_r));
    }
  }
}
    
    

void::Bindings::saveSites(TF& tf)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    saveSites(gene, tf);
  }
}

void::Bindings::saveSites(Gene& gene, TF& tf)
{
  saved_sites[&gene][&tf] = sites[&gene][&tf];
}

void::Bindings::restoreSites(TF& tf)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    restoreSites(gene, tf);
  }
    
}

void::Bindings::restoreSites(Gene& gene, TF& tf)
{
  //for (int i=0; i<tfs->size(); i++)
  //{
  //  TF& tf_test = tfs->getTF(i);
  // // printSites(gene,tf_test,cerr);
  //  //cerr << "verify before restore for " << tf_test.getName() << endl;
  //  verify_order(gene, tf_test);
  //}
  sites[&gene][&tf] = saved_sites[&gene][&tf];
  order_sites(gene); 
  //cerr << "fsites.size() = " << ordered_sites_f[&gene].size() << endl;
  //cerr << "rsites.size() = " << ordered_sites_r[&gene].size() << endl;
  
  //for (int i=0; i<tfs->size(); i++)
  //{
  //  TF& tf_test = tfs->getTF(i);
  //  //printSites(gene,tf_test,cerr);
  //  //cerr << "verify after restore for " << tf_test.getName() << endl;
  //  verify_order(gene, tf_test);
  //}
}


void Bindings::saveOccupancy()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      site_ptr_vector& tfsites = gsites[&tf];
      int nsites = tfsites.size();
      for (int k=0; k<nsites; k++)
      {
        tfsites[k]->saved_total_occupancy = tfsites[k]->total_occupancy;
      }
    }
  }
}


void Bindings::saveOccupancy(Gene& gene)
{
  boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
      tfsites[j]->saved_total_occupancy = tfsites[j]->total_occupancy;
  }
}
        

void Bindings::restoreOccupancy()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      site_ptr_vector& tfsites = gsites[&tf];
      int nsites = tfsites.size();
      for (int k=0; k<nsites; k++)
      {
        tfsites[k]->total_occupancy = tfsites[k]->saved_total_occupancy;
      }
    }
  }
}

void Bindings::restoreOccupancy(Gene& gene)
{
  boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
      tfsites[j]->total_occupancy = tfsites[j]->saved_total_occupancy;
  }
}


/* this function saves effective occupancy and resets occupancy to 
simple fractional occupancy */
void Bindings::saveEffectiveOccupancy()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      site_ptr_vector& tfsites = gsites[&tf];
      int nsites = tfsites.size();
      for (int k=0; k<nsites; k++)
      {
        tfsites[k]->saved_effective_occupancy = tfsites[k]->effective_occupancy;
      }
    }
  }
}

void Bindings::saveEffectiveOccupancy(Gene& gene)
{
  boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
      tfsites[j]->saved_effective_occupancy = tfsites[j]->effective_occupancy;
  }
}

void Bindings::restoreEffectiveOccupancy()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      site_ptr_vector& tfsites = gsites[&tf];
      int nsites = tfsites.size();
      for (int k=0; k<nsites; k++)
      {
        tfsites[k]->effective_occupancy = tfsites[k]->saved_effective_occupancy;
      }
    }
  }
}

void Bindings::restoreEffectiveOccupancy(Gene& gene)
{
  boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
      tfsites[j]->effective_occupancy = tfsites[j]->saved_effective_occupancy;
  }
}

void Bindings::saveModeOccupancy()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      site_ptr_vector& tfsites = gsites[&tf];
      int nsites = tfsites.size();
      for (int k=0; k<nsites; k++)
      {
        tfsites[k]->saved_effective_occupancy = tfsites[k]->effective_occupancy;
        tfsites[k]->saved_mode_occupancy      = tfsites[k]->mode_occupancy;
      }
    }
  }
}

void Bindings::saveModeOccupancy(Gene& gene)
{
  boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
    {
      tfsites[j]->saved_effective_occupancy = tfsites[j]->effective_occupancy;
      tfsites[j]->saved_mode_occupancy      = tfsites[j]->mode_occupancy;
    }
  }
}


void Bindings::restoreModeOccupancy()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      site_ptr_vector& tfsites = gsites[&tf];
      int nsites = tfsites.size();
      for (int k=0; k<nsites; k++)
      {
        tfsites[k]->effective_occupancy = tfsites[k]->saved_effective_occupancy;
        tfsites[k]->mode_occupancy      = tfsites[k]->saved_mode_occupancy;
      }
    }
  }
}

void Bindings::restoreModeOccupancy(Gene& gene)
{
  boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
    {
      tfsites[j]->effective_occupancy = tfsites[j]->saved_effective_occupancy;
      tfsites[j]->mode_occupancy      = tfsites[j]->saved_mode_occupancy;
    }
  }
}


void Bindings::saveAllOccupancy()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      site_ptr_vector& tfsites = gsites[&tf];
      int nsites = tfsites.size();
      for (int k=0; k<nsites; k++)
      {
        site_ptr b = tfsites[k];
        b->saved_kv                  = b->kv;
        b->saved_effective_occupancy = b->effective_occupancy;
        b->saved_mode_occupancy      = b->mode_occupancy;
        b->saved_total_occupancy     = b->total_occupancy;
      }
    }
  }
}

void Bindings::saveAllOccupancy(Gene& gene)
{
  boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
    {
      site_ptr b = tfsites[j];
      b->saved_kv                  = b->kv;
      b->saved_effective_occupancy = b->effective_occupancy;
      b->saved_mode_occupancy      = b->mode_occupancy;
      b->saved_total_occupancy     = b->total_occupancy;
    }
  }
}

void Bindings::restoreAllOccupancy()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      site_ptr_vector& tfsites = gsites[&tf];
      int nsites = tfsites.size();
      for (int k=0; k<nsites; k++)
      {
        site_ptr b = tfsites[k];
        b->kv                  = b->saved_kv;
        b->mode_occupancy      = b->saved_mode_occupancy;
        b->effective_occupancy = b->saved_effective_occupancy;
        b->total_occupancy     = b->saved_total_occupancy;
      }
    }
  }
}

void Bindings::restoreAllOccupancy(Gene& gene)
{
  boost::unordered_map<TF*, site_ptr_vector>& gsites = sites[&gene];
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
    {
      site_ptr b = tfsites[j];
      b->kv                  = b->saved_kv;
      b->mode_occupancy      = b->saved_mode_occupancy;
      b->effective_occupancy = b->saved_effective_occupancy;
      b->total_occupancy     = b->saved_total_occupancy;
    }
  }
}


bool Bindings::isEqual(Bindings& test_bindings)
{
  // make sure we have the same number of tfs and genes
  genes_ptr test_genes = test_bindings.getGenes();
  tfs_ptr   test_tfs   = test_bindings.getTFs();
  
  int n_test_genes = test_genes->size();
  int n_test_tfs   = test_tfs->size();
  
  int n_genes = genes->size();
  int n_tfs   = genes->size();
  
  if (n_test_genes != n_genes) warning("Number of genes not equal!");
  if (n_test_tfs != n_tfs) warning("Number of tfs not equal!");
  
  // make sure they are the same tfs and genes
  for (int i=0; i<n_genes; i++)
  {
    Gene& test_gene = test_genes->getGene(i);
    Gene& gene      = genes->getGene(i);
    if (test_gene.getName() != gene.getName()) warning("Genes not the same!");
    for (int j=0; j<n_tfs; j++)
    {
      TF& test_tf = test_tfs->getTF(j);
      TF& tf      = tfs->getTF(j);
      if (test_tf.getName() != tf.getName()) warning("TFs not the same!");
      TFscore& test_score = test_bindings.getScores(test_gene, test_tf);
      TFscore& score      = getScores(gene, tf);
      int n = test_score.fscore.size();
      for (int k=0; k<n; k++)
      {
        if (test_score.fscore[k] != score.fscore[k]) warning("fscores not equal");
        if (test_score.rscore[k] != score.rscore[k]) warning("rscores not equal");
      }
    }
    vector<BindingSite*>& test_sites_f = test_bindings.getFsites(test_gene);
    vector<BindingSite*>& test_sites_r = test_bindings.getRsites(test_gene);
    
    int n = ordered_sites_f.size();
    //for (int j=0; j<n; j++)
    //{
    //}
      
  }
    
  // first check if scores are the same
}
  

void Bindings::printSites(ostream& os)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    printSites(gene, os);
  }
}

void Bindings::printSites(Gene& gene, ostream& os)
{
  os << gene.getName() << endl;
  int print_source = 0;
  int ntfs = tfs->size();
  
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    const string& pwm_source = tf.getPWMSource();
    if (pwm_source != string(""))
      print_source = max(print_source, (int) pwm_source.size());
  }
  
  int p = mode->getPrecision();
  int w = p+7;
  os << setprecision(p)
     << setw(w) << "name"
     << setw(w) << "index"
     << setw(w) << "index_f"
     << setw(w) << "index_r"
     << setw(w) << "start"
     << setw(w) << "end"
     << setw(w) << "score"
     << setw(w) << "K"
     << setw(w) << "K*kmax"
     << setw(w) << "strand";
     
  if (print_source)
    os << setw(print_source+2) << "source";
  
  os << endl;
  for (int i=0; i<ntfs; i++)
    printSites(gene, tfs->getTF(i), os, p, print_source);
}

void Bindings::printSites(TF& tf, ostream& os)
{

  int print_source = 0;
  const string& pwm_source = tf.getPWMSource();
  if (pwm_source != string(""))
    print_source = max(print_source, (int) pwm_source.size());
  
  int p = mode->getPrecision();
  int w = p+7;
  
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    os << gene.getName() << endl;
    os << setprecision(p)
     << setw(w) << "name"
     << setw(w) << "index"
     << setw(w) << "index_f"
     << setw(w) << "index_r"
     << setw(w) << "start"
     << setw(w) << "end"
     << setw(w) << "score"
     << setw(w) << "K"
     << setw(w) << "K*kmax"
     << setw(w) << "strand";
     
    if (print_source)
      os << setw(print_source+2) << "source";
    
    os << endl;
    printSites(gene, tf, os, p, print_source);
  }
}

void Bindings::printSites(Gene& gene, TF& tf, ostream& os)
{
  int print_source = 0;
  int ntfs = tfs->size();
  
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    const string& pwm_source = tf.getPWMSource();
    if (pwm_source != string(""))
      print_source = max(print_source, (int) pwm_source.size());
  }
  
  int p = mode->getPrecision();
  int w = p+7;
  os << setprecision(p)
     << setw(w) << "name"
     << setw(w) << "index"
     << setw(w) << "index_f"
     << setw(w) << "index_r"
     << setw(w) << "start"
     << setw(w) << "end"
     << setw(w) << "score"
     << setw(w) << "K"
     << setw(w) << "K*kmax"
     << setw(w) << "strand";
     
  if (print_source)
    os << setw(print_source+2) << "source";
  
  os << endl;
  printSites(gene, tf, os, p, print_source);
}
  
  
void Bindings::printSites(Gene& gene, TF& tf, ostream& os, int p, int print_source)
{ 
  int w = p+7;
  site_ptr_vector& tmp_sites = sites[&gene][&tf];
  double maxscore   = tf.getMaxScore();
  int    left_bound = gene.getLeftBound();
  int    nsites     = tmp_sites.size();
  for (int i=0; i<nsites; i++)
  {
    BindingSite& b = *tmp_sites[i];
    os << setprecision(p)
       << setw(w) << b.tf->getName()
       << setw(w) << b.index_in_site_map
       << setw(w) << b.index_in_ordered_f
       << setw(w) << b.index_in_ordered_r
       << setw(w) << b.m + left_bound
       << setw(w) << b.n + left_bound
       << setw(w) << b.score
       << setw(w) << b.K_exp_part
       << setw(w) << b.K_exp_part_times_kmax;
     if (b.orientation == 'F')
       os << setw(w) << "+";
     else 
       os << setw(w) << "-";
     
     if (print_source)
        os << setw(print_source+2) << tmp_sites[i]->tf->getPWMSource();
      os << endl;
  }
}

void Bindings::printScores(Gene& gene, ostream& os)
{
  int p = mode->getPrecision();
  int w = p+7;
  
  map<TF*, TFscore> gscores = scores[&gene];
  
  int ntfs = tfs->size();
  
  os << setw(w) << "base";
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    os << setw(w) << tf.getName() + "_f";
    os << setw(w) << tf.getName() + "_r";
  }
  os << endl;
  
  int length     = gene.length();
  int left_bound = gene.getLeftBound();
  
  for (int bp=0; bp<length; bp++)
  {
    os << setw(w) << bp + left_bound;
    for (int i=0; i<ntfs; i++)
    {
      TF& tf = tfs->getTF(i);
      os << setw(w) << gscores[&tf].fscore[bp];
      os << setw(w) << gscores[&tf].rscore[bp];
    }
    os << endl;
  }
}
      

void Bindings::printTotalOccupancy(Gene& gene, ostream& os, bool invert)
{
  int p = mode->getPrecision();
  int w = p + 7;
  os << setprecision(p);
  
  int ntfs = tfs->size();
  if (invert == false)
  {
    // loop through everything and figure out the header
    os << setw(w) << "id";
  
    // find all the ids, using map to prevent duplicates
    // simultaneously print binding site headers
    for (int i = 0; i<ntfs; i++)
    {
      TF& tf = tfs->getTF(i);
      site_ptr_vector& tmp_sites = sites[&gene][&tf];
      int nsites = tmp_sites.size();
      for (int i=0; i<nsites; i++)
      {
        string header = tf.getName();
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
      string& id = IDs[nuc];
      os << setw(w) << id;
    
      for (int i = 0; i<ntfs; i++)
      {
        TF& tf = tfs->getTF(i);
        site_ptr_vector& tmp_sites = sites[&gene][&tf];
        int nsites = tmp_sites.size();
        for (int i=0; i<nsites; i++)
          os << setw(w) << tmp_sites[i]->total_occupancy[nuc];
      }
      os << endl;
    }
  } 
  else 
  {
    os << setw(w) << "site";
    for (int nuc=0; nuc<nnuc; nuc++)
    {
      string& id = IDs[nuc];
      os << setw(w) << id;
    }
    os << endl;
    for (int i = 0; i<ntfs; i++)
    {
      TF& tf = tfs->getTF(i);
      site_ptr_vector& tmp_sites = sites[&gene][&tf];
      int nsites = tmp_sites.size();
      for (int i=0; i<nsites; i++)
      {
        string header = tf.getName();
        stringstream convert;
        convert << i;
        string num = convert.str();
        header += num;
        //header << i;
        os << setw(w) <<  header;
        for (int nuc=0; nuc<nnuc; nuc++)
        {
          os << setw(w) << tmp_sites[i]->total_occupancy[nuc];
        }
      os << endl;
      }
    }
  }
}

void Bindings::printEffectiveOccupancy(Gene& gene, ostream& os, bool invert)
{
  int w = 12;
  os << setprecision(3);
  
  int ntfs = tfs->size();
  int ngenes = genes->size();
  if (!invert)
  {
    
    // loop through everything and figure out the header
    os << setw(w) << "id";
  
    // find all the ids, using map to prevent duplicates
    // simultaneously print binding site headers
    for (int k=0; k<ntfs; k++)
    {
      TF& tf = tfs->getTF(k);
      site_ptr_vector& tmp_sites = sites[&gene][&tf];
      int nsites = tmp_sites.size();
      int nmodes = tf.getNumModes();
      for (int i=0; i<nsites; i++)
      {      
        for (int j=0; j<nmodes; j++)
        {
          string header = tf.getName();
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
      string& id = IDs[nuc];
      os << setw(w) << id;
   
      for (int k=0; k<ntfs; k++)
      {
        TF& tf = tfs->getTF(k);
        site_ptr_vector& tmp_sites = sites[&gene][&tf];
        int nsites = tmp_sites.size();
        int nmodes = tf.getNumModes();
        for (int i=0; i<nsites; i++)
        {
          for (int j=0; j<nmodes; j++)
            os << setw(w) << tmp_sites[i]->effective_occupancy[j][nuc];
        }
      }
      os << endl;
    }
  }
  else 
  {
    os << setw(w) << "site";
    for (int nuc=0; nuc<nnuc; nuc++)
    {
      string& id = IDs[nuc];
      os << setw(w) << id;
    }
    os << endl;

    for (int k=0; k<ntfs; k++)
    {
      TF& tf = tfs->getTF(k);
      site_ptr_vector& tmp_sites = sites[&gene][&tf];
      int nsites = tmp_sites.size();
      int nmodes = tf.getNumModes();
      for (int i=0; i<nsites; i++)
      {
        for (int j=0; j<nmodes; j++)
        {
          string header = tf.getName();
          stringstream convert;
          convert << i;
          if (j>0) convert << "_" << j+1;
          string num = convert.str();
          header += num;
          //header << i;
          os << setw(w) <<  header;
          for (int nuc=0; nuc<nnuc; nuc++)
          {
            os << setw(w) << tmp_sites[i]->effective_occupancy[j][nuc];
          }
        os << endl;
        }
      }
    }
  }
}

void Bindings::printModeOccupancy(Gene& gene, ostream& os, bool invert)
{
  int w = 12;
  os << setprecision(3);
  int ntfs = tfs->size();
  int ngenes = genes->size();
  if (!invert)
  {
    
    // loop through everything and figure out the header
    os << setw(w) << "id";
  
    // find all the ids, using map to prevent duplicates
    // simultaneously print binding site headers
    for (int k=0; k<ntfs; k++)
    {
      TF& tf = tfs->getTF(k);
      site_ptr_vector& tmp_sites = sites[&gene][&tf];
      int nsites = tmp_sites.size();
      int nmodes = tf.getNumModes();
      for (int i=0; i<nsites; i++)
      {      
        for (int j=0; j<nmodes; j++)
        {
          string header = tf.getName();
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
      string& id = IDs[nuc];
      os << setw(w) << id;
   
      for (int k=0; k<ntfs; k++)
      {
        TF& tf = tfs->getTF(k);
        site_ptr_vector& tmp_sites = sites[&gene][&tf];
        int nsites = tmp_sites.size();
        int nmodes = tf.getNumModes();
        for (int i=0; i<nsites; i++)
        {
          for (int j=0; j<nmodes; j++)
            os << setw(w) << tmp_sites[i]->mode_occupancy[j][nuc];
        }
      }
      os << endl;
    }
  }
  else 
  {
    os << setw(w) << "site";
    for (int nuc=0; nuc<nnuc; nuc++)
    {
      string& id = IDs[nuc];
      os << setw(w) << id;
    }
    os << endl;

    for (int k=0; k<ntfs; k++)
    {
      TF& tf = tfs->getTF(k);
      site_ptr_vector& tmp_sites = sites[&gene][&tf];
      int nsites = tmp_sites.size();
      int nmodes = tf.getNumModes();
      for (int i=0; i<nsites; i++)
      {
        for (int j=0; j<nmodes; j++)
        {
          string header = tf.getName();
          stringstream convert;
          convert << i;
          if (j>0) convert << "_" << j+1;
          string num = convert.str();
          header += num;
          //header << i;
          os << setw(w) <<  header;
          for (int nuc=0; nuc<nnuc; nuc++)
          {
            os << setw(w) << tmp_sites[i]->mode_occupancy[j][nuc];
          }
        os << endl;
        }
      }
    }
  }
}
