/*********************************************************************************
*                                                                                *
*     quenching.cpp                                                              *
*                                                                                *
*     Contains structure and methods to hold quenching interactions              *
*                                                                                *
*********************************************************************************/

#include "quenching.h"
#include <boost/foreach.hpp>
#include <limits>

# define foreach_ BOOST_FOREACH

QuenchingInteractions::QuenchingInteractions() {}

void QuenchingInteractions
::create(genes_ptr g, tfs_ptr t, bindings_ptr b, distances_ptr d)
{
  genes     = g;
  tfs       = t;
  bindings  = b;
  distances = d;
  
  dist = distances->getDistance("Quenching");
  
  int ngenes = genes->size();
  for (int k=0; k<ngenes; k++)
  {
    Gene& gene = genes->getGene(k);
    int ntfs   = tfs->size();
    for (int i = 0; i<ntfs; i++) // loop through actors
    {
      TF& tf1 = tfs->getTF(i);
      if (tf1.neverQuenches()) continue; // skip if never a quencher
      for (int j=0; j<ntfs; j++) // loop through targets
      {
        TF& tf2 = tfs->getTF(j);
        if (tf2.neverActivates()) continue; // skip if never an activator
        set(gene, tf1, tf2);
      }
    }   
  }
  //printSummary();
}

void QuenchingInteractions
::set(Gene& gene, TF& actor, TF& target)
{
  site_ptr_vector& actorsites  = bindings->getSites(gene, actor);
  site_ptr_vector& targetsites = bindings->getSites(gene, target);
  
  vector<QuenchingInteraction>& q = quenches[&gene][&actor][&target];
  double max_dist     = dist->getMaxDistance();
  int    nactors      = actorsites.size();
  int    ntargets     = targetsites.size();
  
  for (int i=0; i<nactors; i++)
  {
    site_ptr actor_ptr  = actorsites[i];
    int m1 = actor_ptr->m;
    int n1 = actor_ptr->n;
    
    for (int j=0; j<ntargets; j++)
    {
      
      site_ptr target_ptr = targetsites[j];
      
      int m2 = target_ptr->m;
      int n2 = target_ptr->n;
      double d;
      
      int dm = abs(m1 - n2);
      int dn = abs(n1 - m2);
      
      if (dn <= dm)
        d = dn;
      else
        d = dm;
      
      bool  overlapped = (m1 < n2 && m2 < n1);
      double df        = dist->getDistFunc(d);
      if ( d < max_dist && !overlapped && df > 0)
      {
        QuenchingInteraction quench;
        quench.actor  = actor_ptr;
        quench.target = target_ptr;
        quench.distcoef = df;
        
        q.push_back(quench);
        
      }
    }
  }
}
  
bool QuenchingInteractions
::hasQuenchingInteractions(Gene& gene, TF& actor, TF& target)
{
  if (quenches.find(&gene) == quenches.end())
    return false;
  else
  {
    if (quenches[&gene].find(&actor) == quenches[&gene].end())
      return false;
    else
    {
      if (quenches[&gene][&actor].find(&target) == quenches[&gene][&actor].end())
        return false;
    }
  }
  return true;
}


void QuenchingInteractions
::calc()
{
  int ngenes = genes->size();
  for (int k=0; k<ngenes; k++)
  {
    Gene& gene = genes->getGene(k);
    calc(gene);
  }
}

void QuenchingInteractions
::calc(Gene& gene)
{
  initialize(gene);
    
  int ntfs   = tfs->size();
  for (int i = 0; i<ntfs; i++) // loop through actors
  {
    TF& tf1 = tfs->getTF(i);
    if (tf1.neverQuenches()) continue; // skip if never a quencher
    for (int j=0; j<ntfs; j++) // loop through targets
    {
      TF& tf2 = tfs->getTF(j);
      if (tf2.neverActivates()) continue; // skip if never an activator
      calc(gene, tf1, tf2);
    }
  } 
}
  
void QuenchingInteractions
::calc(Gene& gene, TF& actor, TF& target)
{
  vector<QuenchingInteraction>& quench_vector = quenches[&gene][&actor][&target];
  
  vector<double> actor_coefs  = actor.getCoefs();
  vector<double> target_coefs = target.getCoefs();

  int n_actor_modes  = actor_coefs.size();
  int n_target_modes = target_coefs.size();
  
  int nquench = quench_vector.size();
  
  for (int i=0; i<nquench; i++)
  {
    QuenchingInteraction& quench      = quench_vector[i];
    BindingSite& actor_site  = *quench.actor;
    BindingSite& target_site = *quench.target;
    double distcoef          = quench.distcoef;

    vector<double>* actor_occupancy;
    vector<double>* target_occupancy;
    
    for (int j=0; j<n_actor_modes; j++)
    {
      double efficiency = -actor_coefs[j];
      if (efficiency <= 0) continue;
      double efd = efficiency * distcoef;
      actor_occupancy = &(actor_site.mode_occupancy[j]);
      for (int k=0; k<n_target_modes; k++)
      {
        if (target_coefs[k] <= 0) continue;
        
        target_occupancy = &(target_site.effective_occupancy[k]);
        quench_f(*actor_occupancy, *target_occupancy, efd);
      }
    }
  }
}


void QuenchingInteractions
::quench_f(vector<double>& actor_vec, vector<double>& target_vec, double efd)
{
  int n = actor_vec.size();
  for (int i=0; i<n; i++)
  {
    double actor_occupancy = actor_vec[i];
    double reduction = 1 - actor_occupancy * efd;
    double* current = &target_vec[i];
    //target_vec[i] *= reduction;
    *current = *current * reduction;
  }
}

void QuenchingInteractions
::save()
{
  saved_quenches = quenches;
}

void QuenchingInteractions
::save(Gene& gene)
{
  saved_quenches[&gene] = quenches[&gene];
}

void QuenchingInteractions
::clear()
{
  quenches.clear();
}

void QuenchingInteractions
::clear(Gene& gene)
{
  quenches[&gene].clear();
}

void QuenchingInteractions
::restore()
{
  quenches = saved_quenches;
}

void QuenchingInteractions
::restore(Gene& gene)
{
  quenches[&gene] = saved_quenches[&gene];
}

void QuenchingInteractions
::update()
{
  int ngenes = genes->size();
 
  for (int k=0; k<ngenes; k++)
  {
    Gene& gene = genes->getGene(k);
    update(gene);
  }
}

void QuenchingInteractions
::update(Gene& gene)
{
  int ntfs   = tfs->size();

  for (int i = 0; i<ntfs; i++) // loop through actors
  {
    TF& tf1 = tfs->getTF(i);
    if (tf1.neverQuenches()) continue; // skip if not a quencher
    for (int j=0; j<ntfs; j++) // loop through targets
    {
      TF& tf2 = tfs->getTF(j);
      if (tf2.neverActivates()) continue; // skip if never an activator
      quenches[&gene][&tf1][&tf2].clear();
      set(gene, tf1, tf2);
    }
  } 
}
  

void QuenchingInteractions
::printSummary()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  int total  = 0;
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    cerr << gene.getName() << endl;
    for (int j=0; j<ntfs; j++)
    {
      TF& tf1 = tfs->getTF(j);
      for (int k=0; k<ntfs; k++)
      {
        TF& tf2 = tfs->getTF(k);
        int nq = quenches[&gene][&tf1][&tf2].size();
        cerr << tf1.getName() << " -> " << tf2.getName() << " : " << nq << endl;
        total += nq;
      }
    }
  }
  cerr << "Total: " << total << endl << endl;
}

void QuenchingInteractions
::initialize()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    initialize(gene);
  }
}

void QuenchingInteractions
::initialize(Gene& gene)
{
  int ntfs   = tfs->size();
  for (int j = 0; j<ntfs; j++) 
  {
    TF& tf = tfs->getTF(j);
    initialize(gene, tf);
  }
}
      

void QuenchingInteractions
::initialize(Gene& gene, TF& tf)
{
  site_ptr_vector& sites = bindings->getSites(gene, tf);

  vector<double> coefs = tf.getCoefs();
  
  int nmodes = coefs.size();
  int nsites = sites.size();

  for (int i=0; i<nsites; i++)
  {
    BindingSite* site = sites[i].get();
    int nnuc = site->total_occupancy.size();
    for (int j=0; j<nmodes; j++)
    {
      if (coefs[j] > 0)
        site->effective_occupancy[j] = site->mode_occupancy[j];
      else
      {
        for (int k=0; k<nnuc; k++)
          site->effective_occupancy[j][k] = 0;
      }
    }
  }
}




ModifyingInteractions::ModifyingInteractions() {}


void ModifyingInteractions
::create(genes_ptr g, tfs_ptr t, bindings_ptr b, distances_ptr d)
{
  genes     = g;
  tfs       = t;
  bindings  = b;
  distances = d;
  
  dist = distances->getDistance("Quenching");

  int ngenes = genes->size();
  for (int k=0; k<ngenes; k++)
  {
    Gene& gene = genes->getGene(k);
    int ntfs   = tfs->size();
    for (int i = 0; i<ntfs; i++) // loop through actors
    {
      TF& tf1 = tfs->getTF(i);
      vector< pair<TF*, coeffect_ptr> >& targets = tf1.getTargets();
      int ntargets = targets.size();
      for (int j=0; j<ntargets; j++)
      {
        pair<TF*, coeffect_ptr>& tmp_pair = targets[j];
        TF&          tf2      = *(tmp_pair.first);
        coeffect_ptr cur_coef = tmp_pair.second;
        set(gene, tf1, tf2, cur_coef);
      }
    }   
  }
}

void ModifyingInteractions
::set(Gene& gene, TF& actor, TF& target, coeffect_ptr coef)
{
  site_ptr_vector& actorsites  = bindings->getSites(gene, actor);
  site_ptr_vector& targetsites = bindings->getSites(gene, target);

  vector<ModifyingInteraction>& v = mods[&gene][&actor][&target];
  
  double max_dist     = coef->getMaxDistance();
  int    nactors      = actorsites.size();
  int    ntargets     = targetsites.size();

  for (int i=0; i<nactors; i++)
  {
    for (int j=0; j<ntargets; j++)
    {
      site_ptr actor_ptr  = actorsites[i];
      site_ptr target_ptr = targetsites[j];
      
      int m1 = actor_ptr->m;
      int n1 = actor_ptr->n;
      int m2 = target_ptr->m;  
      int n2 = target_ptr->n;  
      double d;
      
      int dm = abs(m1 - n2);
      int dn = abs(n1 - m2);
      
      if (dn <= dm)
        d = dn;
      else
        d = dm;
      
      bool  overlapped = (m1 < n2 && m2 < n1);
      double df        = coef->distFunc(d);
      if ( d < max_dist && !overlapped && df > 0)
      {
        ModifyingInteraction mod;
        mod.actor  = actor_ptr;
        mod.target = target_ptr;
        mod.distcoef = df;
        mod.coef = coef;
        
        v.push_back(mod);
      }
    }
  }
}
        
void ModifyingInteractions
::calc()
{
  initialize();
  int ngenes = genes->size();
  for (int k=0; k<ngenes; k++)
  {
    Gene& gene = genes->getGene(k);
    int ntfs   = tfs->size();
    for (int i = 0; i<ntfs; i++) // loop through actors
    {
      TF& tf1 = tfs->getTF(i);
      vector< pair<TF*, coeffect_ptr> >& targets = tf1.getTargets();
      int ntargets = targets.size();
      for (int j=0; j<ntargets; j++)
      {
        pair<TF*, coeffect_ptr>& tmp_pair = targets[j];
        TF&          tf2      = *(tmp_pair.first);
        coeffect_ptr cur_coef = tmp_pair.second;
        calc(gene, tf1, tf2, cur_coef);
      }
    }   
  }
}

void ModifyingInteractions
::calc(Gene& gene)
{
  initialize(gene);

  int ntfs   = tfs->size();
  for (int i = 0; i<ntfs; i++) // loop through actors
  {
    TF& tf1 = tfs->getTF(i);
    vector< pair<TF*, coeffect_ptr> >& targets = tf1.getTargets();
    int ntargets = targets.size();
    for (int j=0; j<ntargets; j++)
    {
      pair<TF*, coeffect_ptr>& tmp_pair = targets[j];
      TF&          tf2      = *(tmp_pair.first);
      coeffect_ptr cur_coef = tmp_pair.second;
      calc(gene, tf1, tf2, cur_coef);
    }
  }   

}

void ModifyingInteractions
::calc(Gene& gene, TF& actor, TF& target, coeffect_ptr coef)
{
  vector<ModifyingInteraction>& mod_vector = mods[&gene][&actor][&target];
  
  double efficiency = coef->getEfficiency();
  int    coef_idx   = coef->getIdx() - 1;
  
  int nmods = mod_vector.size();
  
  for (int i=0; i<nmods; i++)
  {

    ModifyingInteraction& mod = mod_vector[i];
    BindingSite& actor_site   = *mod.actor;
    BindingSite& target_site  = *mod.target;
    double distcoef           = mod.distcoef;
    
    vector<double>& actor_occupancy = actor_site.total_occupancy;
    vector<double>& start_occupancy = target_site.mode_occupancy[0];
    vector<double>& end_occupancy   = target_site.mode_occupancy[coef_idx];
    
    mod_f(actor_occupancy, start_occupancy, end_occupancy, efficiency, distcoef);
  } 
}

void ModifyingInteractions
::mod_f(vector<double>& actor_vec, 
        vector<double>& start_vec,
        vector<double>& end_vec,
        double ef, double d)
{
  int n = actor_vec.size();

  for (int i=0; i<n; i++)
  {
    double actor_occupancy = actor_vec[i];
    double reduction = actor_occupancy * ef * d;
    double current   = start_vec[i];
    double reduced   = current*(1-reduction);
    double increased = current - reduced;
    start_vec[i]     = reduced;
    end_vec[i]      += increased;
  }
}
    
void ModifyingInteractions
::initialize()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    initialize(gene);
  }
}

void ModifyingInteractions
::initialize(Gene& gene)
{
  int ntfs   = tfs->size();
  for (int j = 0; j<ntfs; j++) 
  {
    TF& tf = tfs->getTF(j);
    initialize(gene, tf);
  }
}
      

void ModifyingInteractions
::initialize(Gene& gene, TF& tf)
{
  int nmodes             = tf.getNumModes();
  site_ptr_vector& sites = bindings->getSites(gene, tf);

  int nsites = sites.size();

  for (int i=0; i<nsites; i++)
  {
    BindingSite* site = sites[i].get();
    vector<double>& total_occupancy = site->total_occupancy;
    site->mode_occupancy[0] = total_occupancy;
    int n = total_occupancy.size();
    for (int j=1; j<nmodes; j++)
    {
      vector<double>& tmp_occ = site->mode_occupancy[j];
      for (int k=0; k<n; k++)
        tmp_occ[k] = 0;
    }
  }
}




void ModifyingInteractions
::save(Gene& gene)
{
  saved_mods[&gene] = mods[&gene];
}

void ModifyingInteractions
::save()
{
  saved_mods = mods;
}

void ModifyingInteractions
::clear(Gene& gene)
{
  mods[&gene].clear();
}

void ModifyingInteractions
::clear()
{
  mods.clear();
}

void ModifyingInteractions
::restore(Gene& gene)
{
  mods[&gene] = saved_mods[&gene];
}

void ModifyingInteractions
::restore()
{
  mods = saved_mods;
}


void ModifyingInteractions
::update()
{
  int ngenes = genes->size();
  for (int k=0; k<ngenes; k++)
  {
    Gene& gene = genes->getGene(k);
    update(gene);
  }
}

void ModifyingInteractions
::update(Gene& gene)
{
  int ntfs   = tfs->size();
  for (int i = 0; i<ntfs; i++) // loop through actors
  {
    TF& tf1 = tfs->getTF(i);
    vector< pair<TF*, coeffect_ptr> >& targets = tf1.getTargets();
    int ntargets = targets.size();
    for (int j=0; j<ntargets; j++)
    {
      pair<TF*, coeffect_ptr>& tmp_pair = targets[j];
      TF&          tf2      = *(tmp_pair.first);
      coeffect_ptr cur_coef = tmp_pair.second;
      mods[&gene][&tf1][&tf2].clear();
      set(gene, tf1, tf2, cur_coef);
    }
  }   
}










