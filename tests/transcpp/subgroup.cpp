/*********************************************************************************
*                                                                                *
*     subgroup.cpp                                                               *
*                                                                                *
*     Contains the subgroup classes and methods.                                 *
*                                                                                *
*********************************************************************************/

#include "subgroup.h"
#include <boost/foreach.hpp>
#include <algorithm>
#include <limits>
#include <math.h>
#include <climits>
#include <boost/unordered_map.hpp>

#define foreach_ BOOST_FOREACH

/********************************   Subgroup    *********************************/


/*
bool compareBindingSiteRight(BindingSite* a, BindingSite* b)
{
  return(a->n < b->n);
}

bool compareBindingSiteLeft(BindingSite* a, BindingSite* b)
{
  return(a->m > b->m);
}*/


/*    Constructors    */

Subgroup::Subgroup() {}

Subgroup::Subgroup(BindingSite* site, bindings_ptr b) 
{
  bindings  = b;
  addSite(site);
}


/*    Setters   */

void Subgroup::addSite(BindingSite* site)
{
  sites_f.push_back(site);
  if ( sites_f.size()==1 )
  {
    right_bound = site->n;
    left_bound  = site->m;
  } else
  {
    left_bound   = min(left_bound,  site->m);
    right_bound  = max(right_bound, site->n);
  }
}

void Subgroup::addSubgroup(Subgroup& sub)
{
  vector<BindingSite*>& subSites = sub.getSites();
  int nsites = subSites.size();
  for (int i=0; i<nsites; i++)
  {
    addSite( subSites[i] );
  }
}


/*    Getters   */

int Subgroup::size() const
{
  return sites_f.size();
}

vector<BindingSite*>& Subgroup::getSites()
{
  return sites_f;
}

int Subgroup::getRightBound() const { return right_bound; }

int Subgroup::getLeftBound() const { return left_bound; }


/*    Methods   */

bool Subgroup::overlaps(BindingSite* site)
{
  return ( site->m < right_bound && left_bound < site->n);
}

// check to see if this site coops with any sites already in group
bool Subgroup::checkCoop(BindingSite* site1)
{
  char o1,o2;
  TF* tf1 = site1->tf;
  int nsites = sites_f.size();
  for (int i=0; i<nsites; i++)
  {
    BindingSite* site2 = sites_f[i];
    TF* tf2 = site2->tf;
    
    if (site1->n <= site2->m)
    {
      o1 = site1->orientation;
      o2 = site2->orientation;
    }
    else
    {
      o2 = site1->orientation;
      o1 = site2->orientation;
    }
    
    if (tf1->checkCoops(tf2, o1, o2))
    {
      coop_ptr cur_coop = tf1->getCoop(tf2);
      int dist1 = abs(site2->m - site1->n);
      int dist2 = abs(site1->m - site2->n);
      int dist  = min(dist1, dist2);
      double distfunc = cur_coop->distFunc(dist);
      if (distfunc > 0)
        return true;
    }
  }
  return false;
}

bool Subgroup::overlaps(BindingSite* site1, BindingSite* site2)
{
  return ( site1->m < site2->n && site2->m < site1->n);
}

void Subgroup::sort()
{
  std::sort(sites_f.begin(), sites_f.end(), compareBindingSiteRight);
  int nsites = sites_f.size();
  sites_r.resize(nsites);
  f2r.resize(nsites);
  for (int i=0; i<nsites; i++) 
    sites_r[i] = sites_f[nsites - i - 1];
  std::sort(sites_r.begin(), sites_r.end(), compareBindingSiteLeft);

  for (int i=0; i<nsites; i++)
  {
    for (int j=0; j<nsites; j++)
    {
      if (sites_f[i] == sites_r[j])
      {
        f2r[i] = j;
        continue;
      }
    }
  }
}

// for each site, find the last site that does not compete
void Subgroup::pre_process()
{
  int nsites = sites_f.size();
  int nnuc = bindings->getNnuc();
  
  ZF.resize(nsites+1);
  ZR.resize(nsites+1);
    
  ZF[0].Z.resize(nnuc);
  ZF[0].Zc.resize(nnuc);
  ZF[0].Znc.resize(nnuc);
  
  ZR[0].Z.resize(nnuc);
  ZR[0].Zc.resize(nnuc);
  ZR[0].Znc.resize(nnuc);
  
  for (int i=0; i<nnuc; i++)
  {
    ZF[0].Z[i] = 1.0;
    ZR[0].Z[i] = 1.0;
  }
    
  
  for (int i=0; i<nsites; i++) // loop over all sites
  {
    int pindex = i + 1;
    
    ZF[pindex].Z.resize(nnuc);
    ZF[pindex].Zc.resize(nnuc, 0);
    ZF[pindex].Znc.resize(nnuc, 0);
    
    ZR[pindex].Z.resize(nnuc);
    ZR[pindex].Zc.resize(nnuc, 0);
    ZR[pindex].Znc.resize(nnuc, 0);
      
    BindingSite* s1f = sites_f[i];
    BindingSite* s1r = sites_r[i];
    
    TF* tf1f = s1f->tf;
    TF* tf1r = s1r->tf;
    
    char o1f = s1f->orientation;
    char o1r = s1r->orientation;
    
    for (int j = 0; j<i; j++) // loop over previous sites
    {
      BindingSite* s2f = sites_f[j];
      BindingSite* s2r = sites_r[j];
      
      TF* tf2f = s2f->tf;
      TF* tf2r = s2r->tf;
      
      char o2f = s2f->orientation;
      char o2r = s2r->orientation;
      
      pre_process_pair(ZF, i, s1f, tf1f, o1f, j, s2f, tf2f, o2f); 
      pre_process_pair(ZR, i, s1r, tf1r, o1r, j, s2r, tf2r, o2r); 
    }
  }
}

void Subgroup::pre_process_pair(vector<Partition>& p, int i, BindingSite* s1, TF* tf1, char o1, int j, BindingSite* s2, TF* tf2, char o2)
{
  int pindex = i + 1;
  if (!overlaps(s1,s2))
  {
    p[pindex].last = j + 1;
    
    if (tf1->checkCoops(tf2, o2, o1))
    {
      coop_ptr cur_coop = tf1->getCoop(tf2);
      int dist1 = abs(s1->m - s2->n);
      int dist2 = abs(s2->m - s1->n);
      int dist = min(dist1,dist2);
      
      double distfunc = cur_coop->distFunc(dist);
      if (distfunc > 0)
      {
        p[pindex].coops = true;
        p[pindex].coop_site.push_back(j);
        p[pindex].coop_past.push_back(p[j+1].last);
        p[pindex].dist_coef.push_back(distfunc);
        p[pindex].coop.push_back(cur_coop);
      }
    }
  }
}
  

void Subgroup::
iterate_partition(vector<Partition>& p, vector<BindingSite*>& sites, int site_index)
{
  int pindex           = site_index + 1;

  Partition& cur_part  = p[pindex];
  Partition& last_part = p[cur_part.last];
  Partition& init_part = p[site_index];
    
  int ncoops = cur_part.coop_site.size();
    
  vector<double>& kv   = sites[site_index]->kv;
  
  int nnuc = bindings->getNnuc();
  
  if (ncoops)
  {
    for (int i=0; i<nnuc; i++)
    {
      double cur_kv = kv[i];
      double init_Z = init_part.Z[i];
      double last_Z = last_part.Z[i];
      double new_Z  = last_Z*cur_kv;
      
      cur_part.Z[i]   = init_Z + new_Z;
      cur_part.Znc[i] = new_Z;
      cur_part.Zc[i]  = 0;
    }
    
    for (int j=0; j<ncoops; j++)
    {
      double kcoop     = cur_part.coop[j]->getK();
      int    coop_site = cur_part.coop_site[j]; // the site it coops with
      int    coop_past = cur_part.coop_past[j];
      double dfunk     = cur_part.dist_coef[j];
      kcoop *= dfunk;
      
      vector<double>& coopkv = sites[coop_site]->kv;
      
      for (int k=0; k<nnuc; k++)
      {
        double cur_coopkv = coopkv[k];
        double cur_kv     = kv[k];
        double last_Z     = p[coop_past].Z[k];
        double weight     = last_Z*cur_coopkv*cur_kv*kcoop;
        cur_part.Z[k]  += weight;
        cur_part.Zc[k] += weight;
      }
    }
  }
  else
  {
    for (int i=0; i<nnuc; i++)
    {
      double cur_kv = kv[i];
      double init_Z = init_part.Z[i];
      double last_Z = last_part.Z[i];
      double new_Z  = last_Z*cur_kv;
      
      cur_part.Z[i]   = init_Z + new_Z;
      cur_part.Znc[i] = new_Z;
      //cur_part.Zc[i]  = 0;
    }
  }
      
}

     
void Subgroup::occupancy()
{
  int nsites = sites_f.size();
  int nnuc   = bindings->getNnuc();
  
  for (int i=0; i<nsites; i++)
  {
    iterate_partition(ZF, sites_f, i); // Z0 now holds all the partitions
    iterate_partition(ZR, sites_r, i);
  }
   
  vector<double>& Z = ZF[nsites].Z;

  for (int i=0; i<nsites; i++)
  { 
    int r_idx = f2r[i];
    
    vector<double>& zfnc = ZF[i+1].Znc;
    vector<double>& zrnc = ZR[r_idx+1].Znc;
    
    vector<double>& zfc = ZF[i+1].Zc;
    vector<double>& zrc = ZR[r_idx+1].Zc;
      
    BindingSite* site = sites_f[i];
    
    for (int j=0; j<nnuc; j++)
    {
      
      double kv = site->kv[j];
      if ( kv == 0)
      {
        site->total_occupancy[j]   = 0;
        site->mode_occupancy[0][j] = 0;
      }
      else
      {
        double denom = (zfnc[j] + zfc[j])*(zrnc[j] + zrc[j]) - zrc[j]*zfc[j];
        double f  = denom/(Z[j]*kv);
        site->total_occupancy[j]   = f;
        site->mode_occupancy[0][j] = f;
      }
    }
  }
}
      

/*    Print   */

void Subgroup::print(ostream& os)
{
  os << setw(10) << "f_index";
  os << setw(10) << "r_index";
  
  printSiteHeader(os);
  
  os << setw(10) << "last_f" << setw(10) << "last_r" << endl;
  
  int nsites = sites_f.size();
  for (int i=0; i<nsites; i++)
  {
    os << setw(10) << i;
    os << setw(10) << f2r[i];
    
    printSite(*sites_f[i],os);
    os << setprecision(3) << setw(10);
    os << ZF[i+1].last;
    os << setw(10) << ZR[f2r[i]+1].last << endl;
  }
}


/******************************   Subgroups   **********************************/

Subgroups::Subgroups() {}

void Subgroups::clear()
{
  groups.clear();
}

void Subgroups::clear(Gene& gene)
{
  groups[&gene].clear();
}


void Subgroups::create(genes_ptr g, tfs_ptr t, bindings_ptr b, mode_ptr m)
{
  genes     = g;
  tfs       = t;
  bindings  = b;
  mode      = m;
  
  update();
}

void Subgroups::addSites(Gene& gene)
{
  boost::unordered_map<TF*, site_ptr_vector>& gene_sites = bindings->getSites(gene);
  list<Subgroup>& gene_groups = groups[&gene];
  
  int ntfs = tfs->size();
  for(int i=0; i<ntfs; i++)
  {  
    TF& tf = tfs->getTF(i);
    site_ptr_vector& sites = gene_sites[&tf];
    int nsites = sites.size();
    for (int j=0; j<nsites; j++)
      addSite(gene_groups, sites[j]);
  }
  
  list<Subgroup>::iterator i;
  for (i=gene_groups.begin(); i != gene_groups.end(); ++i)
  {
    i->sort();
    i->pre_process();
  }
}

void Subgroups::addSite(list<Subgroup>& gene_groups, site_ptr site)
{
  list<Subgroup>::iterator i;  // iterate through list
  list<Subgroup>::iterator j;  // point to group that site was added to
   
  bool added = false;  
  i = gene_groups.begin();
  while (i != gene_groups.end())
  {
    /* I check every subgroup. I add the site to the first one that it overlaps or
    coops with. If it coops with additional ones I merge them. */
    bool condition = i->overlaps(site.get()) || i->checkCoop(site.get());
    
    if (condition && added==false)
    {
      i->addSite(site.get());
      added = true;
      j = i;
      i++;
    } else if ( condition && added==true)
    {
      j->addSubgroup(*i);
      i = gene_groups.erase(i);
    } else {
      i++;
    }
  }
  
  if (added==false)
    gene_groups.push_back(Subgroup(site.get(), bindings));
}
    
  

void Subgroups::update()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    update(gene);
  }
}

void Subgroups::update(Gene& gene)
{
  groups[&gene].clear();
  addSites(gene);
}
 

void Subgroups::calc_f()
{
  int ngenes = genes->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    list<Subgroup>& gene_groups = groups[&gene];
    list<Subgroup>::iterator j;
    for (j=gene_groups.begin(); j != gene_groups.end(); ++j)
      j->occupancy();
  }
}

void Subgroups::calc_f(Gene& gene)
{
  list<Subgroup>& gene_groups = groups[&gene];
  list<Subgroup>::iterator i;
  for (i=gene_groups.begin(); i != gene_groups.end(); ++i)
    i->occupancy();
}


void Subgroups::save()
{
  saved_groups = groups;
}

void Subgroups::save(Gene& gene)
{
  saved_groups[&gene] = groups[&gene];
}

void Subgroups::restore()
{
  groups = saved_groups;
}

void Subgroups::restore(Gene& gene)
{
  groups[&gene] = saved_groups[&gene];
}

void Subgroups::print(ostream& os)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    os << gene.getName();
    print(gene, os);
  }
}

void Subgroups::print(Gene& gene, ostream& os)
{
  list<Subgroup>& gene_groups = groups[&gene];
  int counter = 1;
  list<Subgroup>::iterator i;
  
  for (i=gene_groups.begin(); i != gene_groups.end(); ++i)
  {
    os << "Subroup " << counter++ << endl;
    os << "Left bound:  " << i->getLeftBound()  << endl;
    os << "Right bound: " << i->getRightBound() << endl;
    i->print(os);
  }
}
