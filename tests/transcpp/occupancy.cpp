/*********************************************************************************
*                                                                                *
*     occupancy.cpp                                                              *
*                                                                                *
*     It has become increasingly clear over the course of writing this code      *
*     that most of the intense calculations are much more efficient if they      *
*     can simply traverse an array of occupancies. This class is designed        *
*     to be a matrix of occupancies and effective occupancies and all subgroup   *
*     and quenching calculations will take a pointer to this object to do        *
*     their calculations                                                         *
*                                                                                *
*********************************************************************************/

#include "occupancy.h"

using namespace std;


Occupancy::Occupancy() {}


/*    Setters    */

void Occupancy::create()
{
  /* we should have the pointer to all the binding sites, as well as the number
  of nuclei, so this is as simple as looping through binding sites and creating
  an array. Because subgroups and quenchers need to do their calculations using
  pointers to binding sites, so we have two options:
  - have binding sites point to their own occupancy structure
  - have subgroups/quenching deal with each seperately
  For now I think the former possibility is the more elegant */

  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      update(gene, tf);
    }
  }
}

void Occupancy::update(Gene& gene, TF& tf)
{
  site_ptr_vector& sites = bindings->getSites(gene, tf);
  int nsites = sites.size();
  OccData& occ_data_ref = occ_data[&gene][&tf];
  
  occ_data_ref.kv.resize(nsites);
  occ_data_ref.occupancy.resize(nsites);
  occ_data_ref.effective_occupancy.resize(nsites);
  
  for (int i=0; i<nsites; i++)
  {
    occ_data_ref.kv[i].resize(nnuc);
    occ_data_ref.occupancy[i].resize(nnuc);
    occ_data_ref.effective_occupancy[i].resize(nnuc);
    
    sites[i]->kv = &occ_data_ref.kv[i];
    sites[i]->occupancy = &occ_data_ref.occupancy[i];
    sites[i]->effective_occupancy= &occ_data_ref.effective_occupancy[i];
  }
  calcKV(gene, tf);
}

// must call this every time binding sites might change
void Occupancy::update(TF& tf)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    update(gene, tf);
  }
}

void Occupancy::updateKV(TF& tf)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    calcKV(gene, tf);
  }
}
    

void Occupancy::calcKV(Gene& gene, TF& tf)
{
  site_ptr_vector& tmp_sites = bindings->getSites(gene, tf);
  OccData& occ_data_ref      = occ_data[&gene][&tf];
  
  const string& tf_name = tf.getName();
  vector<double>& conc  = tf_data->getConc(tf_name);
  
  int nsites = tmp_sites.size();
  for (int i=0; i<nsites; i++)
  {
    BindingSite* site = tmp_sites[i].get();
    for (int j=0; j<nnuc; j++)
      occ_data_ref.kv[i][j] = conc[j]*site->K_exp_part_times_kmax;
  }
}


/*    Move Functions    */


void Occupancy::saveOccupancy()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      OccData& tmp_occ_data = occ_data[&gene][&tf];
      tmp_occ_data.saved_occupancy = tmp_occ_data.occupancy;
    }
  }
}
        

void Occupancy::restoreOccupancy()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      OccData& tmp_occ_data = occ_data[&gene][&tf];
      tmp_occ_data.occupancy = tmp_occ_data.saved_occupancy;
    }
  }
}


/* this function saves effective occupancy and resets occupancy to 
simple fractional occupancy */
void Occupancy::saveEffectiveOccupancy()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      OccData& tmp_occ_data = occ_data[&gene][&tf];
      tmp_occ_data.saved_effective_occupancy = tmp_occ_data.effective_occupancy;
      tmp_occ_data.effective_occupancy       = tmp_occ_data.occupancy;
    }
  }
}


void Occupancy::restoreEffectiveOccupancy()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      OccData& tmp_occ_data = occ_data[&gene][&tf];
      tmp_occ_data.effective_occupancy = tmp_occ_data.saved_effective_occupancy;
    }
  }
}


void Occupancy::saveAllOccupancy()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      OccData& tmp_occ_data = occ_data[&gene][&tf];
      tmp_occ_data.saved_kv                  = tmp_occ_data.kv;
      tmp_occ_data.saved_effective_occupancy = tmp_occ_data.effective_occupancy;
      tmp_occ_data.saved_occupancy           = tmp_occ_data.occupancy;
      tmp_occ_data.effective_occupancy       = tmp_occ_data.occupancy;
    }
  }
}

void Occupancy::restoreAllOccupancy()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      OccData& tmp_occ_data = occ_data[&gene][&tf];
      tmp_occ_data.kv                  = tmp_occ_data.saved_kv;
      tmp_occ_data.effective_occupancy = tmp_occ_data.saved_effective_occupancy;
      tmp_occ_data.occupancy           = tmp_occ_data.saved_occupancy;
    }
  }
}


void Occupancy::printOccupancy(Gene& gene, ostream& os)
{
  int w = 10;
  os << setprecision(3);
  
  int ntfs = tfs->size();
  // loop through everything and figure out the header
  os << setw(w) << "id";
  
  // find all the ids, using map to prevent duplicates
  // simultaneously print binding site headers
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tmp_sites = bindings->getSites(gene, tf);
    int nsites = tmp_sites.size();
    for (int j=0; j<nsites; j++)
    {
      string header = tf.getName();
      stringstream convert;
      convert << j;
      string num = convert.str();
      header += num;
      //header << i;
      os << setw(w) <<  header;
    }
  }
  
  os << endl;
  
  for (int nuc=0; nuc<nnuc; nuc++)
  {
    //int id = idx_2_id[nuc];
    os << setw(w) << nuc;
   
    for (int i=0; i<ntfs; i++)
    {
      TF& tf = tfs->getTF(i);
      OccData& tmp_occ = occ_data[&gene][&tf];
      int nsites = tmp_occ.occupancy.size();
      for (int j=0; j<nsites; j++)
        os << setw(w) << tmp_occ.occupancy[j][nuc];
    }
    os << endl;
  }
}


void Occupancy::printEffectiveOccupancy(Gene& gene, ostream& os)
{
  int w = 10;
  os << setprecision(3);
  
  int ntfs = tfs->size();
  // loop through everything and figure out the header
  os << setw(w) << "id";
  
  // find all the ids, using map to prevent duplicates
  // simultaneously print binding site headers
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tmp_sites = bindings->getSites(gene, tf);
    int nsites = tmp_sites.size();
    for (int j=0; j<nsites; j++)
    {
      string header = tf.getName();
      stringstream convert;
      convert << j;
      string num = convert.str();
      header += num;
      //header << i;
      os << setw(w) <<  header;
    }
  }
  
  os << endl;
  
  for (int nuc=0; nuc<nnuc; nuc++)
  {
    //int id = idx_2_id[nuc];
    os << setw(w) << nuc;
   
    for (int i=0; i<ntfs; i++)
    {
      TF& tf = tfs->getTF(i);
      OccData& tmp_occ = occ_data[&gene][&tf];
      int nsites = tmp_occ.effective_occupancy.size();
      for (int j=0; j<nsites; j++)
        os << setw(w) << tmp_occ.effective_occupancy[j][nuc];
    }
    os << endl;
  }
}


