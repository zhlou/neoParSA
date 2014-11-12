/*********************************************************************************
*                                                                                *
*     bindingsite.h                                                              *
*                                                                                *
*     Contains the structure definition for a binding site.                      *
*     Pretty dull for now                                                        *
*                                                                                *
*********************************************************************************/

#ifndef BINDINGSITE_H
#define BINDINGSITE_H

#include "TF.h"
#include <cstdlib>
#include <map>


struct OccData; // the information about occupancy over nuclei

struct BindingSite
{
  TF* tf;
  char orientation;
  int  m;
  int  n;
  double score;
  double K_exp_part;  // a number from 0-1 where 0 is worst, 1 is best
  double K_exp_part_times_kmax;
  
  /* we need to have information very accessible from binding sites so it is more
  convenient to store some information here on binding site creation */
  vector<double> kv;
  vector<double> total_occupancy;
  
  vector< vector<double> > mode_occupancy;
  vector< vector<double> > effective_occupancy;
  
  vector<double> saved_kv;
  vector<double> saved_total_occupancy;
  
  vector< vector<double> > saved_mode_occupancy;
  vector< vector<double> > saved_effective_occupancy;
};

void printSiteHeader(ostream& os);

void printSite(BindingSite& b, ostream& os);

void printSite(BindingSite* b, ostream& os);

typedef boost::shared_ptr<BindingSite> site_ptr;
typedef vector<site_ptr> site_ptr_vector;







#endif
