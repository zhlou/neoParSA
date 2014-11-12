/*********************************************************************************
*                                                                                *
*     bindingsite.cpp                                                            *
*                                                                                *
*     Contains the structure definition for a binding site.                      *
*     Pretty dull for now                                                        *
*                                                                                *
*********************************************************************************/

#include "bindingsite.h"
#include <iostream>
#include <iomanip>
#include "TF.h"

void printSiteHeader(ostream& os)
{
  int w = 10;
  os << setprecision(3)
     << setw(w) << "name"
     << setw(w) << "start"
     << setw(w) << "end"
     << setw(w) << "score"
     << setw(w) << "K"
     << setw(w) << "K*kmax";
}

void printSite(BindingSite& b, ostream& os)
{
  int w = 10;
  os << setprecision(3)
     << setw(w) << b.tf->getName()
     << setw(w) << b.m 
     << setw(w) << b.n 
     << setw(w) << b.score 
     << setw(w) << b.K_exp_part
     << setw(w) << b.K_exp_part_times_kmax;
}

void printSite(BindingSite* b, ostream& os) { printSite(*b, os); }
