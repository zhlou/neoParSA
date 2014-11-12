/*********************************************************************************
*                                                                                *
*     getsequence.cpp                                                            *
*                                                                                *
*     gets a sequence from a 2bit file using coordinates                         *
*                                                                                *
*********************************************************************************/

# include "twobit.h"

int main(int argc, char* argv[])
{

  string fname = argv[1];
  TwoBit genome(fname);
  
  cout << endl;
  
  genome.printHeader(cout);
  genome.printIndex(cout);
 
  cout << endl;
  
  return 0;
}
