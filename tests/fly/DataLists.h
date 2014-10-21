/*
 * DataLists.h
 *
 *  Created on: Feb 5, 2013
 *      Author: zhlou
 */

#ifndef DATALISTS_H_
#define DATALISTS_H_

#include <cstdio>
using namespace std;
/* linked lists for reading bicoid (Blist) and bias/facts (Dlist) */

struct Blist
{ /* linked list used to read bicoid */
    unsigned int lineage;
    double conc;
    struct Blist *next;
};

struct Dlist
{ /* linked list used to read bias & facts */
    unsigned int lineage;
    double *d;
    struct Dlist *next;
};

/* linked list used to read in genotype specific sections from datafile */

struct Slist
{
    char *genotype; /* genotype string (see dataformatX.X) */
    char *bcd_section; /* bcd section title */
    char *bias_section; /* bias section title */
    char *fact_section; /* fact section title */
    char *hist_section; /* fact section title */
    char *ext_section; /* fact section title */
    struct Slist *next;
};

/* Following functions are utility functions for different linked lists ****
 * which are used to read in data of unknown size. All utility functions   *
 * follow the same scheme (X stands for the different letters below):      *
 *                                                                         *
 * - init_Xlist:         allocates first element and returns pointer to it *
 * - addto_Xlist:        adds the adduct to an already existing linkd list *
 * - free_Xlist:         frees memory of linked list                       *
 * - count_Xlist:        counts elements in the linked list                *
 *                                                                         *
 ***************************************************************************/

/* Utility functions for Blist (for reading bicoid) */

Blist        *init_Blist    (void);
Blist        *addto_Blist   (Blist *start, Blist *adduct);
void         free_Blist     (Blist *start);
int          count_Blist    (Blist *start);


/* Utility functions for Dlist (for reading bias and facts) */

Dlist        *init_Dlist    (int size);
Dlist        *addto_Dlist   (Dlist *start, Dlist *adduct);
void         free_Dlist     (Dlist *start);
int          count_Dlist    (Dlist *start);


/* Utility functions for Slist (for reading genotypes ) */

Slist        *init_Slist    (void);
Slist        *addto_Slist   (Slist *start, Slist *adduct);
void         free_Slist     (Slist *start);
int          count_Slist    (Slist *start);

const int IGNORE = -1;
Dlist *ReadData(FILE *fp , char *section, int *ndp, int ngenes);
Dlist *ReadInterpData(FILE *fp , char *section, int num_genes, int *ndp);

#endif
