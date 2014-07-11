/*
 * DataLists.cpp
 *
 *  Created on: Feb 5, 2013
 *      Author: zhlou
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include "DataLists.h"
#include "flyData.h"
using namespace std;

/* FOLLOWING FUNCTIONS ARE UTILITY FUNCTIONS FOR DIFFERENT LINKED LISTS ****
 * which are used to read in data of unknown size. All utility functions   *
 * follow the same scheme (X stands for the different letters below):      *
 *                                                                         *
 * - init_Xlist:         allocates first element and returns pointer to it *
 * - addto_Xlist:        adds the adduct to an already existing linkd list *
 * - free_Xlist:         frees memory of linked list                       *
 * - count_Xlist:        counts elements in the linked list                *
 *                                                                         *
 ***************************************************************************/

/*** Utility functions for Blist *******************************************/

Blist *init_Blist(void)
{
  Blist     *p;

  if( (p =(Blist *)malloc(sizeof(Blist))) ) {

    p->lineage = 0;
    p->conc    = 0;
    p->next    = NULL;

  }
  else
    error("init_Blist: couldn't allocate!");

  return p;
}



Blist *addto_Blist(Blist *start, Blist *adduct)
{
  Blist     *current;

  if( !start )
    return adduct;

  current = start;

  while (current->next) {
    current = current->next;
  }

  current->next = adduct;
  return start;
}



void free_Blist(Blist *start)
{
  if(start->next)
    free_Blist(start->next);

  free(start);
}



int count_Blist(Blist *start)
{
  int        n = 0;

  while( start != NULL ) {
    n++;
    start = start->next;
  }

  return n;
}



/*** Utility functions for Dlist *******************************************/

Dlist *init_Dlist(int size)
{
  Dlist    *p;

  if( (p=(Dlist *)malloc(sizeof(Dlist))) ) {

    if ( !(p->d=(double *)calloc(size, sizeof(double))) )
      error("initDlist: couldn't allocate!");

    p->lineage = 0;
    p->next    = NULL;

}
  else
    error("initDlist: couldn't allocate!");

  return p;
}



Dlist *addto_Dlist(Dlist *start, Dlist *adduct)
{
  Dlist         *current;

  if(!start)
    return adduct;

  current = start;

  while (current->next) {
    current = current->next;
  }

  current->next = adduct;
  return start;
}


void free_Dlist(Dlist *start)
{
  if(start->next)
    free_Dlist(start->next);

  free(start->d);
  free(start);
}



int count_Dlist(Dlist *start)
{
  int n = 0;
  while(start != NULL) {
    n++;
    start = start->next;
  }

  return n;
}



/*** Utility functions for Slist *******************************************/

Slist *init_Slist(void)
{
  Slist     *p;

  if( (p=(Slist *)malloc(sizeof(Slist))) ) {

    p->bias_section =       NULL;
    p->fact_section =       NULL;
    p->bcd_section  =       NULL;
    p->hist_section  =       NULL;
    p->ext_section  =       NULL;
    p->genotype     =       NULL;
    p->next         =       NULL;

  }
  else
    error("initSlist: couldn't allocate");

  return p;
}



Slist *addto_Slist(Slist *start, Slist *adduct)
{
  Slist      *current;

  if(!start)
    return adduct;

  current = start;

  while (current->next) {
    current = current->next;
  }

  current->next = adduct;
  return start;
}



void free_Slist(Slist *start)
{
  if(start->next)
    free_Slist(start->next);

  free(start->bias_section);
  free(start->fact_section);
  free(start->bcd_section);
  free(start->hist_section);
  free(start->ext_section);
  free(start->genotype);
  free(start);
}



int count_Slist(Slist *start)
{
  int n = 0;
  while(start != NULL) {
    n++;
    start = start->next;
  }

  return n;
}

/*** ReadData: reads in a data or bias section and puts it in a linked *****
 *             list of arrays, one line per array; ndp is used to count    *
 *             the number of data points in a data file (ndp), which is    *
 *             used to calculate the root mean square (RMS) if required    *
 *                                                                         *
 *             ReadData allows comments that start with a letter or punc-  *
 *             tuation mark, and treats lines starting with a number or    *
 *             .num or -.num as data. It expects to read an int and        *
 *             ngenes + 1 data points: lineage, time and ngenes protein    *
 *             concentrations. It expects data to be in increasing spatial *
 *             and temporal order.                                         *
 ***************************************************************************/

Dlist *ReadData(FILE *fp , char *section, int *ndp, int ngenes)
{
  int     c;                           /* holds char for parser            */
  int     lead_punct;                  /* flag for leading punctuation     */
  int     i;                           /* loop counter                     */

  char    *base;                       /* pointer to beginning of string   */
  char    *record;                     /* pointer to string, used as       */
                                       /* counter                          */
  Dlist   *current;                    /* holds current element of Dlist   */
  Dlist   *inlist;                     /* holds whole read Dlist           */

/* the following chars are used for parsing lines of data (record)         */
/* using sscanf                                                            */

  char       *fmt  = NULL;             /* used as format for sscanf        */
  char       *skip = NULL;             /* skip this stuff!                 */

  const char init_fmt[] = "%*d ";      /* initial format string            */
  const char skip_fmt[] = "%*lg ";     /* skip one more float              */
  const char read_fmt[] = "%lg ";      /* skip one more float              */

  if ( (fp=FindSection(fp,section)) ) {                 /* position the fp */

    base    = (char *)calloc(MAX_RECORD, sizeof(char));
    fmt     = (char *)calloc(MAX_RECORD, sizeof(char));
    skip    = (char *)calloc(MAX_RECORD, sizeof(char));

    current = NULL;
    inlist  = NULL;

/* while loop: reads and processes lines from file until sections ends *****/

    while ( strncmp(( base=fgets(base, MAX_RECORD, fp)), "$$", 2)) {

      record     = base;                  /* base always points to start   */
      lead_punct = 0;                     /* of string                     */

/* while loop: parses and processes each line from the data file *************/

      c=(int)*record;

      while ( c != '\0' ) {

    if( isdigit(c) ) {                            /* number means data */

      record = base;                  /* reset pointer to start of str */
      current = init_Dlist(ngenes+1);

      if ( 1 != sscanf(record, "%d ", &(current->lineage)) )
        error("ReadData: error reading %s", base);

/* the following loop reads a line of data of variable length into a d-    */
/* array, skipping lineage number and all previously read data values      */

      skip = strcpy(skip, init_fmt); /* format to be skipped by sscanf */

      for(i=0; i < (ngenes + 1); i++) {

        fmt  = strcpy(fmt, skip);     /* all this stuff here is to get */
        fmt  = strcat(fmt, read_fmt);  /* the fmt string for sscanf to */
        skip = strcat(skip, skip_fmt);           /* the correct format */

        if ( 1 != sscanf(record, (const char *)fmt, &(current->d[i])) )
          error("ReadData: error reading %s", base);

                                           /* update number of data points */
        if ( ( i != 0) && (current->d[i] != IGNORE) )
          (*ndp)++;

      }
                                      /* now add this to the lnkd list */
      inlist = addto_Dlist(inlist, current);
      break;
    }

    else if ( isalpha(c) ){                    /* letter means comment */
      break;
    }
    else if ( c == '-' ) {                /* next two elsifs for punct */
      if ( ((int)*(record+1)) == '.')
        record++;
      lead_punct = 1;
      c=(int)*(++record);
    }
    else if ( c == '.' ) {
      lead_punct = 1;
      c=(int)*(++record);
    }
    else if ( ispunct(c) ){               /* other punct means comment */
      break;
    }
    else if ( isspace(c) ) {             /* ignore leading white space */
      if ( lead_punct )               /* white space after punct means */
        break;                                              /* comment */
      else {
        c=(int)*(++record);            /* get next character in record */
      }
    }
    else {
      error ("ReadData: Illegal character in %s", base);
    }
      }
    }

    free(base);
    free(fmt);
    free(skip);

    return inlist;

  } else {

    return NULL;

  }
}

Dlist *ReadInterpData(FILE *fp , char *section, int num_genes, int *ndp)
{
  int     c;                           /* holds char for parser            */
  int     lead_punct;                  /* flag for leading punctuation     */
  int     i;                           /* loop counter                     */

  char    *base;                       /* pointer to beginning of string   */
  char    *record;                     /* pointer to string, used as       */
                                       /* counter                          */
  Dlist   *current;                    /* holds current element of Dlist   */
  Dlist   *inlist;                     /* holds whole read Dlist           */

/* the following chars are used for parsing lines of data (record)         */
/* using sscanf                                                            */

  char       *fmt  = NULL;             /* used as format for sscanf        */
  char       *skip = NULL;             /* skip this stuff!                 */

  const char init_fmt[] = "%*d ";      /* initial format string            */
  const char skip_fmt[] = "%*lg ";     /* skip one more float              */
  const char read_fmt[] = "%lg ";      /* skip one more float              */

  if ( (fp=FindSection(fp,section)) ) {                 /* position the fp */

    base    = (char *)calloc(MAX_RECORD, sizeof(char));
    fmt     = (char *)calloc(MAX_RECORD, sizeof(char));
    skip    = (char *)calloc(MAX_RECORD, sizeof(char));

    current = NULL;
    inlist  = NULL;

/* while loop: reads and processes lines from file until sections ends *****/

    while ( strncmp(( base=fgets(base, MAX_RECORD, fp)), "$$", 2)) {

      record     = base;                  /* base always points to start   */
      lead_punct = 0;                     /* of string                     */

/* while loop: parses and processes each line from the data file *************/

      c=(int)*record;

      while ( c != '\0' ) {

    if( isdigit(c) ) {                            /* number means data */

      record = base;                  /* reset pointer to start of str */
      current = init_Dlist(num_genes+1);

      if ( 1 != sscanf(record, "%d ", &(current->lineage)) )
        error("ReadInterpData: error reading %s", base);

/* the following loop reads a line of data of variable length into a d-    */
/* array, skipping lineage number and all previously read data values      */

      skip = strcpy(skip, init_fmt); /* format to be skipped by sscanf */

      for(i=0; i < (num_genes + 1); i++) {

        fmt  = strcpy(fmt, skip);     /* all this stuff here is to get */
        fmt  = strcat(fmt, read_fmt);  /* the fmt string for sscanf to */
        skip = strcat(skip, skip_fmt);           /* the correct format */

        if ( 1 != sscanf(record, (const char *)fmt, &(current->d[i])) )
          error("ReadInterpData: error reading %s", base);

                                           /* update number of data points */
        if ( ( i != 0) && (current->d[i] != IGNORE) )
          (*ndp)++;

      }
                                      /* now add this to the lnkd list */
      inlist = addto_Dlist(inlist, current);
      break;
    }

    else if ( isalpha(c) ){                    /* letter means comment */
      break;
    }
    else if ( c == '-' ) {                /* next two elsifs for punct */
      if ( ((int)*(record+1)) == '.')
        record++;
      lead_punct = 1;
      c=(int)*(++record);
    }
    else if ( c == '.' ) {
      lead_punct = 1;
      c=(int)*(++record);
    }
    else if ( ispunct(c) ){               /* other punct means comment */
      break;
    }
    else if ( isspace(c) ) {             /* ignore leading white space */
      if ( lead_punct )               /* white space after punct means */
        break;                                              /* comment */
      else {
        c=(int)*(++record);            /* get next character in record */
      }
    }
    else {
      error ("ReadInterpData: Illegal character in %s", base);
    }
      }
    }

    free(base);
    free(fmt);
    free(skip);

    return inlist;

  } else {

    return NULL;

  }
}
