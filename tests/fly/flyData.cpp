/*
 * flyData.cpp
 *
 *  Created on: Feb 14, 2013
 *      Author: zhlou
 */

#include "flyData.h"

/*** FreeSolution: frees memory of the solution structure created by *******
 *                 Blastoderm() or gut functions                           *
 ***************************************************************************/
void FreeSolution(NArrPtr *solution)
{
    int i;                                                   /* loop counter */

    for(i=0; i < solution->size; i++)
        free(solution->array[i].state.array);
    free(solution->array);
}

void FreeFacts(DataTable *D)
{

    int i,j;

    for (i=0; i < D->size; i++)
                free(D->record[i].array);

    free(D->record);
    free(D);

    return;

}

int descend(const int *x, const int *y)
{

    if (*x < *y) return 1;
    else if (*x > *y) return -1;
    else return 0;

}

/* Converts a data table into NarrPtr and also returns the index when
the number of nuclei is maximum */

NArrPtr Dat2NArrPtr(DataTable *table, int *maxind)
{

    int i,j,max;
    NArrPtr input;

    max = 0;
    *maxind = 0;
    input.size = table->size;
    input.array = (NucState *) calloc(input.size, sizeof(NucState));

    for (i=0; i < table->size; i++)
    {

/*      printf("Heloo: %d %f\n",table->record[i].size,
                                    table->record[i].time);*/

        input.array[i].time=table->record[i].time+EPSILON;
        input.array[i].state.size=table->record[i].size;

        if (max < input.array[i].state.size) {
            max = input.array[i].state.size;
            *maxind = i;
        }

        input.array[i].state.array=(double *)
                    calloc(input.array[i].state.size, sizeof(double));
        for (j=0; j < table->record[i].size; j++)
             input.array[i].state.array[j]=table->record[i].array[j].conc;

    }

    return input;
}

/*** ConvertAnswer: little function that gets rid of bias times, division **
 *                  times and such and only returns the times in the tab-  *
 *                  times struct as its output; this is used to produce    *
 *                  unfold output that contains only the requested times;  *
 *                  it also makes sure that we return the right solution   *
 *                  at cell division, i.e. the solution right *after* the  *
 *                  cells have divided                                     *
 ***************************************************************************/

NArrPtr ConvertAnswer(NArrPtr answer, DArrPtr tabtimes)
{
    int       i;                                  /* loop counter for answer */
    int       outindex = 0;                     /* loop counter for tabtimes */
    double    nexttime;                 /* used below to avoid memory errors */
    NArrPtr   outtab;                               /* answer to be returned */

    /* allocate the outtab array */

    outtab.size = tabtimes.size;
    outtab.array = (NucState *)calloc(outtab.size, sizeof(NucState));

    /* the following loop goes through all answer elements and copies only     */
    /* those which correspond to a tabtime into the outtab struct              */

    for ( i=0; (i<answer.size) && (outindex<tabtimes.size); i++ ) {

        /* this kludge makes sure that we don't bomb below when there is no next   *
         * time anymore                                                            */

        if ( i == (answer.size-1) )
            nexttime = 999999999;
        else                                                   /* if statement */
            nexttime = answer.array[i+1].time;

        /* this if makes sure that we only return the requested times and that     *
         * we return the solution *after* the cells have divided for times that    *
         * correspond exactly to the cell division time (i.e. the end of mitosis); *
         * note that times before and after cell divisions are separated by EPIS-  *
         * LON and that BIG_EPSILON is always (much) bigger than EPSILON (but      *
         * still small enough for all practical purposes...)                       */

        if ( (fabs(answer.array[i].time - tabtimes.array[outindex])
                < BIG_EPSILON)
                && (fabs(answer.array[i].time - nexttime)
                        > BIG_EPSILON) ) {
            outtab.array[outindex].time = answer.array[i].time;
            outtab.array[outindex].state.size = answer.array[i].state.size;
            outtab.array[outindex].state.array = answer.array[i].state.array;
            outindex++;
        }
    }
    return outtab;
}
