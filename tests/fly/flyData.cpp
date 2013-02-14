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
