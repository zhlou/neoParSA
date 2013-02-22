/*
 * flyData.cpp
 *
 *  Created on: Feb 14, 2013
 *      Author: zhlou
 */

#include "flyData.h"

#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <cmath>
#include <cstring>

using namespace std;
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

FILE* FindSection(FILE* fp, const char* input_section)
{
    int c; /* input happens character by character */
    int nsought; /* holds length of section title */
    char *base; /* string for section title */
    long looksite; /* file position indicator */

    rewind(fp); /* start looking at the beginning of the file */

    nsought = strlen(input_section);
    base = (char *) calloc(MAX_RECORD, sizeof(char));

    /*** while loop goes through the file character by character... ************/
    while ((c = getc(fp)) != EOF) { /* ...until if finds a '$' */

        if (c == '$') { /* found a sectioning control string */
            looksite = ftell(fp); /* where are we? */
            base = fgets(base, MAX_RECORD, fp); /* get sect title (without $)*/
            if (!(strncmp(base, input_section, nsought))) {
                fseek(fp, looksite, 0); /* found the sought string: reposi- */
                fscanf(fp, "%*s\n"); /* tion the pointer to '$', then */
                free(base);
                return (fp); /* advance the pointer to after the */
            } /* section title and return it */
            else { /* didn't find it: skip this control */
                fseek(fp, looksite, 0); /* record, keep looking */
                fscanf(fp, "%*s"); /* NOTE: "%*s" advances pointer */
            } /* without assignment */
        }
    }

    free(base);

    return (NULL); /* couldn't find the right section */
}

/*** ERROR HANDLING FUNCTIONS **********************************************/

/*** The following two routines print error messages with value. 'error' ***
 *   then exits, while 'warning' returns to the calling function.          *
 *                                                                         *
 *   Both functions are called analogous to the way you would call printf. *
 *   The only conversion specs that are legal are %c, %d, %f, %g and %s.   *
 *   No modifiers (h,l or float precision or size specifiers are allowed). *
 *                                                                         *
 *   Both functions can only handle %c, %d, %f, %g, %s (char, string, int  *
 *   (includes long integers too) & double); it does not yet work for      *
 *   floats or long doubles or the more esoteric format specifiers of      *
 *   printf.                                                               *
 *                                                                         *
 *   No % sign means just print msg and exit.                              *
 *                                                                         *
 ***************************************************************************/

void error(const char *format, ...)
{
    int i; /* loop counter */
    char *msg; /* same as format but with appended /n */
    char *msg_ptr; /* used to parse format string */
    char **msg_array; /* holds different parts of the parsed msg */
    char *msg_types; /* holds types of the messages to be printed */
    size_t length; /* length of the format string */
    int arg_count = 0; /* number of additional arguments */
    va_list arg_ptr; /* pointer to additional arguments */

    /* initialize variable argument pointer */

    va_start(arg_ptr, format); /* initialize the argument pointer */

    /* allocate memory for message arrays */

    msg_array = (char **) calloc(MAX_ARGS, sizeof(char *));
    for (i = 0; i < MAX_ARGS; i++)
        msg_array[i] = (char *) calloc(MAX_RECORD, sizeof(char));
    msg_types = (char *) calloc(MAX_ARGS, sizeof(char));

    /* copy the format string to msg and attach a newline; note that we allo-  *
     * cate length+5 bytes for the string; I dunno exactly why that is nece-   *
     * ssary but 5 is the minimum number for which Third Degree will not com-  *
     * plain                                                                   */

    length = strlen(format);
    msg = (char *) calloc((length + 5), sizeof(char));
    msg = strcpy(msg, format);
    msg = strcat(msg, "\n");

    /* finds '%' in msg by reverse search; then copies subportions of msg to   *
     * the msg_array and the corresponding types to msg_types; these arrays    *
     * contain the parts of the error message *in reverse*, from end to be-    *
     * ginning; a new /0 is then inserted at the position of the % and the     *
     * search is repeated until no more % can be found in the format string    *
     * (the upper limit of %'s is given by MAX_ARGS in error.h)                */

    while ((msg_ptr = strrchr(msg, '%')) != NULL) {
        if (arg_count >= MAX_ARGS) {
            fprintf(stderr, "error: too many arguments (max. %d)!\n", MAX_ARGS);
            exit(1);
        }
        msg_array[arg_count] = strcpy(msg_array[arg_count], msg_ptr);
        msg_types[arg_count] = *(msg_ptr + 1);
        ++arg_count;
        *msg_ptr = '\0';
    }

    /* then print error message according to format string; note that only     *
     * self-promoting types are allowed for va_arg(); that's why we use int    *
     * for char, since char would get promoted to an int in the process and    *
     * some compilers (like gcc) have a problem with that                      */

    fprintf(stderr, msg);
    if (arg_count == 0) {
        exit(1);
    } else
        for (i = 1; i <= arg_count; i++)
            switch (msg_types[arg_count - i]) {
                case 's':
                    fprintf(stderr, msg_array[arg_count-i], va_arg(arg_ptr, char *));
                    break;
                case 'c':
                case 'd':
                    fprintf(stderr, msg_array[arg_count-i], va_arg(arg_ptr, int));
                    break;
                case 'f':
                case 'g':
                    fprintf(stderr, msg_array[arg_count-i], va_arg(arg_ptr, double));
                    break;
                default:
                    fprintf(stderr,
                            "error: error in message parsing, bailing out!\n");
                    exit(1);
            }

    /* clean up and go home */

    va_end(arg_ptr);

    for (i = 0; i < MAX_ARGS; i++)
        free(msg_array[i]);
    free(msg_array);
    free(msg_types);
    free(msg);

    exit(1);

}

void warning(const char *format, ...)
{
    int i; /* loop counter */
    char *msg; /* same as format but with appended /n */
    char *msg_ptr; /* used to parse format string */
    char **msg_array; /* holds different parts of the parsed msg */
    char *msg_types; /* holds types of the messages to be printed */
    size_t length; /* length of the format string */
    int arg_count = 0; /* number of additional arguments */
    va_list arg_ptr; /* pointer to additional arguments */

    /* initialize variable argument pointer */

    va_start(arg_ptr, format); /* initialize the argument pointer */

    /* allocate memory for message arrays */

    msg_array = (char **) calloc(MAX_ARGS, sizeof(char *));
    for (i = 0; i < MAX_ARGS; i++)
        msg_array[i] = (char *) calloc(MAX_RECORD, sizeof(char));
    msg_types = (char *) calloc(MAX_ARGS, sizeof(char));

    /* copy the format string to msg and attach a newline; note that we allo-  *
     * cate length+5 bytes for the string; I dunno exactly why that is nece-   *
     * ssary but 5 is the minimum number for which Third Degree will not com-  *
     * plain                                                                   */

    length = strlen(format);
    msg = (char *) calloc((length + 5), sizeof(char));
    msg = strcpy(msg, format);
    msg = strcat(msg, "\n");

    /* finds '%' in msg by reverse search; then copies subportions of msg to   *
     * the msg_array and the corresponding types to msg_types; these arrays    *
     * contain the parts of the error message *in reverse*, from end to be-    *
     * ginning; a new /0 is then inserted at the position of the % and the     *
     * search is repeated until no more % can be found in the format string    *
     * (the upper limit of %'s is given by MAX_ARGS in error.h)                */

    while ((msg_ptr = strrchr(msg, '%')) != NULL) {
        if (arg_count >= MAX_ARGS) {
            fprintf(stderr, "warning: too many arguments (max. %d)!\n",
                    MAX_ARGS);
            exit(1);
        }
        msg_array[arg_count] = strcpy(msg_array[arg_count], msg_ptr);
        msg_types[arg_count] = *(msg_ptr + 1);
        ++arg_count;
        *msg_ptr = '\0';
    }

    /* then print error message according to format string; note that only     *
     * self-promoting types are allowed for va_arg(); that's why we use int    *
     * for char, since char would get promoted to an int in the process and    *
     * some compilers (like gcc) have a problem with that                      */

    fprintf(stderr, msg);
    if (arg_count == 0) {
        return;
    } else
        for (i = 1; i <= arg_count; i++)
            switch (msg_types[arg_count - i]) {
                case 's':
                    fprintf(stderr, msg_array[arg_count-i], va_arg(arg_ptr, char *));
                    break;
                case 'c':
                case 'd':
                    fprintf(stderr, msg_array[arg_count-i], va_arg(arg_ptr, int));
                    break;
                case 'f':
                case 'g':
                    fprintf(stderr, msg_array[arg_count-i], va_arg(arg_ptr, double));
                    break;
                default:
                    fprintf(stderr,
                            "warning: error in message parsing ... bailing out!\n");
                    exit(1);
            }

    /* clean up and go home */

    va_end(arg_ptr);

    for (i = 0; i < MAX_ARGS; i++)
        free(msg_array[i]);
    free(msg_array);
    free(msg_types);
    free(msg);

    return;

}

