/*
 * fly.h
 *
 *  Created on: Feb 1, 2013
 *      Author: zhlou
 */

#ifndef FLY_H_
#define FLY_H_

#include <string>

using namespace std;

class maternal;
class zygotic;
class scoring;

struct fly_params
{
    fly_params(); // initialize parameters to defaults
    FILE          *infile;                     /* pointer to input data file */

    FILE          *dumpptr;          /* pointer for dumping raw model output */

    FILE          *slog;                          /* solver log file pointer */

  /* the follwoing few variables are read as arguments to command line opts  */

    int           ndigits;     /* precision of output, default = 12 */
    int           gutndigits;      /* precision for gut output, def = 6 */
    int           penaltyflag;              /* flag for printing penalty */
    int           rmsflag;     /* flag for printing root mean square */
    int           gutflag;         /* flag for root square diff guts */

    double        stepsize;                   /* stepsize for solver */
    double        accuracy;                /* accuracy for solver */

    string section_title;                  /* parameter section name */
    string solver_name; // name of solver used
    string infile_name;
    int debug;


};

class fly // This class implements the interfaces and has old translate stuff
{
private:
    /* Tweak struct is for tweaking individual parameters or not; each pointer *
     * below points to an array of ints which represent each paramter          *
     * 1 means tweak, 0 means leave it alone                                   *
     * see ../doc/dataformatX.X for further details on the tweak section       */

    struct Tweak {
        int *Rtweak; /* which Rs to be tweaked */
        int *Ttweak; /* which Ts to be tweaked */
        int *Etweak; /* which Es to be tweaked */
        int *mtweak; /* which ms to be tweaked */
        int *htweak; /* which hs to be tweaked */
        int *dtweak; /* which ds to be tweaked */
        int *lambdatweak; /* which lambdas to be tweaked */
        int *tautweak; /* which taus to be tweaked */
    };

    maternal *TheMaternal;
    zygotic *zygote;
    scoring *score;
    TheProblem &defs;

    Tweak tweak; /* tells the annealer which parameters to tweak */
public:
    fly(const fly_params &params);

};

#endif /* FLY_H_ */
