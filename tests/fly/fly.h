/*
 * fly.h
 *
 *  Created on: Feb 1, 2013
 *      Author: zhlou
 */

#ifndef FLY_H_
#define FLY_H_

#include "maternal.h"
#include "zygotic.h"
#include "scoring.h"

#include <string>
#include <vector>
#include <libxml/tree.h>

using namespace std;

struct fly_params
{
    FILE          *infile;                     /* pointer to input data file */

    FILE          *dumpptr;          /* pointer for dumping raw model output */

    FILE          *slog;                          /* solver log file pointer */

  /* the follwoing few variables are read as arguments to command line opts  */

    int           ndigits;     /* precision of output, default = 12 */
    int           gutndigits;      /* precision for gut output, def = 6 */
    //int           penaltyflag;              /* flag for printing penalty */
    // int           rmsflag;     /* flag for printing root mean square */
    int           gutflag;         /* flag for root square diff guts */
    int olddivstyle;
    GFunc GofU;

    double        stepsize;                   /* stepsize for solver */
    double        accuracy;                /* accuracy for solver */

    string section_title;                  /* parameter section name */
    string solver_name; // name of solver used
    string infile_name;
    string outfile_name;
    int debug;


};

// set default parameters and read from xml for any parameters that get override
fly_params readFlyParams(xmlNode *docroot, const char* default_section="eqparms");

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

    struct ParamList{
      double    *param; /* pointers to parameters to be tweaked */
      Range     *param_range; /* pointers to corresponding range limits */
      double previous; // new member: store previous value for restore
      bool canRestore; // new member: store if can be restored
    };
    // ptab is now a vector so Translate is easier to implement. It can
    // automatically be deleted also. In addition, since all items in ptab
    // are actually pointing to other parts of the fly problem, ParamList
    // does not need a destructor.
    vector<ParamList> ptab;

    int nparams; /* number of parameters to be tweaked */


    void Translate(vector<ParamList> &tab);


    string infile;
    string outfile;
    int ndigits;
    maternal TheMaternal;
    zygotic zygote;
    scoring score;
    const TheProblem &defs;
    double chisq;
    bool score_valid;
    double updateChisq() {chisq = score.Score(); score_valid = true; return chisq;}

    Tweak ReadTweak(const fly_params& params);

    Tweak tweak; /* tells the annealer which parameters to tweak */

    void WriteParameters(const char *filename, EqParms *p, const char *title,
                         int ndigits);

public:
    fly(const fly_params &params);
    ~fly();
    int getDimension(){return nparams;}; // returns the dimension of the problem
    double get_score(); // returns the energy (score) of the problem
    double get_rms(); // returns the RMS
    void generateMove(int index, double delta);
    // make the perturbation on parameter index with theta_bar
    void restoreMove(int index);
    int getStateSize(){return nparams * sizeof(double);}; // returns the byte count of the state
    void serialize(void *buf) const; // serialize its state to buf
    void deserialize(void const *buf); // inflate buf to a new state. calculation
                                 // of new score is not required here.
    void state2theta(void const *buf, double *theta);
    void writeAnswer(const char* title);
    void scramble();
};

#endif /* FLY_H_ */
