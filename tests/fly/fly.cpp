/*
 * fly.cpp
 *
 *  Created on: Feb 18, 2013
 *      Author: zhlou
 */

#include "fly.h"
#include "maternal.h"
#include "zygotic.h"
#include "scoring.h"

fly_params::fly_params() : // initialize default parameters
        section_title("eqparms"), solver_name("Rk4")
{
    ndigits = 12;
    gutndigits = 6;
    penaltyflag = 0;
    rmsflag = 1;
    gutflag = 0;
    stepsize = 1.;
    accuracy = 0.001;
    infile = NULL;
    dumpptr = NULL;
    slog = NULL;
    debug = 0;
}

fly::fly(const fly_params &params) {
    TheMaternal = new maternal(params.infile);
    zygote = new zygotic(*TheMaternal, params.infile,
                         params.section_title.c_str(), params.debug,
                         params.solver_name.c_str());
    score = new scoring(params.infile, *zygote, params.stepsize,
                        params.accuracy,params.slog, params.infile_name.c_str(),
                        params.debug);

}
