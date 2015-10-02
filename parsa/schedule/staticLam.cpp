
#include <stdexcept>
#include <cmath>
#include "xmlUtils.h"
#include "utils/vectorUtils.h"

#include "staticLam.h"

const char * staticLam::name = "staticLam";

staticLam::Param::Param(xmlNode *root, debugStatus in_st, const char *name) :
        st(in_st), logname(name)
{
    xmlNode *section = getSectionByName(root, "staticLam");
    segLength = 100;
    try {
        segLength = getPropInt(section, "segLength");
    } catch (const std::exception &e) {
        // ignored
    }
    adjustAlpha=0;
    try {
        adjustAlpha = getPropInt(section, "adjustAlpha");
    } catch (const std::exception &e) {
        // ignored
    }
    lambda = getPropDouble(section, "lambda");
    minRate = 1e-15;
    try {
        minRate = getPropDouble(section, "minRate");
    } catch (const std::exception &e) {
        // ignored
    }
    filename = (char *)xmlGetProp(section, (xmlChar *)"filename");
    if (NULL == filename) {
        throw std::runtime_error("Error: fail to find filename in staticLam");
    }
}

staticLam::staticLam(Param &param) : 
        debugOut(param.st, param.logname), segLength(param.segLength),
        lambda(param.lambda), minRate(param.minRate), count(0), i(0), 
        alpha(0.23), success(0), adjustAlpha(param.adjustAlpha)
{
    readDoubleVectorFromText(betaVec, variance, param.filename);
    size = betaVec.size();
    b0=betaVec[0];
    bEnd=betaVec.back();
    cEnd=bEnd*bEnd*variance.back();
}

double staticLam::getVar(double beta)
{
    if (beta < b0)
        return v0;
    if (beta >= bEnd) // for Rastrigin, heat capacity b^2*var = Constant for large b
        return cEnd/(beta*beta);
    while(beta>betaVec[i+1]) // find the appropiate section to interpolate
        ++i;
    // now do interpolation
    double b1 = std::log(betaVec[i]);
    double b2 = std::log(betaVec[i+1]);
    double t = (std::log(beta)-b1)/(b2-b1);
    return (1-t)*variance[i]+t*variance[i+1];
}

double staticLam::updateS(const aState &state)
{
    double b = state.s;
    double var = getVar(b);
    double delta = lambda * alpha / (var * std::sqrt(var) * b * b);
    if (delta < minRate * b)
        return b + minRate * b;
    else
        return b + delta;
}

void staticLam::updateStep(bool accept, const aState &)
{
    if (accept)
        ++ success;
}
void staticLam::updateStats(const aState &state)
{
    ++ count;
    if (segLength == count) {
        count = 0;
        calcStats(segLength, state);
    }
}

void staticLam::calcStats(unsigned nsteps, const aState &state)
{
    double ar = (double)success / (double)nsteps;
    success = 0;
    if (adjustAlpha) {
        double d = (1.0 - ar) / (2.0 - ar);
        alpha = 4.0 * ar * d * d;
    }
    debugOut << state.step_cnt << " " << state.s << " " 
             << state.energy << " " << getVar(state.s) 
             << " " << ar << std::endl;

}

void staticLam::initStats(double, double, double initAccRatio, const aState& state)
{
    if (adjustAlpha) {
        double d = (1.0 - initAccRatio) / (2.0 - initAccRatio);
        alpha = 4.0 * initAccRatio * d * d;
    }
}

void staticLam::initStats(const aState &state)
{
    calcStats(state.step_cnt, state);
}
