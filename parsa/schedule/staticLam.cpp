
#include <stdexcept>
#include <cmath>
#include "xmlUtils.h"
#include "utils/vectorUtils.h"

#include "staticLam.h"

const char * staticLam::name = "staticLam"

staticLam::Param:Param(xmlNode *root, debugStatus in_st, const char *name) :
        st(in_st), logname(name)
{
    xmlNode *section = getSectionByName(root, "staticLam");
    segLength = 1;
    try {
        segLength = getPropInt(xmlsection, "segLength");
    } catch (const set::exception &e) {
        // ignored
    }
    lambda = getPropDouble(section, "lambda");
    filename = (char *)xmlGetProp(section, (xmlChar *)"filename");
    if (NULL == filename) {
        throw runtime_error(string("Error: fail to find filename in staticLam"));
    }
}

staticLam::staticLam(Param &param) : 
        debugOut(param.st, param.logname), segLength(param.segLength),
        lambda(param.lambda), count(0), i(0), alpha(0.23), success(0)
{
    readDoubleVectorFromText(betaVec, variance, filename);
    size = betaVec.size();
    b0=betaVec[0];
    bEnd=betaVec.back();
    cEnd=bEnd*bEnd*betaVec.back();
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
    return b + lambda * alpha / (var * std::sqrt(var) * b * b);
}

void staticLam::updateStep(bool accept, const aState &)
{
    if (accept)
        ++ success;
}
