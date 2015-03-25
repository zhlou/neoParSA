/*
 * feedbackMove.hpp
 *
 *  Created on: Jan 21, 2013
 *      Author: zhlou
 */
#include <limits>
#include <cstdlib>
#include <cmath>
using namespace std;

//template<class Problem>
//const double feedbackMove<Problem>::theta_min = 0.;

template<class Problem>
const char *feedbackMove<Problem>::name = "feedbackMove";

template<class Problem>
feedbackMove<Problem>::feedbackMove(Problem& in_problem,
                                    unirandom& in_rand,
                                    xmlNode* root) :
        problem(in_problem), rnd(in_rand)
{
    nparams = problem.getDimension();
    index = -1;
    success = new long[nparams];
    moves = new long[nparams];
    theta_bars = new double[nparams];
    theta_mins = new double[nparams];
    theta_maxs = new double[nparams];
    for (int i = 0; i < nparams; ++i) {
        success[i] = 0;
        moves[i] = 0;
        theta_bars[i] = 1.0;
        theta_mins[i] = 0.;
        theta_maxs[i] = numeric_limits<double>::max();
    }
    energy = problem.get_score();
    prev_energy = energy;
    sweep = 0;

    xmlNode *section = NULL;
    if (root != NULL) {
        section = root->children;
        while (section != NULL) {
            if (!xmlStrcmp(section->name, (xmlChar *)"move"))
                break;
            section = section->next;
        }
    }
    if (section == NULL) {
        move_gain = 0.03;
        move_interval = 100;
    } else {
        xmlChar *prop = NULL;
        if ((prop = xmlGetProp(section, (xmlChar *)"gain")) != NULL) {
            move_gain = strtod((char *)prop, NULL);
            xmlFree(prop);
            prop = NULL;
        } else
            move_gain = 0.03;
        if ((prop = xmlGetProp(section, (xmlChar *)"interval")) != NULL) {
            move_interval = atoi((char *)prop);
            xmlFree(prop);
            prop = NULL;
        } else
            move_interval = 100;
        if ((prop = xmlGetProp(section, (xmlChar *)"target")) != NULL) {
            target = strtod((char *)prop, NULL);
            xmlFree(prop);
            prop = NULL;
        } else
            target = 0.44;
        if ((prop = xmlGetProp(section, BAD_CAST "init_theta_bar")) != NULL){
            double initTheta = strtod((char *)prop, NULL);
            xmlFree(prop);
            prop = NULL;
            for (size_t i = 0; i < nparams; ++i) {
                theta_bars[i] = initTheta;
            }
        }
        // check theta min and max
        int i;
        for (xmlNodePtr thetaNode = section->children; thetaNode != NULL;
        		thetaNode = thetaNode->next) {
        	if (xmlStrcmp(thetaNode->name, BAD_CAST "theta"))
        		continue;
        	if ((prop = xmlNodeGetContent(thetaNode)) == NULL)
        		continue;
        	i = atoi((char *)prop);
        	xmlFree(prop);
        	prop = NULL;
        	if (i < 0 || i >= nparams)
        		continue;
        	if ((prop = xmlGetProp(thetaNode, BAD_CAST "max")) != NULL){
        		theta_maxs[i] = strtod((char *)prop, NULL);
        		xmlFree(prop);
        		prop = NULL;
        	}
        	if ((prop = xmlGetProp(thetaNode, BAD_CAST "min")) != NULL){
        		theta_mins[i] = strtod((char *)prop, NULL);
        		xmlFree(prop);
        		prop = NULL;
        	}
        }

    }
}

template<class Problem>
inline double feedbackMove<Problem>::get_score()
{
    return energy;
}

template<class Problem>
double feedbackMove<Problem>::propose()
{
    index++;
    index %= nparams;
    if (index == 0) {
        ++sweep;
        if (sweep % move_interval == 0) {
            collectMoveStats();
            move_control();
        }

    }

    prev_energy = energy;
    // generate theta from theta bar here and pass only theta
    // to the problem so problem doesn't have to have random
    // number generator
    double uniform = 2.0 * rnd.random() - 1.0;
    double theta;
    if (uniform >= 0.)
        theta = -1 * theta_bars[index] * log(abs(uniform));
    else
        theta = theta_bars[index] * log(abs(uniform));
    problem.generateMove(index, theta);
    ++ moves[index]; // prefix increment has lower precedence than [], right?
    energy = problem.get_score();
    return (energy - prev_energy);
}

template<class Problem>
void feedbackMove<Problem>::accept()
{
    ++ success[index];
}

template<class Problem>
feedbackMove<Problem>::~feedbackMove()
{
    delete[] success;
    delete[] moves;
    delete[] theta_bars;
    delete[] theta_mins;
    delete[] theta_maxs;
}

template<class Problem>
void feedbackMove<Problem>::reject()
{
    problem.restoreMove(index);
    energy = prev_energy;
}

template<class Problem>
void feedbackMove<Problem>::move_control()
{
    debugOut << sweep;
    for (int i = 0; i < nparams; ++i) {
        double acc_ratio = (double) success[i] / (double) moves[i];
        double x = log(theta_bars[i]);
        debugOut << "\t" << acc_ratio;
        x += move_gain * (acc_ratio - target);
        theta_bars[i] = exp(x);
        if (theta_bars[i] < theta_mins[i])
            theta_bars[i] = theta_mins[i];
        if (theta_bars[i] > theta_maxs[i])
            theta_bars[i] = theta_maxs[i];
        success[i] = 0;
        moves[i] = 0;
    }
    if (!debugOut.isIgnore()){
        for (int i = 0; i < nparams; ++i) {
            debugOut << "\t" << theta_bars[i];
        }
    }

    debugOut << endl;
}

template<class Problem>
void feedbackMove<Problem>::writeState(xmlNodePtr docroot) const
{
    xmlNodePtr moveNode = xmlNewChild(docroot, NULL, BAD_CAST"moveSize", NULL);
    char *paramNumString = NULL;
    char *thetaString = NULL;
    xmlNodePtr paramIter = NULL;
    int i;
    for (i = 0; i < nparams; ++i) {
        asprintf(&paramNumString, "%d", i);
        paramIter = xmlNewChild(moveNode, NULL, BAD_CAST "param", 
                BAD_CAST paramNumString);
        free(paramNumString);
        asprintf(&thetaString, "%.15g", theta_bars[i]);
        xmlNewProp(paramIter, BAD_CAST "theta", BAD_CAST thetaString);
        free(thetaString);
    }
}

template<class Problem>
void feedbackMove<Problem>::readState(xmlNodePtr docroot)
{
    xmlNodePtr moveNode = getSectionByName(docroot, "moveSize");
    if (!moveNode) {
        cerr << "Warning: fail to read moveSize from state file" << endl;
        return;
    }
    int i;
    xmlNodePtr paramIter;
    char *paramNumString;
    for (paramIter = moveNode->children; paramIter != NULL; 
            paramIter = paramIter->next) {
        if (xmlStrcmp(paramIter->name, BAD_CAST "param"))
            continue;
        paramNumString = (char *)xmlNodeGetContent(paramIter);
        sscanf(paramNumString,"%d",&i);
        xmlFree(paramNumString);
        if (i < nparams) {
            theta_bars[i] = getPropDouble(paramIter,"theta");
        }
    }
        
    
}