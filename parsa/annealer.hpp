/*
 * annealer.hpp
 *
 *  Created on: Jan 16, 2013
 *      Author: zhlou
 *  The implementation of the annealer template
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <stdexcept>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/format.hpp>

using namespace std;
/*
 * When the constructor is called, all other parts should already be set up.
 */

template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
annealer<Problem, Schedule, FrozenCnd, Move>::annealer(Problem &problem,
        unirandom& in_rand, typename Schedule::Param scheParam,
        typename FrozenCnd::Param frozenParam,
        const ptree &root) :
        problem(problem),
        rand(in_rand)
{
    initState(root);
    cooling = new Schedule(scheParam);
    move = new Move<Problem>(problem, in_rand, root);
    frozen = new FrozenCnd(frozenParam);
    state.energy = move->get_score();
    tlaps = -1;
}
/*
 * This construtor is used only by derived class that have there own way
 * of generating cooling schedule and move scheme
 */

template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
annealer<Problem, Schedule, FrozenCnd, Move>::annealer(Problem &problem,
        unirandom& in_rand,
        const ptree &root) :
        problem(problem),
        rand(in_rand)
{
    initState(root);
    cooling = NULL;
    move = NULL;
    frozen = NULL;
    tlaps = -1;
}

template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
annealer<Problem, Schedule, FrozenCnd, Move>::~annealer()
{
    delete cooling;
    delete move;
    delete frozen;
}

template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
double annealer<Problem, Schedule, FrozenCnd, Move>::loop()
{
    clock_t start = clock();
    if (!is_init)
        initMoves();
    bool accepted;
    // state.step_cnt = 0;
    do { // while not frozen
//        cooling->resetSegmentStats();

        //do { // while in segment
            accepted = step();

            cooling->updateStep(accepted, state);
            frozen->updateStep(accepted, state);
            state.s = cooling->updateS(state);
            ++ (state.step_cnt);
            updateStats(state);
//            cooling->updateStats(state);
//            move->updateStats(state);
        //} while (cooling->inSegment(state));

//        updateSegment(state);
    } while (!frozen->frozen(state));

    clock_t end = clock();
    tlaps = (double)(end - start)/CLOCKS_PER_SEC;

    //std::cout << "Annealing stopped at s = " << state.s << std::endl
    //        << "Total steps is " << state.step_cnt << std::endl;
    return move->get_score();
}

template<class Problem, class Schedule, class FrozenCnd, template<class > class Move>
double annealer<Problem, Schedule, FrozenCnd, Move>::fixedTMoves(double S, long steps)
{
    double prev_s = state.s;
    state.s = S;
    for (int i = 0; i < steps; ++i) {
        step();
    }
    state.s = prev_s;
    return move->get_score();
}

template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
double annealer<Problem, Schedule, FrozenCnd, Move>::initMoves()
{
    fixedTMoves(state.s, initLoop);
    bool accepted;
    int nAccept = 0;
    onePassMeanVar energyMeanVar;
    for (state.step_cnt = 0; state.step_cnt < initLoop; ++ state.step_cnt) {
        nAccept += (accepted = step());
        energyMeanVar.update(state.energy);
        cooling->updateInitStep(accepted, state);
    }
    cooling->initStats(state);
    initMean = energyMeanVar.getMean();
    initVar = energyMeanVar.getVar();
    initAccRatio = double(nAccept)/initLoop;

    is_init = true;
    return move->get_score();

}

template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
double annealer<Problem, Schedule, FrozenCnd, Move>::initMovesOnly()
{
    fixedTMoves(state.s, initLoop);
    return move->get_score();

}

template<class Problem, class Schedule, class FrozenCnd, template<class > class Move>
void annealer<Problem, Schedule, FrozenCnd, Move>::writeResultData(ptree &result) const
{
    std::ostringstream s_energy, s_step, s_time;
    s_energy << state.energy;
    s_step << state.step_cnt;
    s_time << tlaps;
    result.put("<xmlattr>.final_energy",s_energy.str());
    result.put("<xmlattr>.max_count", s_step.str());
    result.put("<xmlattr>.time", s_time.str());

}

template<class Problem, class Schedule, class FrozenCnd, template<class > class Move>
void annealer<Problem, Schedule, FrozenCnd, Move>::writeMethodText(ptree &method) const
{
    method.put("<xmlattr>.cooling-schedule", std::string(Schedule::name));
    method.put("<xmlattr>.move-generation", std::string(Move<Problem>::name));
}

template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
void annealer<Problem, Schedule, FrozenCnd, Move>::ptreeGetResult(ptree &pt)
{
    std::ostringstream s_energy, s_step, s_time;
    s_energy << state.energy;
    s_step << state.step_cnt;
    s_time << tlaps;
    pt.put("annealing_result.final_energy",s_energy.str());
    pt.put("annealing_result.max_count",s_step.str());
    pt.put("annealing_result.time",s_time.str());
    pt.put("annealing_method.cooling-schedule",std::string(Schedule::name));
    pt.put("annealing_method.move-generation",std::string(Move<Problem>::name));
}

template<class Problem, class Schedule, class FrozenCnd, template<class > class Move>
void annealer<Problem, Schedule, FrozenCnd, Move>::writeResult(ptree &root) const
{
    ptree result, method;
    writeResultData(result);
    writeMethodText(method);
    root.put_child("annealing_result", result);
    root.put_child("annealing_method", method);

}

template<class Problem, class Schedule, class FrozenCnd, template<class > class Move>
void annealer<Problem, Schedule, FrozenCnd, Move>::initState(const ptree &root)
{
    const ptree &sec_attr = root.get_child("annealer_input.<xmlattr>");
    initS = 1.0 / sec_attr.get<double>("init_T");
    initLoop = sec_attr.get<int>("init_loop");
    state.s = initS;
    is_init = false;
    state.step_cnt = 0;
}

template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
bool annealer<Problem, Schedule, FrozenCnd, Move>::step()
{
    double delta, crit, ran_n;
    bool flag;
    delta = move->propose(state.s);
    state.proposed = state.energy + delta;
    // cout << state.energy + delta << "@" << state.step_cnt << endl;
    crit = std::exp(-state.s * delta);
    ran_n = rand.random();
    if ((delta <= 0.0) || crit > ran_n) {
        move->accept();
        state.energy += delta;

        flag = true;
    } else {
        move->reject();
        flag = false;
    }
    debugOut << state.step_cnt << " " << state.s << " "
             << state.energy << " " << state.proposed << endl;
    return flag;
}

/*
template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
void annealer<Problem, Schedule, FrozenCnd, Move>::saveUnifiedInitState(const char* filename)
{
    char * stateName;
    char * xmlName;
    asprintf(&stateName, "%s.state", filename);
    asprintf(&xmlName, "%s.xml", filename);

    int buf_size = problem.getStateSize();
    char * state_buf = new char[buf_size];
    problem.serialize(state_buf);
    ofstream stateFile(stateName, ios::out | ios::binary);
    stateFile.write(state_buf, buf_size);
    stateFile.close();
    delete[] state_buf;
    free(stateName);

    xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
    xmlNodePtr docroot = xmlNewNode(NULL, BAD_CAST "root");
    xmlDocSetRootElement(doc, docroot);

    xmlNodePtr initStatNode = xmlNewChild(docroot, NULL, BAD_CAST"initStat", NULL);
    char * initMeanString;
    char * initVarString;
    char * initAccRatioString;
    asprintf(&initMeanString, "%.16g", initMean);
    asprintf(&initVarString, "%.16g", initVar);
    asprintf(&initAccRatioString, "%.16g", initAccRatio);
    xmlNewProp(initStatNode, BAD_CAST "initMean", BAD_CAST initMeanString);
    xmlNewProp(initStatNode, BAD_CAST "initVar", BAD_CAST initVarString);
    xmlNewProp(initStatNode, BAD_CAST "initAccRatio", BAD_CAST initAccRatioString);

    move->writeState(docroot);

    xmlSaveFormatFile(xmlName, doc, 1);
    xmlFreeDoc(doc);
    free(initMeanString);
    free(initVarString);
    free(initAccRatioString);
    free(xmlName);
}
*/

template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
void annealer<Problem, Schedule, FrozenCnd, Move>
        ::saveUnifiedInitState(const std::string &statePrefix)
{
    std::string stateName = statePrefix + ".state";
    std::string xmlName = statePrefix + ".xml";
    int buf_size = problem.getStateSize();
    char * state_buf = new char[buf_size];
    problem.serialize(state_buf);
    ofstream stateFile(stateName.c_str(), ios::out | ios::binary);
    stateFile.write(state_buf, buf_size);
    stateFile.close();
    delete[] state_buf;

    ptree pt;
    ptree sec_attr;
    sec_attr.put("initMean", (boost::format("%.16g") % initMean).str());
    sec_attr.put("initVar", (boost::format("%.16g") % initVar).str());
    sec_attr.put("initAccRatio", (boost::format("%.16g") % initAccRatio).str());
    pt.put_child("root.initStat.<xmlattr>", sec_attr);
    ptree &root = pt.get_child("root");
    move->writeState(root);
    boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
    write_xml(xmlName, pt, std::locale(), settings);

}

/*
template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
double annealer<Problem, Schedule, FrozenCnd, Move>::readUnifiedInitState(const char* filename)
{
    char * stateName;
    char * xmlName;
    asprintf(&stateName, "%s.state", filename);
    asprintf(&xmlName, "%s.xml", filename);

    int bufSize = problem.getStateSize();
    char * stateBuf = new char[bufSize];
    ifstream stateFile(stateName, ios::in | ios::binary);
    if (!stateFile) {
        throw std::runtime_error("Fail to open state file");
    }
    stateFile.read(stateBuf, bufSize);
    stateFile.close();
    problem.deserialize(stateBuf);
    state.energy = move->forceUpdateEnergy();

    delete[] stateBuf;
    free(stateName);

    xmlDocPtr doc = xmlParseFile(xmlName);
    xmlNodePtr docroot = xmlDocGetRootElement(doc);

    xmlNodePtr initStateNode = getSectionByName(docroot, "initStat");
    initMean = getPropDouble(initStateNode, "initMean");
    initVar = getPropDouble(initStateNode, "initVar");
    initAccRatio = getPropDouble(initStateNode, "initAccRatio");

    move->readState(docroot);

    xmlFreeDoc(doc);
    free(xmlName);

    cooling->initStats(initMean, initVar, initAccRatio, state);
    is_init = true;
    return state.energy;
}
*/

template<class Problem, class Schedule, class FrozenCnd, template<class> class Move>
double annealer<Problem, Schedule, FrozenCnd, Move>::readUnifiedInitState(const std::string &prefix)
{
    std::string stateName = prefix + ".state";
    std::string xmlName = prefix + ".xml";

    int bufSize = problem.getStateSize();
    char * stateBuf = new char[bufSize];
    ifstream stateFile(stateName.c_str(), ios::in | ios::binary);
    if (!stateFile) {
        throw std::runtime_error("Fail to open state file");
    }
    stateFile.read(stateBuf, bufSize);
    stateFile.close();
    problem.deserialize(stateBuf);
    state.energy = move->forceUpdateEnergy();

    delete[] stateBuf;

    ptree pt;
    read_xml(xmlName, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &root = pt.begin()->second;

    const ptree &sec_attr = root.get_child("initStat.<xmlattr>");
    initMean = sec_attr.get<double>("initMean");
    initVar = sec_attr.get<double>("initVar");
    initAccRatio = sec_attr.get<double>("initAccRatio");

    move->readState(root);

    cooling->initStats(initMean, initVar, initAccRatio, state);
    is_init = true;
    return state.energy;
}
