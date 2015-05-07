/*********************************************************************************
*                                                                                *
*     organism.cpp                                                               *
*                                                                                *
*     An organism is defined here as a collection of nuclei. This is the master  *
*     class which holds most of the data. It creates and array of nuclei         *
*     dynamicaly based on the data present to be scored against.                 *
*                                                                                *
*********************************************************************************/


#include "organism.h"
#include <boost/foreach.hpp>
#include <limits>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/lexical_cast.hpp>

#define foreach_ BOOST_FOREACH
#define to_ boost::lexical_cast
#define to_string_ boost::lexical_cast<string>


/*    Constructors   */

Organism::Organism() :
    //mode(mode_ptr(new Mode)),
    distances(distances_ptr(new DistanceContainer)),
    master_tfs(tfs_ptr(new TFContainer)),
    master_genes(genes_ptr(new GeneContainer)),
    tfdata(table_ptr(new DataTable<double>)),
    ratedata(table_ptr(new DataTable<double>)),
    promoters(promoters_ptr(new PromoterContainer)),
    scale_factors(scale_factors_ptr(new ScaleFactorContainer)),
    coeffects(coeffects_ptr(new CoeffectContainer)),
    score_class(score_ptr(new Score)),
    coops(coops_ptr(new CooperativityContainer)),
    competition(competition_ptr(new Competition)),
    nuclei(nuclei_ptr(new Nuclei))
{
  move_count = 0;
  test_int = 0;
}

Organism::Organism(string fname, string section) :
    //mode(mode_ptr(new Mode)),
    distances(distances_ptr(new DistanceContainer)),
    master_tfs(tfs_ptr(new TFContainer)),
    master_genes(genes_ptr(new GeneContainer)),
    tfdata(table_ptr(new DataTable<double>)),
    ratedata(table_ptr(new DataTable<double>)),
    promoters(promoters_ptr(new PromoterContainer)),
    scale_factors(scale_factors_ptr(new ScaleFactorContainer)),
    coeffects(coeffects_ptr(new CoeffectContainer)),
    score_class(score_ptr(new Score)),
    coops(coops_ptr(new CooperativityContainer)),
    competition(competition_ptr(new Competition)),
    nuclei(nuclei_ptr(new Nuclei))
{
  move_count = 0;
  fstream infile(fname.c_str());
  ptree pt;
  read_xml(infile, pt, boost::property_tree::xml_parser::trim_whitespace);

  ptree& root_node  = pt.get_child("Root");
  ptree& mode_node  = root_node.get_child("Mode");
  ptree& input_node = root_node.get_child(section);
  mode_ptr m(new Mode(fname, mode_node));
  mode = m;
  initialize(input_node);
}

Organism::Organism(ptree& pt, mode_ptr m) :
    //mode(mode_ptr(new Mode)),
    distances(distances_ptr(new DistanceContainer)),
    master_tfs(tfs_ptr(new TFContainer)),
    master_genes(genes_ptr(new GeneContainer)),
    tfdata(table_ptr(new DataTable<double>)),
    ratedata(table_ptr(new DataTable<double>)),
    promoters(promoters_ptr(new PromoterContainer)),
    scale_factors(scale_factors_ptr(new ScaleFactorContainer)),
    coeffects(coeffects_ptr(new CoeffectContainer)),
    score_class(score_ptr(new Score)),
    coops(coops_ptr(new CooperativityContainer)),
    competition(competition_ptr(new Competition)),
    nuclei(nuclei_ptr(new Nuclei))
{
  move_count = 0;
  mode = m;
  initialize(pt);
}

void Organism::initialize(ptree& pt)
{
  move_count = 0;
  test_int = 0;

  //mode->read(pt);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized the problem..." << endl;


  if (mode->getCompetition())
    competition->set(pt, mode);

  distances->setMode(mode);
  distances->add(pt);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized distance functions" << endl;

  promoters->setMode(mode);
  promoters->read(pt);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized promoter functions" << endl;

  master_tfs->add(pt, mode);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized transcription factors" << endl;

  scale_factors->setMode(mode);
  scale_factors->read(pt);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized scale factors" << endl;

  master_genes->setPromoters(promoters);
  master_genes->setScaleFactors(scale_factors);
  master_genes->read(pt);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized genes" << endl;

  ratedata->read(pt, "RateData");
  if (mode->getVerbose() >= 2)
    cerr << "Initialized rate data" << endl;

  tfdata->read(pt, "TFData");
  if (mode->getVerbose() >= 2)
    cerr << "Initialized TF data" << endl;

  coops->setMode(mode);
  coops->read(pt, distances);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized cooperativity" << endl;

  coeffects->setMode(mode);
  coeffects->read(pt,distances);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized coactivation and corepression" << endl;

  master_tfs->setCoops(coops);
  master_tfs->setCoeffects(coeffects);

  competition->getParameters(params);
  distances->getParameters(params);
  master_tfs->getParameters(params);
  promoters->getParameters(params);
  scale_factors->getParameters(params);
  coops->getParameters(params);
  coeffects->getParameters(params);
  master_genes->getParameters(params);

  competition->getAllParameters(all_params);
  distances->getAllParameters(all_params);
  master_tfs->getAllParameters(all_params);
  promoters->getAllParameters(all_params);
  scale_factors->getAllParameters(all_params);
  coops->getAllParameters(all_params);
  coeffects->getAllParameters(all_params);

  if (mode->getVerbose() >= 2)
  {
    cerr << "Initialized parameters" << endl;
    cerr << endl;
  }

  if (mode->getVerbose() >= 1)
    printParameters(cerr);

  populate_nuclei();
  if (mode->getVerbose() >= 2)
    cerr << "Created Nuclei" << endl;

  score_class->set(this);

  score();

  setMoves();
}

/*    Setters   */

void Organism::populate_nuclei()
{
  nuclei->setParent(this);

  ids        = ratedata->getNames("ID");
  int nnuc   = ids.size();

  for (int i = 0; i<nnuc; i++)
  {
    string& id = ids[i];
    nuclei->addNuc(id);
  }

  nuclei->create();

  if (nnuc == 0)
  {
    stringstream err;
    err << "ERROR: populate_nuclei() did not find any nuclei!" << endl;
    error(err.str());
  }
}

double* Organism::getPrediction(Gene& gene, string& id)
{
  vector<string>& tmp_ids = nuclei->getIDs();
  int nids = tmp_ids.size();
  for (int j=0; j<nids; j++)
  {
    if (tmp_ids[j] == id)
      return &(nuclei->getRate(gene, id));
  }

  stringstream err;
  err << "ERROR: could not get prediction for " << gene.getName() << " at " << id << endl;
  error(err.str());
  return(0);
}

double* Organism::getData(Gene& gene, string& id)
{
  const string& gname = gene.getName();

  return &(ratedata->getDataPoint("gene", gname, "ID", id));
}

bindings_ptr Organism::getBindings()
{
  return nuclei->getBindings();
}

/*    Methods   */

/* scoring is actually going to be somewhat tricky in the future. Calculations
over nuclei need to be vectorized, which means that nuclei loses information about
IDs. However, at the same time we don't necessarily want nuclei and ratedata
to hold everything in the same order. This means we will need to create some map
of rate to model in scoring */

void Organism::score()
{
  score_out = score_class->getScore();
}

void Organism::printScore(ostream& os)
{
  score_class->print(os);
}

void Organism::checkScale(ostream& os)
{
  score_class->checkScale(os);
}


/*    Output    */

void Organism::write(string node, ptree& pt)
{
  ptree& output = pt.add(node,"");
  if (mode->getVerbose() >= 2)
    cerr << "writing output..." << endl;
  if (mode->getCompetition())
    competition->write(output);
  distances->write(output);
  if (mode->getVerbose() >= 2)
    cerr << "wrote distances" << endl;
  promoters->write(output);
  if (mode->getVerbose() >= 2)
    cerr << "wrote promoters" << endl;
  master_tfs->write(output);
  if (mode->getVerbose() >= 2)
    cerr << "wrote master_tfs" << endl;
  coops->write(output);
  if (mode->getVerbose() >= 2)
    cerr << "wrote coops" << endl;
  coeffects->write(output);
  if (mode->getVerbose() >= 2)
    cerr << "wrote coeffects" << endl;
  master_genes->write(output);
  if (mode->getVerbose() >= 2)
    cerr << "wrote master_genes" << endl;
  scale_factors->write(output);
  if (mode->getVerbose() >= 2)
    cerr << "wrote scale_factors" << endl;
  ratedata->write(output, "RateData", mode->getPrecision());
  if (mode->getVerbose() >= 2)
    cerr << "wrote ratedata" << endl;
  tfdata->write(output, "TFData", mode->getPrecision());
  if (mode->getVerbose() >= 2)
    cerr << "wrote tfdata" << endl;
}

void Organism::printRate(ostream& os, bool invert)
{
  int p = mode->getPrecision();   // the precision to print with
  int w = p+7; // set the minimum width, must be 6+precision to ensure columns dont merge with scientific notation
  int namew=0;
  int ngenes = master_genes->size();
  for (int i=0; i<ngenes; i++)
    namew = max(namew,(int) master_genes->getGene(i).getName().size());

  namew++;

  if (!invert)
  {
    w = max(w,namew);
    os << setw(w) << setprecision(p) << "id";
    for (int i=0; i<ngenes; i++)
      os << setw(w) << master_genes->getGene(i).getName();
    os << endl;

    vector<string>& IDs = nuclei->getIDs();
    int nids            = IDs.size();
    for (int j=0; j<nids; j++)
    {
      os << setw(w) << IDs[j];
      for (int k=0; k<ngenes; k++)
      {
        Gene& gene = master_genes->getGene(k);
        vector<double>& rate = nuclei->getRate(gene);
        os << setw(w) << rate[j];
      }
      os << endl;
    }
  }
  else
  {
    os << setw(namew) << setprecision(p) << "id";
    vector<string>& IDs = nuclei->getIDs();
    int nids            = IDs.size();
    for (int j=0; j<nids; j++)
      os << setw(w) << IDs[j];
    os << endl;

    for (int k=0; k<ngenes; k++)
    {
      Gene& gene = master_genes->getGene(k);
      os << setw(namew) << gene.getName();
      vector<string>& IDs = nuclei->getIDs();
      int nids         = IDs.size();
      vector<double>& rate = nuclei->getRate(gene);
      for (int j=0; j<nids; j++)
        os << setw(w) << rate[j];
      os << endl;
    }
  }
}


void Organism::printR2D(ostream& os)
{
  if (mode->getCompetition())
  {
    int ngenes = master_genes->size();
    for (int j=0; j<ngenes; j++)
    {
      Gene& gene = master_genes->getGene(j);
      nuclei->printR2D(gene, os);
    }
  }
  else
    error("Cannot print 2D rate with competition mode off");

}

void Organism::printN2D(ostream& os)
{
  if (mode->getCompetition())
  {
    int ngenes = master_genes->size();
    for (int j=0; j<ngenes; j++)
    {
      Gene& gene = master_genes->getGene(j);
      nuclei->printN2D(gene, os);
    }
  }
  else
    error("Cannot print 2D rate with competition mode off");
}


void Organism::printRateData(ostream& os, bool invert)
{
  int p = mode->getPrecision();   // the precision to print with
  int w = p+7; // set the minimum width, must be 6+precision to ensure columns dont merge with scientific notation
  int namew=0;
  int ngenes = master_genes->size();
  for (int i=0; i<ngenes; i++)
    namew = max(namew,(int) master_genes->getGene(i).getName().size());

  namew++;


  if (!invert)
  {
    w = max(w,namew);
    os << setw(w) << setprecision(p) << "id";
    for (int i=0; i<ngenes; i++)
    {
      Gene& gene = master_genes->getGene(i);
      if (gene.getInclude())
        os << setw(w) << gene.getName();
    }
    os << endl;

    vector<string>& IDs = nuclei->getIDs();
    int nids         = IDs.size();
    for (int j=0; j<nids; j++)
    {
      string& id = IDs[j];
      os << setw(w) << id;
      for (int k=0; k<ngenes; k++)
      {
        Gene& gene = master_genes->getGene(k);
        if (!gene.getInclude()) continue;
        scale_factor_ptr scale = gene.getScale();
        double datapoint = ratedata->getDataPoint("gene",gene.getName(),"ID",id);
        datapoint = scale->scale(datapoint);
        os << setw(w) << datapoint;
      }
      os << endl;
    }
  }
  else
  {
    os << setw(namew) << setprecision(p) << "id";

    vector<string>& IDs = nuclei->getIDs();
    int nids         = IDs.size();
    for (int j=0; j<nids; j++)
      os << setw(w) << IDs[j];

    os << endl;

    for (int k=0; k<ngenes; k++)
    {
      Gene& gene = master_genes->getGene(k);
      if (!gene.getInclude()) continue;
      scale_factor_ptr scale = gene.getScale();
      os << setw(namew) << gene.getName();

      vector<string>& IDs = nuclei->getIDs();
      int nids         = IDs.size();
      for (int j=0; j<nids; j++)
      {
        string& id = IDs[j];
        double datapoint = ratedata->getDataPoint("gene",gene.getName(),"ID",id);
        datapoint = scale->scale(datapoint);
        os << setw(w) << datapoint;
      }
      os << endl;
    }
  }
}

void Organism::printSites(ostream& os)
{
  nuclei->printSites(os);
}

void Organism::printSites(Gene& gene, ostream& os)
{
  nuclei->printSites(gene, os);
}

void Organism::printSites(TF& tf, ostream& os)
{
  nuclei->printSites(tf, os);
}

void Organism::printSites(Gene& gene, TF& tf, ostream& os)
{
  nuclei->printSites(gene, tf, os);
}

void Organism::printScores(Gene& gene, ostream& os)
{
  nuclei->printScores(gene, os);
}

void Organism::printSubgroups(Gene& gene, ostream& os)
{
  nuclei->printSubgroups(gene, os);
}

void Organism::printOccupancy(Gene& gene, ostream& os, bool invert)
{
  nuclei->printOccupancy(gene, os, invert);
}

void Organism::printModeOccupancy(Gene& gene, ostream& os, bool invert)
{
  nuclei->printModeOccupancy(gene, os, invert);
}

void Organism::printEffectiveOccupancy(Gene& gene, ostream& os, bool invert)
{
  nuclei->printEffectiveOccupancy(gene, os, invert);
}

/*    Move    */

/* scramble now changes the input node as well as the parameter value, so that
the node can just be reprinted without reprinting anything else */

void Organism::scramble()
{
  boost::minstd_rand baseGen(getpid());
  boost::uniform_real<> uniDblUnit(0,1);
  boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > uniDblGen(baseGen, uniDblUnit);

  uniDblGen();
  int nparams = params.size();
  for (int i=0; i<nparams; i++)
    params[i]->scramble(uniDblGen());
}

void Organism::permute(string& table, string& by)
{
  if (table == string("RateData"))
    ratedata->permute(by, mode->getPrecision());
  else if (table == string("TFData"))
    tfdata->permute(by, mode->getPrecision());
  else
  {
    stringstream err;
    err << "ERROR: no data table with name " << table << endl;
    error(err.str());
  }
}

/***************  Move Generation  **********************************************/

void Organism::setMoves()
{
  /* First, we set up maps for every type of move and restore we might make.
  This is really just here to make it cleaner and not have a bunch of if blocks
  for move types */
  move_map[string("ResetAll"        )] = &Organism::ResetAll;
  move_map[string("PWM"             )] = &Organism::movePWM;
  move_map[string("Scores"          )] = &Organism::moveScores;
  move_map[string("Sites"           )] = &Organism::moveSites;
  move_map[string("Lambda"          )] = &Organism::moveLambda;
  move_map[string("Kmax"            )] = &Organism::moveKmax;
  move_map[string("CoopD"           )] = &Organism::moveCoopD;
  move_map[string("Kcoop"           )] = &Organism::moveKcoop;
  move_map[string("Coef"            )] = &Organism::moveCoef;
  move_map[string("Quenching"       )] = &Organism::moveQuenching;
  move_map[string("QuenchingCoef"   )] = &Organism::moveQuenchingCoef;
  move_map[string("Coeffect"        )] = &Organism::moveCoeffect;
  move_map[string("CoeffectEff"     )] = &Organism::moveCoeffectEff;
  move_map[string("Promoter"        )] = &Organism::movePromoter;
  move_map[string("Null"            )] = &Organism::null_function;

  restore_map[string("ResetAll"     )] = &Organism::ResetAll;
  restore_map[string("PWM"          )] = &Organism::restorePWM;
  restore_map[string("Scores"       )] = &Organism::restoreScores;
  restore_map[string("Sites"        )] = &Organism::restoreSites;
  restore_map[string("Lambda"       )] = &Organism::restoreLambda;
  restore_map[string("Kmax"         )] = &Organism::restoreKmax;
  restore_map[string("CoopD"        )] = &Organism::restoreCoopD;
  restore_map[string("Kcoop"        )] = &Organism::restoreKcoop;
  restore_map[string("Coef"         )] = &Organism::restoreCoef;
  restore_map[string("Quenching"    )] = &Organism::restoreQuenching;
  restore_map[string("QuenchingCoef")] = &Organism::restoreQuenchingCoef;
  restore_map[string("Coeffect"     )] = &Organism::restoreCoeffect;
  restore_map[string("CoeffectEff"  )] = &Organism::restoreCoeffectEff;
  restore_map[string("Promoter"     )] = &Organism::movePromoter;
  restore_map[string("Null"         )] = &Organism::null_function;

  string move;
  int nparams = params.size();
  for (int i=0; i<nparams; i++)
  {
    move = params[i]->getMove();

    // make sure we have these functions available
    if (move_map.find(move) == move_map.end())
      error("setMoves() Could not find move function with name " + move + "for parameter " + params[i]->getParamName());

    if (restore_map.find(move) == restore_map.end())
      error("setMoves() Could not find restore function with name " + move + "for parameter " + params[i]->getParamName());

    moves.push_back(move_map[move]);
    restores.push_back(restore_map[move]);
  }
}


/*
Now we define the move and restore functions. I have tried to order these
according to how much is reset in each case. I thought this might be easier to
maintain if I set a "reset level" for each parameter, however there are no clear
rules for what needs to be reset depending on what parameter we tweak. Instead,
this requires pretty good knowledge of the transcription equations themselves, 
so be careful when adding to these. 
 
There are four types of functions:
 
save     - saves old data structures
update   - implements any preprocessing of interactions and memory managements
calc     - preprocessing is unnecessary, so simply rerun the calculation
restore  - copies over old data structures 
*/

/* The slowest way to reset everything, this function deletes all calculated objects
and repopulates them from scratch! */
void Organism::Recalculate()
{
  nuclei = nuclei_ptr(new Nuclei);
  
  master_tfs->setCoops(coops);
  master_tfs->setCoeffects(coeffects);

  params.clear();
  all_params.clear();
  
  competition->getParameters(params);
  distances->getParameters(params);
  master_tfs->getParameters(params);
  promoters->getParameters(params);
  scale_factors->getParameters(params);
  coops->getParameters(params);
  coeffects->getParameters(params);
  master_genes->getParameters(params);

  competition->getAllParameters(all_params);
  distances->getAllParameters(all_params);
  master_tfs->getAllParameters(all_params);
  promoters->getAllParameters(all_params);
  scale_factors->getAllParameters(all_params);
  coops->getAllParameters(all_params);
  coeffects->getAllParameters(all_params);
  
  populate_nuclei();
  
  score_class->set(this);

  score();

  setMoves();
}

/* clear all data and reset from scratch. If you suspect a move function
is not working, this is the easiest way to check, however it should generally
be avoided as it will be very slow */

void Organism::ResetAll(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Reseting everything" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->updateScores(gene);
    nuclei->updateSites(gene);
    nuclei->updateSubgroups(gene);
    nuclei->updateCoeffects(gene);
    nuclei->updateQuenching(gene);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);

  }
}

/* this will be necessary if something changes the way we score sequence, for instance
tweaking pwms. In general we believe this is a bad idea, but it could be useful
for comparing to other results (segal) or in very careful applications. If you
decide to use this feature, be ready to defend your decision! */
void Organism::moveScores(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving sequence scores" << endl;

  TF& tf = master_tfs->getTF(params[idx]->getTFName());

  int ngenes  = master_genes->size();

#ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    nuclei->saveSites(gene,tf);
    nuclei->saveScores(gene,tf);
    nuclei->saveSubgroups(gene);
    nuclei->saveCoeffects(gene);
    nuclei->saveQuenching(gene);

    nuclei->updateScores(gene,tf);
    nuclei->updateSites(gene,tf);
    nuclei->updateSubgroups(gene);
    nuclei->updateCoeffects(gene);
    nuclei->updateQuenching(gene);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

void Organism::restoreScores(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring sequence scores" << endl;

  TF& tf = master_tfs->getTF(params[idx]->getTFName());

  int ngenes  = master_genes->size();

#ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreScores(gene,tf);
    nuclei->restoreSites(gene,tf);
    nuclei->restoreAllOccupancy(gene);
    nuclei->restoreSubgroups(gene);
    nuclei->restoreCoeffects(gene);
    nuclei->restoreQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

/* this will be necessary if something changes the way we score sequence, for instance
tweaking pwms. In general we believe this is a bad idea, but it could be useful
for comparing to other results (segal) or in very careful applications. If you
decide to use this feature, be ready to defend your decision! */
void Organism::movePWM(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving PWM" << endl;

  TF& tf = master_tfs->getTF(params[idx]->getTFName());
  int ngenes  = master_genes->size();

#ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    nuclei->saveSites(gene, tf);
    nuclei->saveScores(gene, tf);
    nuclei->saveSubgroups(gene);
    nuclei->saveCoeffects(gene);
    nuclei->saveQuenching(gene);

    nuclei->updateScores(gene,tf);
    nuclei->updateSites(gene,tf);
    nuclei->updateSubgroups(gene);
    nuclei->updateCoeffects(gene);
    nuclei->updateQuenching(gene);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

void Organism::restorePWM(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring PWM" << endl;

  TF& tf = master_tfs->getTF(params[idx]->getTFName());
  int ngenes  = master_genes->size();

#ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreScores(gene,tf);
    nuclei->restoreSites(gene,tf);
    nuclei->restoreAllOccupancy(gene);
    nuclei->restoreSubgroups(gene);
    nuclei->restoreCoeffects(gene);
    nuclei->restoreQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}


/* if thresholds are changed, the sites will need to be repopulated for the tf
changed */
void Organism::moveSites(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving Sites" << endl;

  TF& tf = master_tfs->getTF(params[idx]->getTFName());


  int ngenes  = master_genes->size();

#ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    nuclei->saveSites(gene, tf);
    nuclei->saveSubgroups(gene);
    nuclei->saveCoeffects(gene);
    nuclei->saveQuenching(gene);

    nuclei->updateSites(gene,tf);
    nuclei->updateSubgroups(gene);
    nuclei->updateCoeffects(gene);
    nuclei->updateQuenching(gene);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

void Organism::restoreSites(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring Sites" << endl;

  TF& tf = master_tfs->getTF(params[idx]->getTFName());


  int ngenes  = master_genes->size();

#ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreSites(gene, tf);
    //nuclei->updateSites(gene);
    nuclei->restoreAllOccupancy(gene);
    nuclei->restoreSubgroups(gene);
    nuclei->restoreCoeffects(gene);
    nuclei->restoreQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

/* if we move cooperativity distance we need to repopulate subgroups, but
quenching and coeffect interactions are unchanged */
void Organism::moveCoopD(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving cooperativity distance" << endl;


  int ngenes  = master_genes->size();

#ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    nuclei->saveSubgroups(gene);

    nuclei->updateSubgroups(gene);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

void Organism::restoreCoopD(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring cooperativity distance" << endl;


  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreAllOccupancy(gene);
    nuclei->restoreSubgroups(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

/* If we change lambda we need dont need to scan for new sites, but we do
need to update K and recalc occupancy and interations */
void Organism::moveLambda(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving lambda" << endl;

  TF& tf = master_tfs->getTF(params[idx]->getTFName());


  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    // dont bother saving K and Lambda since that takes nearly as long as calculating
    nuclei->updateKandLambda(gene, tf);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

void Organism::restoreLambda(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring lambda" << endl;

  TF& tf = master_tfs->getTF(params[idx]->getTFName());

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreAllOccupancy(gene);
    nuclei->updateKandLambda(gene, tf);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

// Moving Kmax is about the same as moving lambda. Maybe I should just ignore..
void Organism::moveKmax(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving kmax" << endl;

  TF& tf = master_tfs->getTF(params[idx]->getTFName());


  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    nuclei->updateK(gene, tf);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

void Organism::restoreKmax(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring kmax" << endl;

  TF& tf = master_tfs->getTF(params[idx]->getTFName());


  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreAllOccupancy(gene);
    nuclei->updateK(gene, tf);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

// if we move cooperativity we simply need to redo occupancy calculations
void Organism::moveKcoop(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving kcoop" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

void Organism::restoreKcoop(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring kcoop" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreAllOccupancy(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

/* If we change coeffect distances we need to repopulate the coeffects*/
void Organism::moveCoeffect(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving coeffects" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveModeOccupancy(gene);
    nuclei->saveQuenching(gene);
    nuclei->saveCoeffects(gene);

    nuclei->updateCoeffects(gene);
    nuclei->updateQuenching(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

void Organism::restoreCoeffect(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring coeffects" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreCoeffects(gene);
    nuclei->restoreQuenching(gene);
    nuclei->restoreModeOccupancy(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

/* If we allow a coefficient to switch from activator to repressor, we use this
function. If the bounds dont allow this, you can simply point the move generator
to move quenching or move activations */
void Organism::moveCoef(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving TF coefficient" << endl;

  double_param_ptr p = boost::any_cast<double_param_ptr>(params[idx]);
  double val  = p->getValue();
  double prev = p->getPrevious();

  if ( val >= 0 && prev >= 0)
  {
    int ngenes  = master_genes->size();

    /*
    #ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif
    */
    for (int j=0; j<ngenes; j++)
    {
      Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
      //nuclei->updateN();
      nuclei->updateR(gene);
    }
  }
  else if ( val <= 0 && prev <= 0)
    moveQuenchingCoef(idx);
  else
    moveQuenching(idx);
}

void Organism::restoreCoef(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring TF coefficient" << endl;

  double_param_ptr p = boost::any_cast<double_param_ptr>(params[idx]);
  double val  = p->getValue();
  double prev = p->getPrevious();

  if ( val >= 0 && prev >= 0)
  {
    int ngenes  = master_genes->size();

    /*
    #ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif
    */
    for (int j=0; j<ngenes; j++)
    {
      Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
      //nuclei->updateN();
      nuclei->updateR(gene);
    }
  }
  else if ( val <= 0 && prev <= 0)
    restoreQuenchingCoef(idx);
  else
    restoreQuenching(idx);
}

/* These functions are used if a quencher has been added or removed */
void Organism::moveQuenching(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving quenching" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveEffectiveOccupancy(gene);
    nuclei->saveQuenching(gene);
    nuclei->updateQuenching(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

void Organism::restoreQuenching(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring quenching" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreQuenching(gene);
    nuclei->restoreEffectiveOccupancy(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}


/* These functions are used if a the number of quenchers is unchanged */
void Organism::moveQuenchingCoef(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving quenching coefficient" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveEffectiveOccupancy(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

void Organism::restoreQuenchingCoef(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring quenching coefficient" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreEffectiveOccupancy(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

/* If we have moved coactivation or corepression efficiency */
void Organism::moveCoeffectEff(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving coeffect coefficient" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveModeOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

void Organism::restoreCoeffectEff(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring coeffect coefficient" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreModeOccupancy(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

/* if we have moved the promoter properties we simply call this */
void Organism::movePromoter(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving promoter parameter" << endl;

  int ngenes  = master_genes->size();

  /*
  #ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif
  */
  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->updateR(gene);
  }
}

/* this is somewhat awkward, but we dont need to do anything but rescore if
scale factors are used, so this is a null funtions */

void Organism::null_function(int idx)
{}


/*    Annealing   */

int Organism::getStateSize()
{
  int nparams = params.size();
  int size = 0;
  for (int i=0; i<nparams; i++)
    size += params[i]->getSize();
  
  return size;
}

void Organism::serialize(void * buf) const
{
  char * cbuf = (char *) buf;
  int nparams = params.size();
  for (int i = 0; i < nparams; ++i)
  {
    params[i]->serialize(cbuf);
    cbuf += params[i]->getSize();
  }   
}

void Organism::deserialize(void const *buf)
{
  char const * cbuf = (char const *) buf;
  int nparams = params.size();
  for (int i = 0; i < nparams; ++i)
  {
    params[i]->deserialize(cbuf);
    cbuf += params[i]->getSize();
  }
  ResetAll(0);
  score();
}



int    Organism::getDimension() const
{
  return params.size();
}

double Organism::get_score()
{
  return score_out;
}

void   Organism::generateMove(int idx, double theta)
{
  move_count++;
  
  previous_score_out = score_out;

  if (mode->getVerbose() >= 3)
    cerr << "generating move..." << endl;

  params[idx]->tweak(theta);

  if (!params[idx]->isOutOfBounds())
  {
    MFP fp = moves[idx];
    (this->*fp)(idx);
    score();
  }
  else
  {
    if (mode->getVerbose() >= 3)
      cerr << "move is out of bounds!" << endl;
    score_out = numeric_limits<double>::max();
    params[idx]->restore();
  }


}

void Organism::move(int idx)
{
  MFP fp = moves[idx];
  (this->*fp)(idx);
  score();
}

void   Organism::restoreMove(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "restoring move" << endl;

  params[idx]->restore();

  if (!params[idx]->isOutOfBounds())
  {
    MFP fp = restores[idx];
    (this->*fp)(idx);
  }
  score();
  if (mode->getVerbose() >= 3)
    cerr << "the score is now: " << setprecision(16) << score_out << endl;

  if (score_out != previous_score_out)
  {
    stringstream err;
    err << setprecision(16);
    err << "Restoring param with index "
    << idx
    << " failed. The Score was "
    << score_out
    << " but should have been "
    << previous_score_out;
    error(err.str());
  }

  if (mode->getVerbose() >= 3)
  {
    double tmp_score = score_out;
    ResetAll(idx);
    score();
    if (score_out != tmp_score)
    {
      stringstream err;
      err << setprecision(16);
      err << "Restoring param with index "
      << idx
      << " failed. The Score was "
      << tmp_score
      << " but should have been "
      << score_out;
      error(err.str());
    }
  }

}

void Organism::printParameters(ostream& os)
{
  os << "Parameters" << endl;
  int nparams = params.size();

  for (int i=0; i<nparams; i++)
  {
    if (i==0)
    {
      os << setw(5) << "idx";
      params[i]->printHeader(os);
    }
    os << setw(5) << i;
    params[i]->print(os);
  }
  os << endl;
}

string Organism::getParamName(int idx)
{
  stringstream out;

  if (params[idx]->is_tf_param())
    out << params[idx]->getTFName() << "_";
  out << params[idx]->getParamName();

  return out.str();

}

vector<double>& Organism::getN(int gidx)
{
  Gene& gene = master_genes->getGene(gidx);
  vector<double>& out = nuclei->getN(gene);
  return out;
}
