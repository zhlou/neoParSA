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
#include <boost/range/adaptor/map.hpp>
#include <limits>

#include <boost/random/linear_congruential.hpp>      
#include <boost/random/uniform_int.hpp>             
#include <boost/random/uniform_real.hpp>            
#include <boost/random/variate_generator.hpp>

using namespace boost::adaptors;

# define foreach_ BOOST_FOREACH


/*    Constructors   */

Organism::Organism() :
    mode(mode_ptr(new Mode)),
    distances(distances_ptr(new DistanceContainer)),
    master_tfs(tfs_ptr(new TFContainer)),
    master_genes(genes_ptr(new GeneContainer)),
    tfdata(conc_ptr(new ConcContainer)),
    ratedata(conc_ptr(new ConcContainer)),
    promoters(promoters_ptr(new PromoterContainer)),
    scale_factors(scale_factors_ptr(new ScaleFactorContainer)),
    coeffects(coeffects_ptr(new CoeffectContainer)),
    score_class(score_ptr(new Score)),
    coops(coops_ptr(new CooperativityContainer))
{}

Organism::Organism(ptree& pt) :
    mode(mode_ptr(new Mode)),
    distances(distances_ptr(new DistanceContainer)),
    master_tfs(tfs_ptr(new TFContainer)),
    master_genes(genes_ptr(new GeneContainer)),
    tfdata(conc_ptr(new ConcContainer)),
    ratedata(conc_ptr(new ConcContainer)),
    promoters(promoters_ptr(new PromoterContainer)),
    scale_factors(scale_factors_ptr(new ScaleFactorContainer)),
    coeffects(coeffects_ptr(new CoeffectContainer)),
    score_class(score_ptr(new Score)),
    coops(coops_ptr(new CooperativityContainer))
{
  int i = 0;
  mode->read(pt);
  readAnnealing(pt);
  distances->add(pt);
  promoters->read(pt);
  master_tfs->add(pt);
  master_genes->setPromoters(promoters);
  master_genes->read(pt);
  scale_factors->read(pt);
  ratedata->read(pt, scale_factors, "RateData");
  tfdata->read(pt, scale_factors, "TFData");
  coops->read(pt, distances);
  coeffects->read(pt,distances);
  master_tfs->setCoops(coops);
  master_tfs->setCoeffects(coeffects);
  
  master_tfs->getParameters(params);
  promoters->getParameters(params);
  scale_factors->getParameters(params);
  coops->getParameters(params);
  coeffects->getParameters(params);
  
  populate_nuclei();
  score_class->set(this);
  
  score();
  
  setMoves();
  checkParameters(); // make sure we have a move function for each param
  
}


/*    Setters   */

 
/*   For annealing we only want to do computations on on nuclei that we have
data for. Additionally, computations on vectors are much more efficient than
those on maps. In light of these two ideas, we first must read the rate data,
then we will populate nuclei for each data point we have. Missing data is 
given the value of NA = numeric_limits<double>::infinity */

void Organism::populate_nuclei()
{
  if (mode->multipleSubgroups())
    populate_nuclei_multiple();
  else
    populate_nuclei_single();
}

/* Makes only one nuclei object, with all subgroups and quenching interactions
in it. This is faster when sites are dense and there are many interactions */
void Organism::populate_nuclei_multiple()
{
  int nnuc   = ratedata->size(); // this is the number of possible datapoints to score
  int ngenes = master_genes->size();
  for (int i = 0; i<nnuc; i++)
  {
    int id = ratedata->index2id(i);
    ids.push_back(id);
    
    //nuclei_ptr nuclei(new Nuclei(this));
    tfs_ptr    nuc_tfs(new TFContainer);     // each nuc may have different tfs

    // what tf data do we have here?
    int ntfs = master_tfs->size();
    for (int j=0; j<ntfs; j++)
    {
      tf_ptr tf = master_tfs->getTFptr(j);
      const string& tfname = tf->getName();
      if (tfdata->getConcByID(id, tfname, false) != 0)
        nuc_tfs->add(tf);
    }
   
    bool added   = false;
    int  nnuclei = nuclei.size();
    for (int j=0; j<nnuclei; j++)
    {
      if (nuclei[j]->compareTFs(nuc_tfs))
      {
        nuclei[j]->addNuc(id);
        added=true;
      }
    }
    
    if (added == false)
    {
      nuclei_ptr tmp_nuclei(new Nuclei(this, nuc_tfs));
      tmp_nuclei->addNuc(id);
      nuclei.push_back(tmp_nuclei);
    }
  }

  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
    nuclei[i]->create();
  
  if (nuclei.size() == 0)
  {
    cerr << "ERROR: populate_nuclei() did not find any nuclei!" << endl;
    exit(1);
  }
}

void Organism::populate_nuclei_single()
{
  
  nuclei_ptr all_nuclei(new Nuclei(this, master_tfs));
  nuclei.push_back(all_nuclei);
  
  int ngenes = master_genes->size();
  int nnuc   = ratedata->size(); // this is the number of possible datapoints to score
  for (int i = 0; i<nnuc; i++)
  {
    int id = ratedata->index2id(i);
    ids.push_back(id);
    all_nuclei->addNuc(id);
  }
      
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
    nuclei[i]->create();
  
  if (nnuc == 0)
  {
    cerr << "ERROR: populate_nuclei() did not find any nuclei!" << endl;
    exit(1);
  }
}

double* Organism::getPrediction(Gene& gene, int id)
{
  int nnuclei = nuclei.size();
  
  for (int i=0; i<nnuclei; i++)
  {
    vector<int>& tmp_ids = nuclei[i]->getIDs();
    for (int j=0; j<tmp_ids.size(); j++)
    {
      if (tmp_ids[j] == id)
        return &(nuclei[i]->getRate(gene, id));
    }
  }
}

double* Organism::getData(Gene& gene, int id, bool scaled)
{
  const string& gname = gene.getName();

  return &(ratedata->getConcByID(gname, id, scaled));
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
  mode->write(output);
  output.add_child("annealer_input",annealer_input);
  output.add_child("move",move);
  output.add_child("count_criterion",count_criterion);
  output.add_child("mix",mix);
  output.add_child("lam",lam);
  distances->write(output);
  promoters->write(output);
  master_tfs->write(output);
  coops->write(output);
  coeffects->write(output);
  master_genes->write(output);
  scale_factors->write(output);
  ratedata->write(output, "RateData");
  tfdata->write(output, "TFData");
}

void Organism::readAnnealing(ptree& pt)
{
  annealer_input  = pt.get_child("annealer_input");
  move            = pt.get_child("move");
  count_criterion = pt.get_child("count_criterion");
  mix             = pt.get_child("mix");
  lam             = pt.get_child("lam");
}
  

void Organism::printRate(ostream& os)
{
  int w = 12;
  int ngenes = master_genes->size();
  os << setw(w) << setprecision(5) << "id";
  for (int i=0; i<ngenes; i++)
    os << setw(w) << master_genes->getGene(i).getName();
  os << endl;
  
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    vector<int>& IDs = nuc->getIDs();
    int nids         = IDs.size();
    for (int j=0; j<nids; j++)
    {
      os << setw(w) << IDs[j];
      for (int k=0; k<ngenes; k++)
      {
        Gene& gene = master_genes->getGene(k);
        vector<double>& rate = nuc->getRate(gene);
        os << setw(w) << rate[j];
      }
      os << endl;
    }
  }
}
      
void Organism::printRateData(ostream& os)
{
  int w = 12;
  int ngenes = master_genes->size();
  os << setw(w) << setprecision(5) << "id";
  for (int i=0; i<ngenes; i++)
    os << setw(w) << master_genes->getGene(i).getName();
  os << endl;
  
  int nnuc = ratedata->size();
  for (int i=0; i<nnuc; i++)
  {
    int id = ratedata->index2id(i);
    os << setw(w) << id;
    for (int i=0; i<ngenes; i++)
      os << setw(w) << ratedata->getConcByID(master_genes->getGene(i).getName(), id, true);
    os << endl;
  }
}

void Organism::printSites(ostream& os)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei[i]->printSites(os);
  }
}

void Organism::printSites(Gene& gene, ostream& os)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei[i]->printSites(gene, os);
  }
}

void Organism::printSites(TF& tf, ostream& os)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei[i]->printSites(tf, os);
  }
}

void Organism::printSites(Gene& gene, TF& tf, ostream& os)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei[i]->printSites(gene, tf, os);
  }
}

void Organism::printSubgroups(Gene& gene, ostream& os)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei[i]->printSubgroups(gene, os);
  }
}

void Organism::printOccupancy(Gene& gene, ostream& os)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei[i]->printOccupancy(gene, os);
  }
}

void Organism::printModeOccupancy(Gene& gene, ostream& os)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei[i]->printModeOccupancy(gene, os);
  }
}

void Organism::printEffectiveOccupancy(Gene& gene, ostream& os)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei[i]->printEffectiveOccupancy(gene, os);
  }
}

/*    Move    */

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

void Organism::setMoves()
{
  move_map["threshold"]      = &Organism::moveThreshold;
  move_map["kmax"]           = &Organism::moveKmax;
  move_map["coef"]           = &Organism::moveCoef;
  move_map["Q"]              = &Organism::moveQ;
  move_map["Rmax"]           = &Organism::moveQ;
  move_map["Theta"]          = &Organism::moveQ;
  move_map["lambda"]         = &Organism::moveLambda;
  move_map["ScaleFactor"]    = &Organism::moveScaleFactor;
  move_map["Kcoop"]          = &Organism::moveKcoop;
  move_map["CoeffectEff"]    = &Organism::moveCoeffect;
  
  restore_map["threshold"]      = &Organism::restoreThreshold;
  restore_map["kmax"]           = &Organism::restoreKmax;
  restore_map["coef"]           = &Organism::restoreCoef;
  restore_map["Q"]              = &Organism::moveQ;
  restore_map["Rmax"]           = &Organism::moveQ;
  restore_map["Theta"]          = &Organism::moveQ;
  restore_map["lambda"]         = &Organism::restoreLambda;
  restore_map["ScaleFactor"]    = &Organism::moveScaleFactor;
  restore_map["Kcoop"]          = &Organism::restoreKcoop;
  restore_map["CoeffectEff"]    = &Organism::restoreCoeffect;
}
  
void Organism::checkParameters()
{
  int nparams = params.size();
  for (int i=0; i<nparams; i++)
  {
    const string& param_name = params[i]->getParamName();
    if (move_map.find(param_name) == move_map.end())
    {
      cerr << "ERROR: Could not find move function for param " << param_name << endl;
      exit(1);
    }
    if (restore_map.find(param_name) == restore_map.end())
    {
      cerr << "ERROR: Could not find restore function for param " << param_name << endl;
      exit(1);
    }
  }
}

void Organism::calc_f()
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei[i]->calc_f();
  }
}

void Organism::moveScaleFactor(int idx)
{
  //ratedata->scale_data();
  score_class->setWeights();
}

void Organism::resetAll()
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->updateScores();
    nuc->updateSites();
    nuc->updateSubgroups();
    nuc->updateQuenching();
    nuc->calc_f();
    nuc->calcQuenching();
    nuc->updateN();
    nuc->updateR();
  }
}

void Organism::moveThreshold(int idx)
{
  TF& tf = master_tfs->getTF(params[idx]->getTFName());
  
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->saveAllOccupancy();
    nuc->saveSites(tf);         
    nuc->saveScores(tf);
    nuc->saveSubgroups();
    nuc->saveCoeffects();
    nuc->saveQuenching();
    
    nuc->updateScores(tf);
    nuc->updateSites(tf);
    
    nuc->updateSubgroups();
    nuc->updateCoeffects();
    nuc->updateQuenching();
    nuc->calc_f();
    nuc->calcCoeffects();
    nuc->calcQuenching();
    nuc->updateN();
    nuc->updateR();
  }
}

void Organism::restoreThreshold(int idx)
{
  TF& tf = master_tfs->getTF(params[idx]->getTFName());
  
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->restoreScores(tf);
    nuc->restoreSites(tf);
    nuc->restoreAllOccupancy();
    nuc->restoreSubgroups();
    nuc->restoreCoeffects();
    nuc->restoreQuenching();
    nuc->updateN();
    nuc->updateR();
  }
  score();
}

void Organism::moveKmax(int idx)
{
  TF& tf = master_tfs->getTF(params[idx]->getTFName());
  
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->saveAllOccupancy();
    nuc->updateK(tf);
    nuc->calc_f();
    nuc->calcCoeffects();
    nuc->calcQuenching();
    nuc->updateN();
    nuc->updateR();
  }

}


void Organism::restoreKmax(int idx)
{
  TF& tf = master_tfs->getTF(params[idx]->getTFName());

  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->restoreAllOccupancy();
    nuc->updateK(tf);
    nuc->updateN();
    nuc->updateR();
  }
}

void Organism::moveKcoop(int idx)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->saveAllOccupancy();
    nuc->calc_f();
    nuc->calcCoeffects();
    nuc->calcQuenching();
    nuc->updateN();
    nuc->updateR();
  }
}

void Organism::restoreKcoop(int idx)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->restoreAllOccupancy();
    nuc->updateN();
    nuc->updateR();
  }
}


void Organism::moveLambda(int idx)
{
  
  TF& tf = master_tfs->getTF(params[idx]->getTFName());

  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->saveAllOccupancy();
    nuc->updateKandLambda(tf);
    nuc->calc_f();
    nuc->calcCoeffects();
    nuc->calcQuenching();
    nuc->updateN();
    nuc->updateR();
  }
  
}


void Organism::restoreLambda(int idx)
{
  
  TF& tf = master_tfs->getTF(params[idx]->getTFName());

  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->restoreAllOccupancy();
    nuc->updateKandLambda(tf);
    nuc->updateN();
    nuc->updateR();
  }

}

void Organism::moveCoef(int idx)
{
  double val  = params[idx]->getValue();
  double prev = params[idx]->getPrevious();
  
  if ( val >= 0 && prev >= 0)
  {
    int nnuclei = nuclei.size();
    for (int i=0; i<nnuclei; i++)
    {
      nuclei_ptr nuc = nuclei[i];
      nuc->updateN();
      nuc->updateR();
    }
  } 
  else if ( val <= 0 && prev <= 0)
    moveQuenchingCoef(idx);
  else
    moveQuenching(idx);
}

void Organism::restoreCoef(int idx)
{
  double val  = params[idx]->getValue();
  double prev = params[idx]->getPrevious();
  
  if ( val >= 0 && prev >= 0)
  {
    int nnuclei = nuclei.size();
    for (int i=0; i<nnuclei; i++)
    {
      nuclei_ptr nuc = nuclei[i];
      nuc->updateN();
      nuc->updateR();
    }
  }
  else if ( val <= 0 && prev <= 0)
    restoreQuenchingCoef(idx);
  else
    restoreQuenching(idx);
}

void Organism::moveQuenchingCoef(int idx)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->saveEffectiveOccupancy();
    nuc->calcQuenching();
    nuc->updateN();
    nuc->updateR();
  }
}

void Organism::restoreQuenchingCoef(int idx)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->restoreEffectiveOccupancy();
    nuc->updateN();
    nuc->updateR();
  }
}

void Organism::moveCoeffect(int idx)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->saveModeOccupancy();
    nuc->calcCoeffects();
    nuc->calcQuenching();
    nuc->updateN();
    nuc->updateR();
  }
}

void Organism::restoreCoeffect(int idx)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->restoreModeOccupancy();
    nuc->updateN();
    nuc->updateR();
  }
}

void Organism::moveQuenching(int idx)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->saveEffectiveOccupancy();
    nuc->saveQuenching();
    nuc->updateQuenching();
    nuc->calcQuenching();
    nuc->updateN();
    nuc->updateR();
  }
}

void Organism::restoreQuenching(int idx)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->restoreQuenching();
    nuc->restoreEffectiveOccupancy();
    nuc->updateN();
    nuc->updateR();
  }
}

void Organism::moveQ(int idx)
{
  int nnuclei = nuclei.size();
  for (int i=0; i<nnuclei; i++)
  {
    nuclei_ptr nuc = nuclei[i];
    nuc->updateR();
  }
}

/*    Annealing   */

void Organism::serialize(void *buf) const
{
  int nparams = params.size();
  double *dest = static_cast<double *>(buf); // new style cast
  for (int i = 0; i < nparams; ++i)
      dest[i] = params[i]->getValue();
}

void Organism::deserialize(void const *buf)
{
  int nparams = params.size();
  double const *from = static_cast<double const *>(buf);
  for (int i = 0; i < nparams; ++i) 
      params[i]->set(from[i]);
  resetAll();
  score();
}

int    Organism::getStateSize() 
{
  return params.size() * sizeof(double);
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
  params[idx]->tweak(theta);

  if (!params[idx]->isOutOfBounds())
  {
    const string& param_name = params[idx]->getParamName();
    MFP fp = move_map[param_name];
    (this->*fp)(idx);
    score();
  }
  else
  {
    score_out = numeric_limits<double>::max();
    params[idx]->restore();
  }
}


void   Organism::restoreMove(int idx)
{
  params[idx]->restore();

  if (!params[idx]->isOutOfBounds())
  {
    const string& param_name = params[idx]->getParamName();
    MFP fp = restore_map[param_name];
    (this->*fp)(idx);
  }
  score();
}

void Organism::printParameters(ostream& os)
{
  int nparams = params.size();
  for (int i=0; i<nparams; i++)
  {
    os << setw(5) << i;
    params[i]->print(os);
  }
}
