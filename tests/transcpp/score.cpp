/*********************************************************************************
*                                                                                *
*     score.cpp                                                                  *
*                                                                                *
*     Presumably we could use any number of functions to score the model.        *
*     We have some set of predictions and outputs, expressed here as pointer     *
*     vectors, as well as some rule for weighting the scores of each             *
*                                                                                *
*********************************************************************************/


#include "score.h"


/*****************  Scoring Functions  ******************************************/

double sse(vector<double*>& x, vector<double*>& y, scale_factor_ptr s)
{
  int length = x.size();
  double sse=0;
  for (unsigned int i=0; i<length; i++)
  { 
    double tx = *x[i];
    double ty = *y[i];
    ty = s->unscale(ty);
    
    double diff = tx-ty;
    sse += diff*diff;
  }
  return sse;
}

/*****************  Weight Functions  *******************************************/

double area(vector<double*>& x)
{
  int n = x.size();
  double sum = 0;
  for (int i=0; i<n; i++)
    sum += *x[i];
  
  return 1/sum;
}

double none(vector<double*>& x)
{
  return 1;
}


/********************************************************************************/



/*    Constructors    */

Score::Score() {};

Score::Score(Organism* parent) { set(parent); }


/*    Setters     */

void Score::set(Organism* parent)
{
  mode  = parent->getMode();
  genes = parent->getGenes();
  ids   = parent->getIDs();
  
  conc_ptr rates = parent->getRateData();
    
  int ngenes = genes->size();
  int nids   = ids.size();
  
  data.resize(ngenes);
  prediction.resize(ngenes);
  scores.resize(ngenes);
  weights.resize(ngenes);
  scale.resize(ngenes);
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    const string& gname = gene.getName();
    scale[i] = rates->getScale(gname); 
    data[i].resize(nids);
    prediction[i].resize(nids);
    
    for (int j=0; j<nids; j++)
    {
      data[i][j]       = parent->getData(gene, ids[j], false);
      prediction[i][j] = parent->getPrediction(gene, ids[j]);
    }
  }
  
  if (mode->scoreFunction() == "sse")
    scoreFunc = bind(&sse, _1, _2, _3);
  else
  {
    cerr << "ERROR: could not find score function with name " << mode->scoreFunction() << endl;
    exit(1);
  }
  
  if (mode->weightFunction() == "area")
    weightFunc = bind(&area, _1);
  else if (mode->weightFunction() == "none")
    weightFunc = bind(&none, _1);
  else
  {
    cerr << "ERROR: could not find score function with name " << mode->weightFunction() << endl;
    exit(1);
  }
  
  setWeights();
}

/* the sum of all weights should always equal 1, that way the optimizer
does not attempt to lower the weights instead of fit the data! */
void Score::setWeights()
{
  int ngenes = genes->size();
  double sum = 0;
  for (int i=0; i<ngenes; i++)
  {
    double tmp = weightFunc(data[i]);
    weights[i] = tmp;
    sum += tmp;
  }
  for (int i=0; i<ngenes; i++)
    weights[i] /= sum;
}
  
double Score::getScore()
{
  int ngenes = genes->size();
  score = 0;
  
  for (int i=0; i<ngenes; i++)
  {   
    double tmpscore = scoreFunc(data[i],prediction[i],scale[i]) * weights[i];
    scores[i] = tmpscore;
    score += tmpscore;
  }
  return score;
}
    
    
   


