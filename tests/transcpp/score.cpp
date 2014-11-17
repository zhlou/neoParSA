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

double sse(vector<double*>& x, vector<double*>& y, scale_factor_ptr s, int min_weight)
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

double chi_sq(vector<double*>& x, vector<double*>& y, scale_factor_ptr s, int min_weight)
{
  double tmin = s->unscale(min_weight);
  int length = x.size();
  double chisq = 0.0;
  
  for (unsigned int i=0; i<length; i++)
  {
    double tx = *x[i];
    double ty = *y[i];
    ty = s->unscale(ty);
    
    tx <- max(tx, tmin);
    ty <- max(ty, tmin);
    
    double diff = tx-ty;
    chisq += (diff*diff)/tx;
  }
  return chisq/length;
}

double pdiff(vector<double*>& x, vector<double*>& y, scale_factor_ptr s,  int min_weight)
{
  double tmin = s->unscale(min_weight);
  int length = x.size();
  double pdiff=0;
  
  for (unsigned int i=0; i<length; i++)
  {
    double tx = *x[i];
    double ty = *y[i];
    ty = s->unscale(ty);
    
    tx <- max(tx, tmin);
    ty <- max(ty, tmin);
    
    double diff = tx-ty;
    pdiff += abs(diff)/tx;
  }
  return pdiff/length;
}


/*****************  Weight Functions  *******************************************/

void Score::area()
{
  double max_area = 0.0;
  int    ngenes   = genes->size();
  vector<double> gene_area;
  gene_area.resize(ngenes,0.0);
  
  for (int i=0; i<ngenes; i++)
  {
    double tmin = scale[i]->unscale(min_weight);
    int ndata = data[i].size();
    for (int j=0; j<ndata; j++)
      gene_area[i] += max(*data[i][j], tmin);
    max_area = max(max_area, gene_area[i]);
  }
  for (int i=0; i<ngenes; i++)
    weights[i] *= max_area/gene_area[i];
}

void Score::height()
{
  double max_height = 0.0;
  int    ngenes   = genes->size();
  vector<double> gene_height;
  gene_height.resize(ngenes,0.0);
  
  
  for (int i=0; i<ngenes; i++)
  {
    int ndata = data[i].size();
    gene_height[i] = scale[i]->unscale(min_weight);
    for (int j=0; j<ndata; j++)
      gene_height[i] = max(gene_height[i],*data[i][j]);
    max_height = max(max_height, gene_height[i]);
  }
  for (int i=0; i<ngenes; i++)
    weights[i] *= 1/gene_height[i];
}

void Score::height2()
{
  double max_height = 0.0;
  int    ngenes   = genes->size();
  vector<double> gene_height;
  gene_height.resize(ngenes,0.0);
  
  
  for (int i=0; i<ngenes; i++)
  {
    int ndata = data[i].size();
    gene_height[i] = scale[i]->unscale(min_weight);
    for (int j=0; j<ndata; j++)
      gene_height[i] = max(gene_height[i],(*data[i][j]));
    max_height = max(max_height, gene_height[i]);
  }
  for (int i=0; i<ngenes; i++)
  {
    weights[i] *= 1/(gene_height[i]*gene_height[i]);
  }
}

void Score::gene_number()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
    weights[i] /= ngenes;
}

void Score::nuc_number()
{
  int ngenes = genes->size();
  int nnuc   = data[0].size();
  for (int i=0; i<ngenes; i++)
    weights[i] /= nnuc;
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
    scoreFunc = bind(&sse, _1, _2, _3, _4);
  else
  {
    cerr << "ERROR: could not find score function with name " << mode->scoreFunction() << endl;
    exit(1);
  }
  
  min_weight = mode->getMinWeight();
  if (mode->getWeightFunction("area"))
    weight_functs.push_back(&Score::area);
  if (mode->getWeightFunction("height"))
    weight_functs.push_back(&Score::height);
  if (mode->getWeightFunction("height2"))
    weight_functs.push_back(&Score::height2);
  if (mode->getWeightFunction("ngenes"))
    weight_functs.push_back(&Score::gene_number);
  if (mode->getWeightFunction("nnuc"))
    weight_functs.push_back(&Score::nuc_number);

  
  setWeights();
}

/* the sum of all weights should always equal 1, that way the optimizer
does not attempt to lower the weights instead of fit the data! */
void Score::setWeights()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
    weights[i] = 1.0;
  
  int nfuncs = weight_functs.size();
  for (int i=0; i<nfuncs; i++)
  {
    WFP fp = weight_functs[i];
    (this->*fp)();
  }
  
  
}
 
  
double Score::getScore()
{
  int ngenes = genes->size();
  score = 0;
  
  for (int i=0; i<ngenes; i++)
  {   
    double tmpscore = scoreFunc(data[i],prediction[i],scale[i],min_weight) * weights[i];
    scores[i] = tmpscore;
    score += tmpscore;
  }
  return score;
}
    
void Score::print(ostream& os)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    os << setw(14) << genes->getGene(i).getName();
    os << setw(14) << setprecision(5) << scores[i] << endl;
  }
  os << "-----------------------------" << endl;
  os << setw(14) << "Total";
  os << setw(14) << setprecision(5) << getScore() << endl;
}
    
   


