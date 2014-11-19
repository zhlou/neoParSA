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

void Score::sse()
{
  int ngenes = genes->size();
  
  for (int i=0; i<ngenes; i++)
  {
    scale_factor_ptr tscale = scale[i];
    double A = tscale->getA()->getValue();
    vector<double*>& gdata  = data[i];
    vector<double*>& gpred  = prediction[i];
    
    int length = gdata.size();
    double sse = 0;
    
    for (int j=0; j<length; j++)
    {
      double tdata = *gdata[j];
      double tpred = *gpred[j];

      tpred = tscale->unscale(tpred);
      
      double diff = tdata - tpred;
      sse += diff*diff;
    }
    sse *= A*A;
    scores[i] = sse;
  }
}

void Score::chisq()
{
  int ngenes = genes->size();
  
  for (int i=0; i<ngenes; i++)
  {
    scale_factor_ptr tscale = scale[i];
    double mw = tscale->unscale(min_weight);
    vector<double*>& gdata  = data[i];
    vector<double*>& gpred  = prediction[i];
    
    int length = gdata.size();
    double chisq = 0;
    
    for (int j=0; j<length; j++)
    {
      double tdata = *gdata[j];
      double tpred = *gpred[j];
      
      tpred = tscale->unscale(tpred);
      tdata = max(tdata, mw);
      tpred = max(tpred, mw);
      
      double diff = tdata - tpred;
      chisq += diff*diff/tdata;
    }
    scores[i] = chisq;
  }
}

void Score::pdiff()
{
  int ngenes = genes->size();
  
  for (int i=0; i<ngenes; i++)
  {
    scale_factor_ptr tscale = scale[i];
    double mw = tscale->unscale(min_weight);
    vector<double*>& gdata  = data[i];
    vector<double*>& gpred  = prediction[i];
    
    int length = gdata.size();
    double pdiff = 0;
    
    for (int j=0; j<length; j++)
    {
      double tdata = *gdata[j];
      double tpred = *gpred[j];

      tpred = tscale->unscale(tpred);
      tdata = max(tdata, mw);
      tpred = max(tpred, mw);
      
      double diff = tdata - tpred;
      pdiff += abs(diff)/tdata;
    }
    scores[i] = pdiff;
  }
}

void Score::rms()
{
  int ngenes = genes->size();
  
  for (int i=0; i<ngenes; i++)
  {
    scale_factor_ptr tscale = scale[i];
    double mw = tscale->unscale(min_weight);
    vector<double*>& gdata  = data[i];
    vector<double*>& gpred  = prediction[i];
    
    int length = gdata.size();
    double rms = 0;
    
    for (int j=0; j<length; j++)
    {
      double tdata = *gdata[j];
      double tpred = *gpred[j];

      tpred = tscale->unscale(tpred);
      tdata = max(tdata, mw);
      tpred = max(tpred, mw);
      
      double diff = tdata - tpred;
      rms += diff*diff;
    }
    rms /= length;
    rms = sqrt(rms);
    scores[i] = rms;
  }
}

// score and weight just as AhRam did in 2013 PLoS Genetics
void Score::arkim()
{
  int ngenes = genes->size();
  
  double max_area = 0;
  vector<double> area;
  area.resize(ngenes);
  
  for (int i=0; i<ngenes; i++)
  {
    scale_factor_ptr tscale = scale[i];
    vector<double*>& gdata  = data[i];
    vector<double*>& gpred  = prediction[i];
    
    int length  = gdata.size();
    double sse  = 0;
    area[i] = 0;
    for (int j=0; j<length; j++)
    {
      double tdata = *gdata[j];
      double tpred = *gpred[j];
      
      tdata = tscale->scale(tdata);
      
      area[i] += tdata;
      
      double diff = tdata - tpred;
      sse += diff*diff;
    }
    max_area = max(max_area, area[i]);
    scores[i] = sse;
  }
  for (int i=0; i<ngenes; i++)
    scores[i] *= max_area/area[i];
}

// score based on slopes
void Score::sum_slope_squares()
{
  int ngenes = genes->size();
  
  for (int i=0; i<ngenes; i++)
  {
    vector<double*>& gdata  = data[i];
    vector<double*>& gpred  = prediction[i];
    
    int length = gdata.size();
    double sss = 0;
    
    for (int j=1; j<length; j++)
    {
      double delta_data = (*gdata[j]) - (*gdata[j-1]);
      double delta_pred = (*gpred[j]) - (*gpred[j-1]);
      
      double diff = delta_data - delta_pred;
      sss += diff*diff;
    }
    scores[i] = sss;
  }
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
    //double tmin = scale[i]->unscale(min_weight);
    
    int ndata = data[i].size();
    for (int j=0; j<ndata; j++)
      gene_area[i] += max(scale[i]->scale(*data[i][j]), min_weight);
    max_area = max(max_area, gene_area[i]);
  }
  for (int i=0; i<ngenes; i++)
  {
    double A = scale[i]->getA()->getValue();
    weights[i] *= (max_area)/(gene_area[i]);
  }
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
    gene_height[i] = min_weight;
    for (int j=0; j<ndata; j++)
      gene_height[i] = max(gene_height[i],scale[i]->scale(*data[i][j]));
    max_height = max(max_height, gene_height[i]);
  }
  for (int i=0; i<ngenes; i++)
    weights[i] *= max_height/gene_height[i];
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
      gene_height[i] = max(gene_height[i], *data[i][j] );
    max_height = max(max_height, gene_height[i]);
  }
  for (int i=0; i<ngenes; i++)
  {
    weights[i] *= 1/(gene_height[i]*gene_height[i]);
    //cerr << "weight for " << genes->getGene(i).getName() << ": " << weights[i] << endl;
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
    scoreFunc = &Score::sse;
  else if (mode->scoreFunction() == "chisq")
    scoreFunc = &Score::chisq;
  else if (mode->scoreFunction() == "pdiff")
    scoreFunc = &Score::pdiff;
  else if (mode->scoreFunction() == "rms")
    scoreFunc = &Score::rms;
  else if (mode->scoreFunction() == "arkim")
    scoreFunc = &Score::arkim;
  else if (mode->scoreFunction() == "sss")
    scoreFunc = &Score::sum_slope_squares;
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
 

void Score::checkScale(ostream& os)
{
  
  int ngenes = genes->size();
  
  // store old data
  vector<vector <double> > saved_data;
  vector<vector <double> > saved_pred;
  saved_data.resize(ngenes);
  saved_pred.resize(ngenes);
  for (int i=0; i<ngenes; i++)
  {
    int length = data[i].size();
    for (int j=0; j<length; j++)
    {
      saved_data[i].push_back(*data[i][j]);
      saved_pred[i].push_back(*prediction[i][j]);
    }
  }
  
  // test mutiplying data and prediction by 1 through 5
  for (int i=1; i<6; i++)
  {
    for (int j=0; j<ngenes; j++)
    {
      scale[j]->getA()->set(i);
      int length = data[j].size();
      for (int k=0; k<length; k++)
        *data[j][k] = scale[j]->unscale(saved_data[j][k]);
    }
    setWeights();
    os << "when scale factor is " << i << " score is " << getScore() << endl;
  }
}
    
  
double Score::getScore()
{
  int ngenes = genes->size();
  score = 0;
  
  (this->*scoreFunc)();
  
  for (int i=0; i<ngenes; i++)
  {   
    double tmpscore = scores[i] * weights[i];
    score += tmpscore;
  }
  return score;
}
    
void Score::print(ostream& os)
{
  int ngenes = genes->size();
  
  os << setw(14) << "Name";
  os << setw(14) << "Score";
  os << setw(14) << "Weight";
  os << setw(14) << "Product" << endl;
  
  double total = 0.0;
  
  for (int i=0; i<ngenes; i++)
  {
    total += scores[i];
    os << setw(14) << genes->getGene(i).getName();
    os << setw(14) << setprecision(5) << scores[i];
    os << setw(14) << setprecision(5) << weights[i];
    os << setw(14) << setprecision(5) << scores[i]*weights[i] << endl;
  }
  double weighted_score = getScore();
  os << "----------------------------------------------------------" << endl;
  os << setw(14) << "Total";
  os << setw(14) << setprecision(5) << total;
  os << setw(14) << " ";
  os << setw(14) << setprecision(5) << weighted_score << endl;
}
    
   


