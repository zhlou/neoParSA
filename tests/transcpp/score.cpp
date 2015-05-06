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
#include "utils.h"


/*****************  Scoring Functions  ******************************************/

/* Be careful when creating scoring functions. These must NEVER SCALE THE DATA!!!
If the data is scaled during anealing, most score functions will prefer the
lowest scale factor possible. This is the case of the score function used in 
Kim2013, which is now provided only for reverse-compatibility */

/* A scale factor says how the data would have to be scaled to fit the rates
predicted. Unscaling means changing the rates to fit the data, which turns
out to be a more useful function */

void Score::sse()
{
  int ngenes = genes->size();
  
  for (int i=0; i<ngenes; i++)
  {
    //cerr << genes->getGene(i).getName() << endl;
    
    scale_factor_ptr tscale = scale[i];
    //double A = tscale->getA()->getValue();
    //vector<double*>& gdata  = data[i];
    
    vector<double>& wgdata  = weighted_data[i];
    vector<double*>& gpred  = prediction[i];
    
    int length = wgdata.size();
    double tweight = weights[i];
    //cerr << tweight << endl;
    double sse = 0;
    
    for (int j=0; j<length; j++)
    {
      //double tdata = *gdata[j];
      
      double tdata = wgdata[j];
      double tpred = *gpred[j];
      
      //cerr << *data[i][j] << " ";
      //cerr << tdata << " ";
      //cerr << tpred << " ";
      
      //tpred = tpred

      tpred = tscale->unscale(tpred * tweight);
      
      //cerr << tpred << endl;
       
      //cerr << tdata << " " << tpred << endl;
      
      double diff = tdata - tpred;
      sse += diff*diff;
    }
    scores[i] = sse;
  }
}

/* a real chisq score! Note that this is unstable when the data approaches 0,
so instead we substitute the minimum weight specified in the mode node for values
lower than minimum weight */

void Score::chisq()
{
  int ngenes = genes->size();
  
  for (int i=0; i<ngenes; i++)
  {
    //cerr << genes->getGene(i).getName() << endl;
    
    scale_factor_ptr tscale = scale[i];
    vector<double*>& gdata  = data[i];
    vector<double*>& gpred  = prediction[i];
    
    int length = gdata.size();
    double chisq = 0;
    
    for (int j=0; j<length; j++)
    {
      double tdata = *gdata[j];
      double tpred = *gpred[j];
      
      tpred = tscale->unscale(tpred);
      
      double diff = tdata - tpred;
      //cerr << "( " << tdata << " - " << tpred << " )^2 / " << max(tdata, min_data) << " = " << diff*diff/max(tdata, min_data) << endl;
      chisq += diff*diff/max(tdata, min_data);
    }
    scores[i] = chisq;
  }
}

/* The percent difference between data and fit. Same as chisq, but with abs() 
instead of square */

void Score::pdiff()
{
  int ngenes = genes->size();
  
  for (int i=0; i<ngenes; i++)
  {
    scale_factor_ptr tscale = scale[i];
    vector<double*>& gdata  = data[i];
    vector<double*>& gpred  = prediction[i];
    
    int length = gdata.size();
    double pdiff = 0;
    
    for (int j=0; j<length; j++)
    {
      double tdata = *gdata[j];
      double tpred = *gpred[j];

      tpred = tscale->unscale(tpred);
      
      double diff = tdata - tpred;
      pdiff += abs(diff)/max(tdata, min_data);
    }
    scores[i] = pdiff;
  }
}

/* root mean squared differences */

void Score::rms()
{
  int ngenes = genes->size();
  
  for (int i=0; i<ngenes; i++)
  {
    scale_factor_ptr tscale = scale[i];
    vector<double*>& gdata  = data[i];
    vector<double*>& gpred  = prediction[i];
    
    int length = gdata.size();
    double rms = 0;
    
    for (int j=0; j<length; j++)
    {
      double tdata = *gdata[j];
      double tpred = *gpred[j];

      tpred = tscale->unscale(tpred);
      
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

/* Takes the "slopes" of each pair of data points and finds the sum
of squared slope differences between data and model. Not sure how well
this will work for the boundary conditions */

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

/* We run into a problem when different constructs have different magnitudes. 
Higher expressing constructs (or higher bit-content data) result in higher
scores with many types of scoring functions. The following functions attempt to
set weights for different genes in the scoring function, each according to 
different rules. */

/* Multiple weight functions can be used. This might be useful to compare scores
for fits with multiple constructs and different numbers of nuclei. If we divide
scores by the number of genes and nuclei, scores will be roughly comparable. */


/* score all constructs as if the area under the curve were equal to max_weight */
void Score::area()
{
  int    ngenes   = genes->size();
  vector<double> gene_area;
  gene_area.resize(ngenes,0.0);
  
  for (int i=0; i<ngenes; i++)
  {  
    int ndata = data[i].size();
    for (int j=0; j<ndata; j++)
      gene_area[i] += max(*data[i][j], min_data);
  }
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    weights[i] = gene.getWeight();
    weights[i] *= scale_to / gene_area[i];
    int ndata = data[i].size();
    for (int j=0; j<ndata; j++)
      weighted_data[i][j] = *data[i][j] * weights[i];
  }
}

/* score all constructs as if the maximum height were equal to max_weight */
void Score::height()
{
  int    ngenes   = genes->size();
  vector<double> gene_height;
  gene_height.resize(ngenes,0.0);

  for (int i=0; i<ngenes; i++)
  {
    int ndata = data[i].size();
    gene_height[i] = min_data;
    for (int j=0; j<ndata; j++)
      gene_height[i] = max(gene_height[i],*data[i][j]);
  }
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    weights[i] = gene.getWeight();
    weights[i] *= scale_to / gene_height[i];
    int ndata  = data[i].size();
    for (int j=0; j<ndata; j++)
      weighted_data[i][j] = *data[i][j] * weights[i];
  }
}

/* if we arent scaling the data we still need to populate weights and 
weighted data for scoring */

void Score::none()
{
  int    ngenes   = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    weights[i] = gene.getWeight();
    int ndata  = data[i].size();
    for (int j=0; j<ndata; j++)
      weighted_data[i][j] = *data[i][j] * weights[i];
  }
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
  
  table_ptr rates = parent->getRateData();
    
  int ngenes = genes->size();
  int nids   = ids.size();
  
  data.resize(ngenes);
  weighted_data.resize(ngenes);
  prediction.resize(ngenes);
  scores.resize(ngenes);
  weights.resize(ngenes);
  scale.resize(ngenes);
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    if (!gene.getInclude()) continue;
    scale[i] = gene.getScale(); 
    data[i].resize(nids);
    weighted_data[i].resize(nids);
    prediction[i].resize(nids);
    
    for (int j=0; j<nids; j++)
    {
      data[i][j]       = parent->getData(gene, ids[j]);
      prediction[i][j] = parent->getPrediction(gene, ids[j]);
    }
  }
  
  if (mode->getScoreFunction() == "sse")
    scoreFunc = &Score::sse;
  else if (mode->getScoreFunction() == "chisq")
    scoreFunc = &Score::chisq;
  else if (mode->getScoreFunction() == "pdiff")
    scoreFunc = &Score::pdiff;
  else if (mode->getScoreFunction() == "rms")
    scoreFunc = &Score::rms;
  else if (mode->getScoreFunction() == "arkim")
    scoreFunc = &Score::arkim;
  else if (mode->getScoreFunction() == "sss")
    scoreFunc = &Score::sum_slope_squares;
  else
  {
    stringstream err;
    err << "ERROR: could not find score function with name " << mode->getScoreFunction() << endl;
    error(err.str());
  }
  
  min_data = mode->getMinData();
  
  if (mode->getScaleData())
  {
    string type = mode->getScaleDataType();
    scale_to    = mode->getScaleTo();
    if (type == string("area"))
      area();
    else if (type == string("height"))
      height();
    else
      none();
  }
  else
    none();
  
  divisor = 1.0;
  if (mode->getPerGene())
  {
    int n = 0;
    for (int i=0; i<genes->size(); i++)
      if (genes->getGene(i).getInclude()) { n++; }
    divisor *= n;
  }
  if (mode->getPerNuc())
    divisor *= ids.size();
  
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
    
  os << "Trying score function when prediction is 0 everywhere" << endl;
  for (int j=0; j<ngenes; j++)
  {
    int length = data[j].size();
    for (int k=0; k<length; k++)
      *prediction[j][k] = 0;
  }
  getScore();
  print(os);
  os << "Trying score function when prediction is twice the data" << endl;
  for (int j=0; j<ngenes; j++)
  {
    int length = data[j].size();
    for (int k=0; k<length; k++)
      *prediction[j][k] = scale[j]->scale(*data[j][k])*2;
  }
  getScore();
  print(os);
  
}
    
  
double Score::getScore()
{
  int ngenes = genes->size();
  score = 0;
  
  (this->*scoreFunc)();
  
  for (int i=0; i<ngenes; i++)
  {   
    double tmpscore = scores[i] / divisor;
    //cerr << "scores[" << i << "] / divisor = " << scores[i] << " / " << divisor << " = " << tmpscore << endl;
    score += tmpscore;
  }
  return score;
}
    
void Score::print(ostream& os)
{
  int ngenes = genes->size();
  
  os << setw(14) << "Name";
  os << setw(14) << "Score" << endl;
  
  double total = 0.0;
  
  for (int i=0; i<ngenes; i++)
  {
    if (!genes->getGene(i).getInclude()) continue; 
    total += scores[i] / divisor;
    os << setw(14) << genes->getGene(i).getName();
    os << setw(14) << setprecision(5) << scores[i] / divisor << endl;
  }
  os << "-----------------------------" << endl;
  os << setw(14) << "Total";
  os << setw(14) << setprecision(5) << total << endl;
}
    
   


