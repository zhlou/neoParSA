/*********************************************************************************
*                                                                                *
*     TF.cpp                                                                     *
*                                                                                *
*     Contains class definition for Transcription Factos                         *
*                                                                                *
*********************************************************************************/

# include "TF.h"
# include "utils.h"

# include <cmath>
# include <iomanip>
# include <sstream>
# include <boost/property_tree/xml_parser.hpp>
# include <boost/property_tree/ptree.hpp>
# include <boost/foreach.hpp>



using namespace std;
using boost::property_tree::ptree;

# define foreach_ BOOST_FOREACH


/*********************************   TF    **************************************/

/*      Constructors    */

TF::TF():
  kmax(param_ptr(new Parameter)),
  threshold(param_ptr(new Parameter)),
  lambda(param_ptr(new Parameter))
{}

TF::TF(ptree& pt):
  kmax(param_ptr(new Parameter)),
  threshold(param_ptr(new Parameter)),
  lambda(param_ptr(new Parameter))
{
  set(pt);
}

void TF::set(ptree& pt)
{
  maxscore = 0.0;
  Nbehavior = 1;
  
  if (!pt.count("PWM")) error("TF::TF(ptree& pt) TF node malformed");

  tfname    = pt.get<string>("<xmlattr>.name");

  bsize     = pt.get<int>("<xmlattr>.bsize");
  gc        = pt.get<double>("<xmlattr>.gc", 0.5);
  
  // read in coefficients
  ptree& coefs_node = pt.get_child("Coefficients");
  foreach_(ptree::value_type const& coef_node, coefs_node)
  {
    if (coef_node.first != "coef") continue;
    
    param_ptr coef(new Parameter());
    coef->read( (ptree&) coef_node.second);
    coef->setParamName("coef");
    coef->setTFName(tfname);
    coefs.push_back(coef);
  }
  
  ptree& kmax_node = pt.get_child("kmax");
  kmax->read(kmax_node);
  kmax->setParamName("kmax");
  kmax->setTFName(tfname);
  
  ptree& threshold_node = pt.get_child("threshold");
  threshold->read(threshold_node);
  threshold->setParamName("threshold");
  threshold->setTFName(tfname);
  
  ptree& lambda_node = pt.get_child("lambda");
  lambda->read(lambda_node);
  lambda->setParamName("lambda");
  lambda->setTFName(tfname);
  
  ptree& pwm_node = pt.get_child("PWM");
  
  readPWM(pwm_node);
  
}

void TF::readPWM(ptree& pwm_node)
{
  // first read in the pwm
  string         token;
  vector<double> line;
  
  foreach_(ptree::value_type const& v, pwm_node)
  {
    if (v.first == "position")
    {
    line.clear();
    
    stringstream s(v.second.data());
     
    while( getline(s, token, ';'))
    {
      line.push_back(atof(token.c_str()));
    }
    energy.push_back(line);
    //maxscore += tmax; 
    }
  }
  
  string type = pwm_node.get<string>("<xmlattr>.type");
  if (type == "PCM")
  {
    double pseudo = pwm_node.get<double>("<xmlattr>.pseudo", 1.0);
    PCM2PFM(pseudo);
    PFM2PSSM();
    PSSM2BEM();
  } else if (type == "PFM")
  {
    PFM2PSSM();
    PSSM2BEM();
  } else if (type == "PSSM")
  {
    PSSM2BEM();
  } else if (type == "BEM")
  {
    maxscore = 0.0;
  } else
  {
    cerr << "ERROR: getPWM did not recognize pwm of type " << type << endl;
    exit(1);
  }
  setNscore();
  pwmlen = energy.size();
}
    

void TF::PCM2PFM(double pseudo)
{
  // we have counts. We need to add pseudo count and divide rows
  int pwmlen = energy.size();
  double gc_adjusted_count;
  
  for (int i=0; i<pwmlen; i++)
  {
    double rowsum = 0;
    
    for (int j=0; j<4; j++)
    {
      
      if (j==0 || j==3) // we have an A or T
        gc_adjusted_count = pseudo*(1-gc)/2;
      else             // we have a G or C
        gc_adjusted_count = pseudo*(gc)/2;

      energy[i][j] += gc_adjusted_count;

      rowsum       += energy[i][j];
    }
    for (int j=0; j<4; j++)
      energy[i][j] /= rowsum;
  }
}

void TF::PFM2PSSM()
{
  int pwmlen = energy.size();
  double bkgd;
  
  for (int i=0; i<pwmlen; i++)
  {
    for (int j=0; j<4; j++)
    {
      if (j==0 || j==3)      // we have an A or T
        bkgd = (1-gc)/2;
      else if (j==1 || j==2) // we have a G or C
        bkgd = (gc)/2;
      else
        bkgd = min((1-gc)/2, gc/2);
      
      energy[i][j] = log( energy[i][j]/bkgd );
    }
  }
}

void TF::PSSM2BEM()
{
  maxscore = 0.0;
  int pwmlen = energy.size();
  for (int i=0; i<pwmlen; i++)
  {
    double tmax = 0.0;
    for (int j=0; j<4; j++)
      tmax = max(tmax, energy[i][j]);
    for (int j=0; j<4; j++)
      energy[i][j] -= tmax;
    maxscore += tmax;
  }
  threshold->set(threshold->getValue() - maxscore);
  threshold->setLimits(threshold->getLimLow()  - maxscore,
                       threshold->getLimHigh() - maxscore);
  maxscore = 0.0;
}
    
void TF::setNscore()
{
  int pwmlen = energy.size();
  if (Nbehavior == 1)
  {
    for (int i=0; i<pwmlen; i++)
    {
      double minscore = energy[i][0];
      for (int j=0; j<4; j++)
        minscore = min(minscore, energy[i][j]);
      energy[i].push_back(minscore);
    }
  }
}
      
    
  
/*    Setters   */

void TF::setName(string n) { tfname = n; }

void TF::setPwm(vector< vector<double> > p)
{
  int i,j;
  double minscore = 0.0;
  double tmax = 0.0;
  maxscore = 0.0;
  energy = p;
  pwmlen = p.size();
  Nbehavior = 1; // default to conservative N;
  for (i=0; i<pwmlen; i++)
  {
    minscore = tmax = 0.0;
    for (j=0; j<4; j++)
    {
      minscore = min(minscore, energy[i][j]);
      tmax = max(tmax, energy[i][j]);
    }
    energy[i].push_back(minscore);
    maxscore += tmax;
  }  
}
 
void TF::setKmax(double k) { kmax->set(k); }

void TF::setLambda(double l) { lambda->set(l); }

void TF::setThreshold(double t) {threshold->set(t); }

void TF::setBindingSize(int b) { bsize = b; }
  

/*    Getters   */

const string& TF::getName() const { return tfname; }

double TF::getThreshold() { return threshold->getValue(); }

double TF::getKmax() { return kmax->getValue(); }

double& TF::getCoef() { return coefs[0]->getValue(); }

double TF::getModifiedCoef() 
{
  if (coefs.size() > 1)
    return coefs[0]->getValue();
  else
    return 0;
}

bool TF::neverActivates()
{
  int ncoefs = coefs.size();
  for (int i=0; i<ncoefs; i++)
  {
    if (coefs[i]->getValue() > 0)
      return false;
  }
  return true;
}

bool TF::neverQuenches()
{
  int ncoefs = coefs.size();
  for (int i=0; i<ncoefs; i++)
  {
    if (coefs[i]->getValue() < 0)
      return false;
  }
  return true;
}

vector<double> TF::getCoefs()
{
  vector<double> out;
  int ncoefs = coefs.size();
  out.resize(ncoefs);
  for (int i=0; i<ncoefs; i++)
    out[i] = coefs[i]->getValue();
  return out;
}

double TF::getLambda() { return lambda->getValue(); }

double TF::getMaxScore() const { return maxscore; }

int TF::getBindingSize() const { return bsize; }

void TF::getParameters(param_ptr_vector& p)
{
  if (kmax->isAnnealed())
    p.push_back(kmax);
  if (threshold->isAnnealed())
    p.push_back(threshold);
  if (lambda->isAnnealed())
    p.push_back(lambda);
  
  int ncoefs = coefs.size();
  for (int i=0; i<ncoefs; i++)
  {
    if (coefs[i]->isAnnealed())
      p.push_back(coefs[i]);
  }
}

/*bool TF::checkCoops(TF* t)
{
  int ncoops = coops.size();
  for (int i=0; i<ncoops; i++)
  {
    if (coops[i].first == t)
      return true;
  }
  return false;
}*/

// check if cooporates with orientation
bool TF::checkCoops(TF* t, char o1, char o2) 
{
  int ncoops = coops.size();
  for (int i=0; i<ncoops; i++)
  {
    if (coops[i].first == t)
    {
      coop_ptr coop = coops[i].second;
      if (o1 == 'F')
      {
        if      (o2 == 'F' && coop->getHT())
          return true;
        else if (o2 == 'R' && coop->getHH())
          return true;
      } 
      else //(o1 == 'R')
      {
        if      (o2 == 'F' && coop->getTT())
          return true;
        else if (o2 == 'R' && coop->getHT())
          return true;
      }
    }
  }
  return false;
}

coop_ptr TF::getCoop(TF* t)
{
  int ncoops = coops.size();
  for (int i=0; i<ncoops; i++)
  {
    if (coops[i].first == t)
      return coops[i].second;
  }
  coop_ptr void_coop;
  return void_coop;
}
  

/*    Methods   */

void TF::subscore(const vector<int> & s, double * out)
{
  int i,j,k;
  out[0]=0.0;
  out[1]=0.0;
  for (i=0, j=(pwmlen-1); i<pwmlen; i++, j--)
  {
    // forward sequence, i iterates forward
    k = s[i];
    out[0] += energy[i][k];
    
    // reverse sequence, j iterates back over subseq, 3-n give complement
    if (s[j] !=4) 
      k = 3-s[j]; 
    out[1] += energy[i][k];
  }
  out[2] = max(out[0],out[1]);
}



void TF::score(const vector<int>& s, TFscore &t)
{
  /* I do something a little unorthodox here. I want to control for boundary
  effects, so I pad the beginning and end of the sequence with Ns and score
  using the boundary behavior. Unlike the transc code, I report the score for the
  middle of the binding site rather than the m position. This means my output is
  of the same length as the sequence and it is the same length for every factor */
  
  int i, j, k;
  int mdist, ndist;
  int start, end;
  int len;
  vector<int> sub(pwmlen);
  
  ndist = pwmlen/2; // returns floor for middle position
  mdist = pwmlen-ndist;
  
  double * out = new double[3];
  len = s.size();

  t.fscore.resize(len);
  t.rscore.resize(len);
  t.mscore.resize(len);
  //t.tfname = tfname;
  
  for (i=0; i<len; i++)
  {
    start = i - mdist;
    end = i + ndist;
    // j = index in s
    // j = index in pwm
    for (j=start, k=0; j<end; j++, k++)
    {
      if (j<0 || j>=len)
      {
        sub[k]=4;
      } else 
      {
        sub[k]=s[j];
      }
    }
    subscore(sub,out);
    t.fscore[i] = out[0];
    t.rscore[i] = out[1];
    t.mscore[i] = out[2];
  } 
  delete[] out;
}
    
TFscore TF::score(const vector<int>& s)
{
  TFscore t;
  score(s,t);
  return(t);
}

TFscore TF::score(const string & s)
{
  
  vector<int> out(s.size());
  out = string2int(s);
  return(score(out));
}
  

/*    print   */

void TF::print(ostream& os) 
{
  int w;
  int i,j;
  w = 8;
  os << setprecision(4);
  os << "Transcription Factor" << endl
     << "name: "   << tfname   << endl
     << "kmax: "   << kmax->getValue() << endl
     << "lambda: " << lambda->getValue()   << endl
     << "bsize: "  << bsize    << endl
     << "threshold: " << threshold->getValue() << endl
     << "maxscore: " << maxscore << endl
     << setw(w) << ""
     << setw(w) << "A"
     << setw(w) << "C"
     << setw(w) << "G"
     << setw(w) << "T"
     << setw(w) << "N"
     << endl;
     
  for (i=0; i<pwmlen; i++)
  {
    os << setprecision(3);
    os << setw(w) << i+1;
    for (j=0; j<5; j++)
    {
      os << setw(w) << energy[i][j];
    }
    os << endl;
  }
}

void TF::write(ostream& os) 
{
  stringstream tmp;
  tmp << setprecision(3);
  
  ptree pt; // boost property tree
  write(pt);
  
  boost::property_tree::xml_writer_settings<char> settings(' ', 2);
  write_xml_element(os, basic_string<ptree::key_type::value_type>(), pt, -1, settings);
}

void TF::write(ptree& tfsnode) 
{
  int i;
  int w = 12;
  stringstream tmp;
  tmp << setprecision(6);
  
  ptree & tfnode = tfsnode.add("TF","");
  tfnode.put("<xmlattr>.name",      tfname);
  
  ptree & kmax_node      = tfnode.add("kmax          ","");
  kmax->write(kmax_node);
  
  ptree & threshold_node = tfnode.add("threshold     ","");
  threshold->write(threshold_node);

  ptree & lambda_node    = tfnode.add("lambda        ","");
  lambda->write(lambda_node);
  
  ptree& coefs_node = tfnode.add("Coefficients", "");
  int ncoefs = coefs.size();
  for (int i=0; i<ncoefs; i++)
  {
    ptree& coef_node = coefs_node.add("coef","");
    coefs[i]->write(coef_node);
  }
  
  tfnode.put("<xmlattr>.bsize",   bsize);
  tfnode.put("<xmlattr>.gc",      gc);
  tfnode.put("<xmlattr>.include", "true");
  ptree & pwmnode = tfnode.add("PWM","");
  pwmnode.put("<xmlattr>.type","BEM");
  tmp.str("");
  tmp << setw(w+1) << "A"
      << setw(w+1) << "C"
      << setw(w+1) << "G"
      << setw(w+1) << "T";
  pwmnode.add("base   ", tmp.str());
  for (i=0; i<pwmlen; i++)
  {
    tmp.str("");
    tmp << setw(w) << energy[i][0] << ";"
        << setw(w) << energy[i][1] << ";"
        << setw(w) << energy[i][2] << ";"
        << setw(w) << energy[i][3];
    pwmnode.add("position",tmp.str());
  }
}
  
    

/*********************************   TFContainer    *********************************/


/*    Constructor   */

TFContainer::TFContainer() {}

TFContainer::TFContainer(ptree& pt) 
{
  add(pt);
}
    

/*    Getters   */

TF& TFContainer::getTF(const string& n)
{
  int i;
  int ntfs = tfs.size();
  
  for (i=0; i<ntfs; i++)
  {
    if(tfs[i]->getName() == n) // they are equal
      return(*tfs[i]);                // return this tf
  }
  // if we are here TF was not found
  cerr << "ERROF: getTF() TF with name " << n << " not found" << endl;
  exit(1);
}

int TFContainer::size() const { return(tfs.size()); }

TF& TFContainer::getTF(int index) { return *tfs[index]; }

tf_ptr TFContainer::getTFptr(int index) { return tfs[index]; }

tf_ptr TFContainer::getTFptr(const string& n)
{
  int i;
  int ntfs = tfs.size();
  
  for (i=0; i<ntfs; i++)
  {
    if(tfs[i]->getName() == n) // they are equal
      return(tfs[i]);                // return this tf
  }
  // if we are here TF was not found
  cerr << "ERROF: getTF() TF with name " << n << " not found" << endl;
  exit(1);
}
  
void TFContainer::getParameters(param_ptr_vector& p)
{
  int ntfs = tfs.size();
  for (int i=0; i<ntfs; i++)
    tfs[i]->getParameters(p);

}
  
  
/*  Setters   */

void TFContainer::setCoops(coops_ptr c)
{
  int ntfs = tfs.size();
  for (int i=0; i<ntfs; i++)
  {
    coop_pairs cur_coops = c->getCoops(tfs[i]->getName());
    // string comparisons are slow, so we want the final object to be TF pointers
    vector< pair<TF*,coop_ptr> > tmp_pairs;
    int ncur_coops = cur_coops.size();
    for (int j=0; j<ncur_coops; j++)
    {
      TF& tf_ref = getTF(cur_coops[j].first);
      TF* tf = &tf_ref;
      tmp_pairs.push_back(pair<TF*,coop_ptr>(tf, cur_coops[j].second));
    }
    tfs[i]->setCoops(tmp_pairs);
  }
}

void TFContainer::setCoeffects(coeffects_ptr c)
{
  int ntfs = tfs.size();
  for (int i=0; i<ntfs; i++)
  {
    coeffect_pairs cur_coefs = c->getTargets(tfs[i]->getName());
    // string comparisons are slow, so we want the final object to be TF pointers
    vector< pair<TF*, coeffect_ptr> > tmp_pairs;
    int ncur_coefs = cur_coefs.size();
    for (int j=0; j<ncur_coefs; j++)
    {
      TF& tf_ref = getTF(cur_coefs[j].first);
      TF* tf     = &tf_ref;
      tmp_pairs.push_back(pair<TF*,coeffect_ptr>(tf,cur_coefs[j].second));
    }
    tfs[i]->setCoeffects(tmp_pairs);
  }
} 
  
  
// can handle getting a TF node or TFs node for single or multiple add
void TFContainer::add(ptree& pt)
{
  if (pt.count("TFs")) // we are at a parent of TFs
  {
    ptree& tfs_node = pt.get_child("TFs");
    foreach_(ptree::value_type const& tf, tfs_node)
    {
      if (tf.first == "TF")
      {
        ptree& tf_node = (ptree&) tf.second;
        if (tf_node.get<bool>("<xmlattr>.include"))
        {
          tf_ptr t(new TF(tf_node) );
          tfs.push_back(t);
        }
      }
    }
  }
  else if (pt.count("TF")) // we are at the TFs node
  {
    foreach_(ptree::value_type const& tf, pt)
    {

      if (tf.first == "TF")
      {
        ptree& tf_node = (ptree&) tf.second;
        if (tf_node.get<bool>("<xmlattr>.include"))
        {
          tf_ptr t(new TF(tf_node) );
          tfs.push_back(t);
        }
      }
    }
  } 
  else // this is a TF node
  {
    if (pt.get<bool>("<xmlattr>.include"))
    {
      tf_ptr t(new TF(pt) );
      tfs.push_back(t);
    }
  }
} 

void TFContainer::add(tf_ptr t) {tfs.push_back(t);}

/*    Output    */
void TFContainer::write(ostream& os) const
{
  ptree pt;
  write(pt);
  
  boost::property_tree::xml_writer_settings<char> settings(' ', 2);
  write_xml_element(os, basic_string<ptree::key_type::value_type>(), pt, -1, settings);
}

void TFContainer::write(ptree& pt) const
{
  int i;
  int ntfs = tfs.size();
  ptree & tfsnode = pt.add("TFs","");
  for (i=0; i<ntfs; i++)
    tfs[i]->write(tfsnode);
  
}
   

