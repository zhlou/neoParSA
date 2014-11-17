/*********************************************************************************
*                                                                                *
*     conc.cpp                                                                   *
*                                                                                *
*     Contains class which holds a data value at each nucleus, whether           *
*     it is the rate data or the TF concentrations. These can accept             *
*     scale factors to adjust between experiements                               *
*                                                                                *
*********************************************************************************/

# include "conc.h"

# define foreach_ BOOST_FOREACH

# include <boost/foreach.hpp>
# include <cstdlib>

/******************************   ConcContainer    ***********************************/

/*    Constructors    */

ConcContainer::ConcContainer() {}

ConcContainer::ConcContainer(ptree& pt, scale_factors_ptr p, string node_name)
{
  read(pt, p, node_name);
}

/*    Setters   */



/*    Getters   */



double& ConcContainer::getConcByIndex(int idx, const string& n)
{
  return scaled_conc[n][idx];
}

// included in case if forget which order
double& ConcContainer::getConcByIndex(const string& n, int idx)
{
  return scaled_conc[n][idx];
}

double& ConcContainer::getConcByID(int id, const string& n, bool scaled)
{
  if (!hasName(n))
  {
    cerr << "ERROR: could not find gene with name " << n << " in rate data" << endl;
    exit(1);
  }
  
  int i;
  int nids = ids.size();
  for (i=0; i<nids; i++)
    if (ids[i] == id) break;
  
  if (scaled)
    return scaled_conc[n][i];
  else
    return conc[n][i];
}

double& ConcContainer::getConcByID(const string& n, int id, bool scaled)
{
  return getConcByID(id, n, scaled);
}
  
int ConcContainer::size() { return ids.size(); } 
  
/*    Methods   */

bool ConcContainer::hasID(int id)
{
  int i;
  int nids = ids.size();
  for (i=0; i<nids; i++)
  {
    if (ids[i] == id)
      return true;
  }
  return false;
}

bool ConcContainer::hasName(const string & n)
{
  int nnames = names.size();
  for (int i=0; i<names.size(); i++)
  {
    if (names[i] == n)
      return true;
  }
  return false;
}

void ConcContainer::scale_data()
{
  int nids   = ids.size();
  int nnames = names.size();

  for (int i=0; i<nnames; i++)
  {
    string& name = names[i];
    
    scale_factor_ptr scale_factor    = scales_map[name];
    vector<double>&  tmp_conc        = conc[name];
    vector<double>&  tmp_scaled_conc = scaled_conc[name];
    
    for (int j=0; j<nids; j++)
      tmp_scaled_conc[j] = scale_factor->scale(tmp_conc[j]);
  }
}
      
/*    I/O    */


void ConcContainer::read(ptree& pt, scale_factors_ptr p, string node_name)
{
  scales = p;
  
  // first we read through and find all the IDs and nuclei we have
  double NA = numeric_limits<double>::infinity();
  
  string xmlattr("<xmlattr>.");
  
  int id;
  string name;
  double value;
  
  ptree& conc_node = pt.get_child(node_name);
  
  foreach_(ptree::value_type const& nuc, conc_node)
  {
    if (nuc.first != "Nucleus") continue;
    
    id = nuc.second.get<int>("<xmlattr>.ID");
    if (!hasID(id)) ids.push_back(id);
    const ptree & attrs = nuc.second.get_child("<xmlattr>");
    foreach_(ptree::value_type const& attr, attrs)
    {
      if (attr.first != "ID")
      {
        name = attr.first;
        if (!hasName(name)) names.push_back(name);
      }
    }   
  }
  // now that we have the names and ID's, go access them specifically
  foreach_(ptree::value_type const& nuc, conc_node)
  {
    if (nuc.first != "Nucleus") continue;
    
    int nnames = names.size();
    for (int i=0; i<nnames; i++)
    {
      value = nuc.second.get<double>(xmlattr + names[i], NA);
      if (value == NA)
        cerr << "WARNING: missing rate value for " << names[i] << " on nucleus " << nuc.second.get<int>("<xmlattr>.ID") << ", setting to NA" << endl;
      conc[names[i]].push_back(value);
    }
  }
  // now get scale factors and set scaled conc, starting with setting it to default in case node is missing
  int nnames = names.size();
  int nids   = ids.size();
  for (int i=0; i<nnames; i++)
  {
    string scalename("default");
    scales_map[names[i]] = scales->getScaleFactor(scalename);
    scaled_conc[names[i]] = conc[names[i]];
  }
  
  foreach_(ptree::value_type const& nuc, conc_node)
  {
    if (nuc.first != "ScaleFactor") continue;
    int nnames = names.size();
    for (int i=0; i<nnames; i++)
    {
      string scalename = nuc.second.get<string>(xmlattr + names[i], "default");
      scales_map[names[i]] = scales->getScaleFactor(scalename);
    }
  }
  
  scale_data();
}

void ConcContainer::write(ostream & os, string node_name)
{
  ptree pt;
  write(pt, node_name);
  boost::property_tree::xml_writer_settings<char> settings(' ', 2);
  write_xml_element(os, basic_string<ptree::key_type::value_type>(), pt, -1, settings);  
}

void ConcContainer::write(ptree & pt, string node_name)
{
  int i, j;
  int maxname = 0;
  int idsize = ids.size();
  int namesize = names.size();
  
  // get the largest name for pretty printing
  for (i=0; i<namesize; i++)
  {

    maxname = max(maxname, (int) names[i].size());
  }
  
  stringstream name;
  stringstream value;
  
  ptree & conc_node = pt.add(node_name,"");
  ptree & scale = conc_node.add("ScaleFactor","");
  
  for (i=0; i<namesize; i++)
  {
    name.str("");
    value.str("");
    name << names[i];
    scale.put("<xmlattr>."+name.str(),scales_map[names[i]]->getName());
  }
  
  for (i = 0; i<idsize; i++)
  {
    ptree & nuc = conc_node.add("Nucleus","");
    nuc.put("<xmlattr>.ID", ids[i]);
    for (j=0; j<namesize; j++)
    {
      name.str("");
      value.str("");
      name << names[j];
      value << setw(6) << setprecision(4) << conc[names[j]][i];
      nuc.put("<xmlattr>."+name.str(),value.str());
    }
  }
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

