#include "rapidxml.hpp"
#include "deep_class.h"
#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"


using namespace rapidxml;
using namespace std;
using namespace rapidjson;



void peptide_lists::xml_parse(my_parameters& my_params) {
  xml_document<> doc;

  // Read the xml file into a vector
  ifstream theFile(my_params.filename);
  vector<char> buffer((istreambuf_iterator<char>(theFile)), istreambuf_iterator<char>());
  buffer.push_back('\0');
  // Parse the buffer using the xml file parsing library into doc 
  doc.parse<0>(&buffer[0]);
  // Find root node
  xml_node<>* root_node = doc.first_node("msms_pipeline_analysis");
  
  //determine if and which prophets were used.
  bool bPProphet=false;
  bool bIProphet=false;
  xml_node<>* anal_node = root_node->first_node("analysis_summary");
  while (anal_node != NULL) {
    string tmp = anal_node->first_attribute("analysis")->value();
    if (tmp.compare("peptideprophet")==0) bPProphet=true;
    if (tmp.compare("interprophet") == 0) bIProphet = true;
    anal_node = anal_node->next_sibling("analysis_summary");
  }
  if(bIProphet) {
    cout << "iProphet deteced. Using iProbability." << endl;
    my_params.iprophet=true;
  } else if(bPProphet){
    cout << "PeptideProphet detected. Using Probability." << endl;
    my_params.iprophet=false;
  } else {
    cout << "PeptideProphet or iProphet validation is required. Please validate PSMs to obtain probability values before using SPACEPro." << endl;
    exit(24);
  }

  xml_node<>* spectrum_node = root_node->first_node("msms_run_summary");
  my_params.mzml = spectrum_node->first_attribute("base_name")->value();
  my_params.mzml += spectrum_node->first_attribute("raw_data")->value();

  total_psms_in_xml = 0;

  for (spectrum_node; spectrum_node; spectrum_node = spectrum_node->next_sibling()){
      
    //grab the spectrum query
    for (xml_node<>* sample_node = spectrum_node->first_node("spectrum_query"); sample_node; sample_node = sample_node->next_sibling()) {
      total_psms_in_xml++;
      
      //grab on the first hit of the first search result.
      xml_node<>* search_hit = sample_node->first_node("search_result")->first_node("search_hit");
      
      //go to next PSM if below probability threshold.
      xml_node<>* analysis_node = search_hit->first_node("analysis_result");
      string aType=analysis_node->first_attribute("analysis")->value();
      double probability;
      if(my_params.iprophet){
        while(aType.compare("interprophet")!=0){
          analysis_node = analysis_node->next_sibling("analysis_result");
          if(analysis_node ==NULL) {
            cout << " iProphet analysis not found in pepXML. Exiting." << endl;
            exit(12);
          }
          aType = analysis_node->first_attribute("analysis")->value();
        }
        probability = atof(analysis_node->first_node("interprophet_result")->first_attribute("probability")->value());
      } else {
        while (aType.compare("peptideprophet") != 0) {
          analysis_node = analysis_node->next_sibling("analysis_result");
          if (analysis_node == NULL) {
            cout << " PeptideProphet analysis not found in pepXML. Exiting." << endl;
            exit(12);
          }
          aType = analysis_node->first_attribute("analysis")->value();
        }
        probability = atof(analysis_node->first_node("peptideprophet_result")->first_attribute("probability")->value());
      }
      if(probability<my_params.probability) continue;

      //Keeping this PSM, so record necessary information
      dsPSM psm;
      psm.pep_seq = string(search_hit->first_attribute("peptide")->value());
      psm.pep_seq_stripped = psm.pep_seq;
      psm.scan = atoi(sample_node->first_attribute("start_scan")->value());
      psm.charge = atoi(sample_node->first_attribute("assumed_charge")->value());
      psm.pre_neutral_mass = atof(sample_node->first_attribute("precursor_neutral_mass")->value());
      psm.xml_rtime = (float)atof(sample_node->first_attribute("retention_time_sec")->value());
      psm.prev_aa = search_hit->first_attribute("peptide_prev_aa")->value()[0];
      psm.next_aa = search_hit->first_attribute("peptide_next_aa")->value()[0];
      psm.prot_seq = string(search_hit->first_attribute("protein")->value());
      psm.calc_neutral_mass = atof(search_hit->first_attribute("calc_neutral_pep_mass")->value());
      psm.tic=0;
      xml_attribute<>* cv_attr = sample_node->first_attribute("compensation_voltage");
      if(cv_attr==NULL) psm.compensation_voltage=0;
      else psm.compensation_voltage = (float)atof(cv_attr->value());
                            
      //check if peptide is proteotypic, decoy, etc.
      if (my_params.decoy.empty())  psm.decoy = false;
      else {
        psm.decoy = true;
        if (psm.prot_seq.find(my_params.decoy) == string::npos)  psm.decoy = false;
      }
      if (atoi(search_hit->first_attribute("num_tot_proteins")->value()) == 1) {
        psm.proteotypic = true;
      } else {
        psm.proteotypic = false;
        xml_node<>* ap_node = search_hit->first_node("alternative_protein");
        while (ap_node != NULL) {
          string tmp = ap_node->first_attribute("protein")->value();
          if (!my_params.decoy.empty() && tmp.find(my_params.decoy) == string::npos) {
            psm.decoy = false;
            break;
          }
          ap_node = ap_node->next_sibling("alternative_protein");
        }
      }

      //update peptide sequence if it is modified
      xml_node<>* mod_node = search_hit->first_node("modification_info");
      if (mod_node != NULL)  psm.pep_seq = mod_node->first_attribute("modified_peptide")->value();

      all_psm.push_back(psm);
    }
  }

}

bool peptide_lists::enzymatic_calc(my_parameters& my_params) {
    
    for (size_t i = 0; i < all_psm.size(); i++) {
      size_t found;
      size_t found1;
      size_t found2;
      if(my_params.enzymeSense){  //n-terminal cleavage
        found = my_params.cleave_loc.find(all_psm[i].next_aa);
        found1 = my_params.cleave_loc.find(all_psm[i].pep_seq_stripped[0]); 
        found2 = my_params.hyphen.find(all_psm[i].prev_aa);
      } else { //c-terminal cleavage
        found = my_params.cleave_loc.find(all_psm[i].prev_aa);
        found1 = my_params.cleave_loc.find(all_psm[i].pep_seq_stripped.back());
        found2 = my_params.hyphen.find(all_psm[i].next_aa);
      }
      if (found!= string::npos && (found1!=string::npos || found2!=string::npos)) {
        all_psm[i].non_enzymatic = false; 
      } else {
        //check to see if the non-specific cut is the leading Met on the protein. If so, this is probably true enzymatic.
        //ALSO: DO THIS SMARTER - LOOK IN THE DB
        if (all_psm[i].prev_aa=='M' || all_psm[i].prev_aa=='-') all_psm[i].non_enzymatic = false;
        else all_psm[i].non_enzymatic = true;
      }
    }
    return true;

}

bool peptide_lists::miss_cleave(my_parameters& my_params) {

  for (int i = 0; i < all_psm.size(); i++) {
    all_psm[i].miss_cleaves=0;
    if (my_params.enzymeSense) {  //n-terminal cleavage
      for (int j = 1; j < all_psm[i].pep_seq_stripped.size(); j++) {
        size_t found = my_params.cleave_loc.find(all_psm[i].pep_seq_stripped[j]);
        size_t found1 = my_params.anti_cleave_loc.find(all_psm[i].pep_seq_stripped[j - 1]);
        if (found != string::npos && found1 == string::npos) {
          all_psm[i].miss_cleaves++;
        }
      }
    } else {
      for (int j = 0; j < all_psm[i].pep_seq_stripped.size() - 1; j++) {
        size_t found = my_params.cleave_loc.find(all_psm[i].pep_seq_stripped[j]);
        size_t found1 = my_params.anti_cleave_loc.find(all_psm [i].pep_seq_stripped[j + 1]);
        if (found!=string::npos && found1==string::npos) {
          all_psm[i].miss_cleaves++;
        }
      }
    }
  }

  return true;

}

dsPeptide peptide_lists::convert_PSM_to_peptide(dsPSM& psm){
  dsPeptide pep;
  pep.pep_seq=psm.pep_seq;
  pep.prot_seq=psm.prot_seq;
  pep.decoy=psm.decoy;
  pep.miss_cleaves=psm.miss_cleaves;
  pep.non_enzymatic=psm.non_enzymatic;
  pep.proteotypic=psm.proteotypic;
  pep.calc_neutral_mass=psm.calc_neutral_mass;
  pep.areaXIC=psm.areaXIC;
  pep.psm_count=psm.psm_count;
  return pep;
}

//This function combines redundant PSMs into a single PSM entry with the lowest retention time.
//MH: change to keep single PSM with highest TIC.
bool peptide_lists::delete_dup() {

    size_t i;

    //All PSMs are sorted so that they are in order of sequence and charge
    sort(all_psm.begin(), all_psm.end(), compareSeqZ);
    vector<dsPSM> v;
    v.push_back(all_psm[0]);
    v.back().psm_count = 1;

    for (i = 1; i < all_psm.size(); i++) {

      //if we have a redundant PSM, increase the count and adjust RT.
      if (v.back().pep_seq.compare(all_psm[i].pep_seq) == 0 &&
        v.back().charge == all_psm[i].charge) {
        v.back().psm_count++;
        ////always store lowest RT
        //if (all_psm[i].xml_rtime < v.back().xml_rtime) v.back().xml_rtime = all_psm[i].xml_rtime;
        //store highest TIC
        if (all_psm[i].tic > v.back().tic) {
          v.back().xml_rtime = all_psm[i].xml_rtime;
          v.back().tic = all_psm[i].tic;
        }

      //otherwise, add this new PSM to the vector
      } else {
        v.push_back(all_psm[i]);
        v.back().psm_count = 1;
      }
    }

    all_psm=v; //copy the vector over the existing PSMs
   
    return true;

}


bool peptide_lists::reader() {

  //iterate over all extracted precursor signals to find the peak and compute area
  for (size_t i = 0; i < all_psm.size(); i++) {
    cleanNoise(all_psm[i].XIC);
    all_psm[i].areaXIC = calcPeakArea(all_psm[i].XIC);
  }

  //MH: Create a peptide list from the PSM list
  for (size_t i = 0; i < all_psm.size(); i++) {
    //check if there is already a peptide for this PSM
    size_t j;
    for (j = 0; j < all_peptides.size(); j++) {
      if (all_peptides[j].pep_seq.compare(all_psm[i].pep_seq) == 0) break;
    }

    if (j == all_peptides.size()) { //add new peptide
      dsPeptide pep=convert_PSM_to_peptide(all_psm[i]);
      all_peptides.push_back(pep);
    } else { //otherwise, sum areas (and psms)
      all_peptides[j].areaXIC += all_psm[i].areaXIC;
      all_peptides[j].psm_count += all_psm[i].psm_count;
    }
  }

  return true;

}


void peptide_lists::calc(metrics& my_metrics) {

  size_t i;
  int mcPep=0;
  int mcCount=0;
  for(i=0;i<all_peptides.size();i++){
    if(all_peptides[i].miss_cleaves>0){
      mcPep++;
      mcCount+= all_peptides[i].miss_cleaves;
    }
  }
  my_metrics.pep_miscleave_per_peptide=(double)mcCount/mcPep;	

}



bool peptide_lists::cleanNoise(std::vector<dsXIC>& v) {
  size_t i;
  size_t a,b; //MH: these will be your array boundaries

  //MH: apex information
  float max=0;
  size_t maxIndex=0;

  //MH: First step, find the maximum intensity (Assumed to be the peak apex)
  for(i=0;i<v.size();i++){
    if(v[i].intensity>max){
      max=v[i].intensity;
      maxIndex=i;
    }
  }

  //MH: set out threshold
  max/=10;

  //MH: Now step to the left of the apex to find your lower bound
  a=maxIndex;
  while(a>0){
    a--;
    if (v[a].intensity < max) break;
  }

  //MH: Now step to the right of the apex to find your upper bound
  b=maxIndex;
  while(b<v.size()-1){
    b++;
    if (v[b].intensity < max) break;
  }

  //MH: With bounds set, create a subset array
  vector<dsXIC> tmp;
  for(i=a;i<=b;i++){
    tmp.push_back(v[i]);
  }

  //MH: overwrite our original bloated array with the trimmed subset array
  v=tmp;

  //MH: report success
  return true;
	
}




//MH: This function replicates computing (and storing) the area between every two points.
// It then returns the sum of the area of the peak. Note that the boundaries are assumed
// to be the first and last points in the array.
float peptide_lists::calcPeakArea(std::vector<dsXIC>& v){
  float w;
  float hRect;
  float hTri;
  float total=0;
  for (size_t i = 0; i < v.size()-1; i++) {
    w = v[i + 1].rTime - v[i].rTime;
    hTri = abs(v[i + 1].intensity - v[i].intensity);
    if (v[i + 1].intensity < v[i].intensity) hRect = v[i + 1].intensity;
    else hRect = v[i].intensity;
    v[i].tot = w * hRect + w * hTri / 2;
    total+=v[i].tot;
  }
  return total;
}


void peptide_lists::print(metrics& my_metrics, my_parameters& params, string marquee) {

  char str[256];
  string s;
  vector<string> v;

  addString(v,"----------- Parameters -----------");
  sprintf(str," File: %s",params.filename.c_str());
  addString(v,str);
  sprintf(str, " Probability Cutoff: %.2lf", params.probability);
  addString(v, str);
  if(params.iprophet) sprintf(str, " Use iProphet Probability: yes");
  else sprintf(str, " Use PeptideProphet Probability: yes");
  addString(v, str);
  sprintf(str, " Retention time tolerance (min): %.2f", params.ret_time);
  addString(v, str);
  sprintf(str, " Precursor mass tolerance (ppm): %.2lf", params.ppm);
  addString(v, str);
  sprintf(str, " Enzyme cut sites: %s", params.cleave_loc.c_str());
  addString(v, str);
  sprintf(str, " Enzyme cut exceptions: %s", params.anti_cleave_loc.c_str());
  addString(v, str);
  if (params.enzymeSense) sprintf(str, " Enzyme sense: N");
  else sprintf(str, " Enzyme sense: C");
  addString(v, str);

  addString(v," ");
  addString(v, "----------- PSM STATS -----------");
  sprintf(str," Total PSMs in file: %d", my_metrics.psm_xml);
  addString(v, str);
  sprintf(str, " PSMs above probability threshold: %d", my_metrics.psm_total);
  addString(v, str);
  sprintf(str, "  Enzymatic PSMs: %d, %.2lf%%", my_metrics.psm_enzymatic, (double)my_metrics.psm_enzymatic/ my_metrics.psm_total*100);
  addString(v, str);
  sprintf(str, "  Mis-cleaved PSMs: %d, %.2lf%%", my_metrics.psm_miscleave, (double)my_metrics.psm_miscleave / my_metrics.psm_total * 100);
  addString(v, str);
  sprintf(str, "  Nonspecific PSMs: %d, %.2lf%%", my_metrics.psm_nonspecific, (double)my_metrics.psm_nonspecific / my_metrics.psm_total * 100);
  addString(v, str);

  addString(v, " ");
  addString(v, "----------- PEPTIDE STATS -----------");
  sprintf(str, " Total Unique Peptides: %d", my_metrics.pep_unique);
  addString(v, str);
  sprintf(str, " Total Peptides Quantified: %d", my_metrics.pep_count);
  addString(v, str);
  sprintf(str, "  Enzymatic Peptide Signal: %.2lf%%", my_metrics.pep_enzymatic / my_metrics.pep_total * 100);
  addString(v, str);
  sprintf(str, "  Mis-cleaved Peptide Signal: %.2lf%%", my_metrics.pep_miscleave / my_metrics.pep_total * 100); 
  addString(v, str);
  sprintf(str, "  Nonspecific Peptide Signal: %.2lf%%", my_metrics.pep_nonspecific / my_metrics.pep_total * 100); 
  addString(v, str);
  sprintf(str, " Average number of mis-cleavages among mis-cleaved peptides: %.2lf", my_metrics.pep_miscleave_per_peptide);
  addString(v, str);

  addString(v, " ");
  addString(v, "----------- PROTEIN STATS -----------");
  sprintf(str, " Total Proteins from Proteotypic Peptides: %d", my_metrics.prot_count);
  addString(v, str);
  sprintf(str, "  Average Enzymatic Protein Signal: %.2lf%%", my_metrics.prot_avg_enzymatic);
  addString(v, str);
  sprintf(str, "  Average Mis-cleaved Protein Signal: %.2lf%%", my_metrics.prot_avg_miscleave);
  addString(v, str);
  sprintf(str, "  Average Nonspecific Protein Signal: %.2lf%%", my_metrics.prot_avg_nonspecific);
  addString(v, str);

  string rf=params.filename;
  rf+=".sp.txt";
  FILE* f=fopen(rf.c_str(),"wt");
  fprintf(f,"%s\n",marquee.c_str());
  for(size_t i=0;i<v.size();i++){
    cout << v[i] << endl;
    fprintf(f,"%s\n",v[i].c_str());
  }
  fclose(f);

  
}

void peptide_lists::json(string fn) {
    
  StringBuffer s;
  PrettyWriter<StringBuffer> writer(s);
     
  writer.StartObject();
  writer.Key("Proteins");
  writer.StartArray();
  for (int i = 0; i <all_proteins.size(); i++) {
    writer.StartObject();
    writer.Key("Protein Name");
            
    writer.String(all_proteins[i].prot_seq.c_str());
    writer.Key("Total Intensity");
    writer.Double(all_proteins[i].total);
    writer.Key("Total Enzymatic Peptide Intensity");
    writer.Double(all_proteins[i].sumEnz);
    writer.Key("Total Mis-cleaved Peptide Intensity");
    writer.Double(all_proteins[i].sumMiss);
    writer.Key("Total Nonspecific Peptide Intensity");
    writer.Double(all_proteins[i].sumNonSp);
    writer.Key("Percent Enzymatic");
    writer.Double(all_proteins[i].sumEnz / all_proteins[i].total * 100);
    writer.Key("Percent Mis-Cleaved");
    writer.Double(all_proteins[i].sumMiss / all_proteins[i].total * 100);
    writer.Key("Percent Nonspecific");
    writer.Double(all_proteins[i].sumNonSp / all_proteins[i].total*100);
    writer.Key("Peptides");
    writer.StartArray();
    for (int j = 0; j < all_proteins[i].peptides.size(); j++) {
      size_t index = all_proteins[i].peptides[j];
      writer.StartObject();
      writer.Key("Sequence");
      writer.String(all_peptides[index].pep_seq.c_str());
      writer.Key("Abundance");
      writer.Double(all_peptides[index].areaXIC);
      writer.Key("PSM Count");
      writer.Int(all_peptides[index].psm_count);
      writer.Key("Mis-cleaved");
      if (all_peptides[index].miss_cleaves > 0) writer.Bool(true);
      else writer.Bool(false);
      writer.Key("Nonspecific");
      if (all_peptides[index].non_enzymatic) writer.Bool(true);
      else writer.Bool(false);
      writer.EndObject();
    }
    writer.EndArray();
    writer.EndObject();
  }
  writer.EndArray();
  writer.EndObject();

  FILE* f=fopen(fn.c_str(),"wt");
  fprintf(f,"%s",s.GetString());
  fclose(f);   
 
}

bool peptide_lists::prot_stats() {

  //Sort peptides by protein name - note that peptides cannot be sorted ever again to maintain indexes.
  sort(all_peptides.begin(), all_peptides.end(), compareProt);

  //Make the first proteotypic protein from the first peptide
  size_t i=0;
  while(!all_peptides[i].proteotypic) i++;
  dsProtein prot;
  prot.prot_seq=all_peptides[i].prot_seq;
  all_proteins.push_back(prot);
  all_proteins.back().peptides.push_back(i);

  //iterate over all remaining peptides, combining them into proteins to store in the protein array
  for (i = i+1; i < all_peptides.size(); i++) {

    //skip peptides that were not quantified
    if(all_peptides[i].areaXIC==0) continue;

    //skip non-proteotypic peptides
    if (!all_peptides[i].proteotypic) continue;

    //if peptide belongs to current protein, add it to that protein's peptide list
    if(all_peptides[i].prot_seq.compare(all_proteins.back().prot_seq)==0){
      all_proteins.back().peptides.push_back(i);
    
    //otherwise, start a new protein
    } else {
      prot.prot_seq=all_peptides[i].prot_seq;
      all_proteins.push_back(prot);
      all_proteins.back().peptides.push_back(i);
    }

  }


  //iterate over protein list to sum up all peptide signals
  //this could be done above, but is done here for ease of reading the code
  for(i=0;i<all_proteins.size();i++){
    all_proteins[i].sumMiss=0;
    all_proteins[i].sumNonSp=0;
    all_proteins[i].sumEnz=0;
    all_proteins[i].sumPSM=0;
    all_proteins[i].total=0;
    for(size_t j=0;j<all_proteins[i].peptides.size();j++){
      size_t index=all_proteins[i].peptides[j];
      if(all_peptides[index].miss_cleaves>0) all_proteins[i].sumMiss += all_peptides[index].areaXIC;
      if (all_peptides[index].non_enzymatic) all_proteins[i].sumNonSp += all_peptides[index].areaXIC;
      if(all_peptides[index].miss_cleaves==0 && !all_peptides[index].non_enzymatic) all_proteins[i].sumEnz+= all_peptides[index].areaXIC;
      all_proteins[i].total += all_peptides[index].areaXIC;
      all_proteins[i].sumPSM += all_peptides[index].psm_count;
    }
  }

  return true;
}
