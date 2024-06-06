/*
Program       : ProDataComponent                                                    
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 
SVN Info      : $Id: ProDataComponent.cpp 8332 2020-12-18 03:52:51Z real_procopio $


Copyright (C) 2003 Andrew Keller

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Andrew Keller
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA
akeller@systemsbiology.org
*/

#include "ProDataComponent.h"
#include "Parsers/mzParser/mzParser.h"

using namespace mzParser;

ProDataComponent::ProDataComponent(int seq, int pk, int data, int xml_ind, int ind, Boolean heavy2light) {
  xpress_ = NULL;
  seq_ = seq;
  peak_ = pk;
  data_ = data;
  xml_index_ = xml_ind;
  data_index_ = ind;
  msms_run_idx_ = -1;

  //pepdata_ = NULL;
  result_ = NULL;
  basename_index_ = -1;
  pepproph_timestamp_index_ = -1;
  iproph_timestamp_index_ = -1;
  database_index_ = -1;
  asap_timestamp_index_ = -1;
  //xpress_ = -1.0;
  heavy2light_ = heavy2light;
  masstype_ = 0;
  fragmasstype_ = 0;
}

ProDataComponent::ProDataComponent(int seq, int pk, int data, int xml_ind, int ind, Boolean heavy2light, int msms_run_idx) {
  xpress_ = NULL;
  seq_ = seq;
  peak_ = pk;
  data_ = data;
  xml_index_ = xml_ind;
  data_index_ = ind;
  msms_run_idx_ = msms_run_idx;

  //pepdata_ = NULL;
  result_ = NULL;
  basename_index_ = -1;
  pepproph_timestamp_index_ = -1;
  iproph_timestamp_index_ = -1;
  database_index_ = -1;
  asap_timestamp_index_ = -1;
  //xpress_ = -1.0;
  heavy2light_ = heavy2light;
  masstype_ = 0;
  fragmasstype_ = 0;
}

ProDataComponent::~ProDataComponent() {
  if(xpress_ != NULL)
    delete xpress_;
  if(result_ != NULL)
    delete result_;
}

void ProDataComponent::enter(pepDataStrct pepdata, SearchResult* result, int basename_index,
			     int pepproph_timestamp_index, int iproph_timestamp_index,
			     int database_index, int asap_timestamp_index, prodatacomponent_struct comp_data, long scan) {

  //double xpress, double asap_mean, double asap_error, 
  //			     int asap_ind, char* score_summary) {
  pepdata_ = pepdata;
  result_ = result;
  basename_index_ = basename_index;
  pepproph_timestamp_index_ = pepproph_timestamp_index;
  iproph_timestamp_index_ = iproph_timestamp_index;
  database_index_ = database_index;
  asap_timestamp_index_ = asap_timestamp_index;
  xpress_ = comp_data.xpressratio;
  asap_mean_ = comp_data.asap_mean;
  asap_error_ = comp_data.asap_err;
  asap_inv_mean_ = comp_data.asap_inv_mean;
  asap_inv_error_ = comp_data.asap_inv_err;
  asap_index_ = comp_data.asapratio_index;
  score_summary_ = new char[strlen(comp_data.score_summary)+1];
  strcpy(score_summary_, comp_data.score_summary);
  LightFirstScan_ = comp_data.lightfirstscan;
  LightLastScan_ = comp_data.lightlastscan;
  HeavyFirstScan_ = comp_data.heavyfirstscan;
  HeavyLastScan_ = comp_data.heavylastscan;
  LightMass_ = comp_data.lightmass;
  HeavyMass_ = comp_data.heavymass;
  MassTol_ = comp_data.masstol;
  bXpressLight_ = comp_data.xpresslight;
  scan_ = scan;
  xpress_ = new char[strlen(comp_data.xpressratio)+1];
  strcpy(xpress_, comp_data.xpressratio);
}

Boolean ProDataComponent::match(int seq, int pk, int data) {
  return seq_ == seq && peak_ == pk && data_ == data;
}

Boolean ProDataComponent::xml_match(int xml_ind) {
  return xml_index_ == xml_ind;
}


// for html display
char* ProDataComponent::display(Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times) {

  if(basename_index_ >= basenames->length()) {
    cout << "error with basenames" << endl;
    exit(1);
  }
  if(pepproph_timestamp_index_ >= pepproph_times->length()) {
    cout << "error with pepproph_times" << endl;
    exit(1);
  }
  if(iproph_timestamp_index_ >= iproph_times->length()) {
    cout << "error with iproph_times" << endl;
    exit(1);
  }
  if(asap_timestamp_index_ >= asap_times->length()) {
    cout << "error with asap_times" << endl;
    exit(1);
  }
  if(database_index_ >= dbs->length()) {
    cout << "error with databases" << endl;
    exit(1);
  }

  char* output = new char[100000];

  strcpy(output, "basename: ");
  strcat(output, (*basenames)[basename_index_]);
  strcat(output, "\n");
  
  if (pepproph_timestamp_index_ >= 0) {
    strcat(output, "pepproph: ");
    strcat(output, (*pepproph_times)[pepproph_timestamp_index_]);
    strcat(output, "\n");
  }
  if (iproph_timestamp_index_ >= 0) {
    strcat(output, "iproph: ");
    strcat(output, (*iproph_times)[iproph_timestamp_index_]);
    strcat(output, "\n");
  }

  strcat(output, "db: ");
  strcat(output, (*dbs)[database_index_]);
  strcat(output, "\n");
  
  strcat(output, "asap: ");
  strcat(output, (*asap_times)[asap_timestamp_index_]);
  strcat(output, "\n");
  
  if(result_ != NULL) {
    
    strcat(output, "spec: ");
    strcat(output, result_->spectrum_);
    strcat(output, "\n");
    
    strcat(output, "engine: ");
    strcat(output, result_->getName());
    strcat(output, "\n");
  }
  else {
    cout << "error: null result" << endl;
    exit(1);
  }

  return output;
}

void ProDataComponent::setMassType(int val, int fragval) {
  masstype_ = val;
  fragmasstype_ = fragval;
}


char* ProDataComponent::getColoredPeptide(const char* peptide, const char* labeled_aas, const char* starttag, const char* endtag) {
  // first count how large will be
  size_t index = 0;
  Boolean color = False;
  int k;
  for(k = 0; peptide[k]; k++) {
    if(! color && strchr(labeled_aas, peptide[k]) != NULL) { // turn on
      color = True;
      index += strlen(starttag);
    }
    else if(color && peptide[k] >= 'A' && peptide[k] <= 'Z' && strchr(labeled_aas, peptide[k]) == NULL) { // turn off
      color = False;
      index += strlen(endtag);
    }
      
    index++; // for the aa

  } // next aa of pep
  if(color) { // still here
    index += strlen(endtag);
  }
  char* output = new char[index+1];
  index = 0;
  output[index] = 0;
  color = False;
  for(k = 0;peptide[k]; k++) {
    if(! color && k == 0 && peptide[k] == 'n' && strchr(labeled_aas, '1') != NULL) {
      color = True;
      strcat(output, starttag);
      index = strlen(output);
    }
    else if(! color && peptide[k] == 'c' && strchr(labeled_aas, '2') != NULL) {
      color = True;
      strcat(output, starttag);
      index = strlen(output);
    }
    else if(! color && strchr(labeled_aas, peptide[k]) != NULL) { // turn on
      color = True;
      strcat(output, starttag);
      index = strlen(output);
    }
    else if(color && peptide[k] >= 'A' && peptide[k] <= 'Z' && strchr(labeled_aas, peptide[k]) == NULL) { // turn off
      color = False;
      strcat(output, endtag);
      index = strlen(output);
    }
    output[index++] = peptide[k];
    output[index] = 0;
  } // next aa of pep
  if(color) {
    strcat(output, endtag);
  }
  return output;
}
#ifdef USE_STD_MODS
  void  ProDataComponent::write(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times, char* radio, Array<char*>* misc_run_conds, char* colored_aas)
#endif
#ifndef USE_STD_MODS
  void  ProDataComponent::write(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times, char* radio, Array<char*>* aa_mods, Array<char*>* term_mods, Array<char*>* misc_run_conds)
#endif
{
#ifdef USE_STD_MODS
   write(os, inputfiles, basenames, pepproph_times, iproph_times, 
	 dbs, asap_times, NULL, NULL, NULL, NULL, radio,
	 misc_run_conds, colored_aas);
#endif
#ifndef USE_STD_MODS
   write(os, inputfiles, basenames, pepproph_times, iproph_times, 
	 dbs,  asap_times, NULL, NULL, NULL, NULL, radio,
	 aa_mods,  term_mods, misc_run_conds);
#endif


}


#ifdef USE_STD_MODS
void ProDataComponent::write(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times, Array<Boolean>* asap_quantHighBGs, Array<Boolean>* asap_zeroBGs, Array<double>* asap_mzBounds, Array<bool>* asap_wavelets, char* radio, Array<char*>* misc_run_conds, char* colored_aas) {
  Boolean color = colored_aas != NULL;
#endif
#ifndef USE_STD_MODS
void ProDataComponent::write(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times,  Array<Boolean>* asap_quantHighBGs, Array<Boolean>* asap_zeroBGs, Array<double>* asap_mzBounds,  Array<bool>* asap_wavelets, char* radio, Array<char*>* aa_mods, Array<char*>* term_mods, Array<char*>* misc_run_conds) {
#endif
  if(basename_index_ >= basenames->length()) {
    cout << "error with basenames" << endl;
    exit(1);
  }
  if(pepproph_timestamp_index_ >= pepproph_times->length()) {
    cout << "error with pepproph_times" << endl;
    exit(1);
  }
  if(iproph_timestamp_index_ >= iproph_times->length()) {
    cout << "error with iproph_times" << endl;
    exit(1);
  }
  if(asap_timestamp_index_ >= asap_times->length()) {
    cout << "error with asap_times" << endl;
    exit(1);
  }
  if(database_index_ >= dbs->length()) {
    cout << "error with databases" << endl;
    exit(1);
  }
#ifndef USE_STD_MODS
  if(database_index_ >= aa_mods->length()) {
    cout << "error with aa mods" << endl;
    exit(1);
  }
  if(database_index_ >= term_mods->length()) {
    cout << "error with term mods" << endl;
    exit(1);
  }
#endif


  char text[500];

  os << "<td>";
  if(result_->probability_ > -4.0) {
    if(result_->incomplete_prob_)
      sprintf(text, "%f", result_->probability_);
    else
      sprintf(text, "%0.4f%s", result_->probability_, result_->adjusted_prob_ ? "a" : "");

    string pepxml((*inputfiles)[xml_index_]);
    string waf = getDataUrl() + pepxml.substr(strlen(getDataPath()));
    string modelsFileNameWeb = waf.substr(0, waf.rfind('.')) + "-MODELS.html";

    os << "<a title='view analysis models' target='pepmodels' href='" << modelsFileNameWeb << "?Scores=" << score_summary_ << "&amp;Spectrum=" << result_->spectrum_ << "&amp;Prob=" << text << "'>";

    if(result_->adjusted_prob_)
      os << "<font color=\"#FF00FF\">";
    os << text;
    if(result_->adjusted_prob_)
      os << "</font>";
    os << "</a>" << endl;
  }
  os << "</td>" << endl;

  os << "<td>" << result_->spectrum_ << "</td>" << endl;


#ifndef USE_STD_MODS
  char* peptide = new char[strlen(result_->peptide_)+1];
  char* stripped = new char[strlen(result_->peptide_)+1];
  int start = 0;
  char prev = '-';
  char next = '-';

  int stop = strlen(result_->peptide_)-1;
  if(stop > start && result_->peptide_[1] == '.') {
    start = 2;
    prev = result_->peptide_[0];
  }
  if(stop > start && result_->peptide_[strlen(result_->peptide_)-2] == '.') {
    stop = strlen(result_->peptide_)-3;
    next = result_->peptide_[strlen(result_->peptide_)-1];
  }
  int index = 0;
  int stripped_index = 0;
  for(int k = start; k <= stop; k++) {
    peptide[index++] = result_->peptide_[k];
    if(result_->peptide_[k] >= 'A' && result_->peptide_[k] <= 'Z')
      stripped[stripped_index++] = result_->peptide_[k];
  }
  peptide[index] = 0;
  stripped[stripped_index] = 0;
#endif

#ifdef USE_STD_MODS
  // get the modification encoding to pass to spectrum view CGI program

  char prev = result_->prev_aa_;
  char next = result_->next_aa_;
  
  const int modquerylen = 1000;
  char modquery[modquerylen];
  modquery[0] = 0;

  const int mod_len = 500;
  char mod_encoding[mod_len];
  result_->setModificationEncoding(mod_encoding, mod_len);
  char mod_peptide[mod_len];
  if(result_->mod_info_ != NULL) {
    setModifiedPeptide(result_->mod_info_->getModifiedPeptide(), mod_peptide, mod_len); // modify for html outpu

    if(! result_->mod_info_->setQueryString(modquery, modquerylen)) {
      cout << "error: not enough length (" << modquerylen << ") for modification query string" << endl;
      exit(1);
    }
  }
  if(result_->mod_info_ == NULL || strlen(mod_peptide) == 0)
    strcpy(mod_peptide, result_->peptide_); // default


  char* colored = NULL;
  if(color)
    colored = getColoredPeptide(mod_peptide, colored_aas, "<font color=\"red\">", "</font>");


  os << "<td style='text-align:right'>" << endl;

  os << "<a title='view MS/MS spectrum' target=\"lorikeet\" href=\"" << getCgiUrl() << "plot-msms-js.cgi?" << modquery << "MassType=" << masstype_ << "&amp;NumAxis=1&amp;Pep=" << result_->peptide_ << "&amp;Dta=" << (*basenames)[basename_index_] << "/" << result_->spectrum_ << ".dta";
  if(! strcasecmp(result_->getName(), "Comet")) 
    os << "&amp;COMET=1";

  os << "\"><nobr>";
#endif
#ifndef USE_STD_MODS
  os << "<td style='text-align:right'>" << endl;

  if(! strcasecmp(result_->getName(), "Comet")) {
    if(database_index_ >= misc_run_conds->length()) {
      cout << "error with COMET md5 checksum" << endl;
      exit(1);
    }
    os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "cometplot.cgi?TarFile=" << (*basenames)[basename_index_] << ".cmt.tar.gz&amp;File=./" << result_->spectrum_ << ".dta&amp;Xmin=0&amp;Xmax=0&amp;Ymin=2&amp;Ymax=3&amp;LabelType=0&amp;NumAxis=1&amp;Pep=" << stripped << "&amp;ConfigFile=comet.def&amp;MD5=" << (*misc_run_conds)[basename_index_] << "&amp;PepMass=" << result_->neutral_mass_ << "&amp;ShowB=1&amp;ShowY=1&amp;AAMods=" << (*aa_mods)[basename_index_] << "&amp;TerminalMods=" << (*term_mods)[basename_index_] << "\"><nobr>";
  }
  else
    os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "sequest-tgz-plot.cgi?MassType=" << masstype_ << "&amp;NumAxis=1&amp;Pep=" << stripped << "&amp;Dta=" << (*basenames)[basename_index_] << "/" << result_->spectrum_ << ".dta\"><nobr>";

#endif

  if(result_->num_matched_ions_ < 100)
    os << " ";
  if(result_->num_matched_ions_ < 10)
    os << " ";
  
  os << result_->num_matched_ions_ << "/";

  if(result_->tot_num_ions_ < 100)
    os << " ";
  if(result_->tot_num_ions_ < 10)
    os << " ";
  os << result_->tot_num_ions_ << "</nobr></a>" << endl;
  os << "</td>" << endl;

  os << "<td>" << endl;

#ifdef USE_STD_MODS
  if(color)
    os << prev << ".<A TARGET=\"Win3\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=" << result_->peptide_ << "\">" << colored << "</A>." << next << endl;
  else
    os << prev << ".<A TARGET=\"Win3\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=" << result_->peptide_ << "\">" << mod_peptide << "</A>." << next << endl;

#endif
#ifndef USE_STD_MODS
  os << prev << ".<A TARGET=\"Win3\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=" << stripped << "\">" << peptide << "</A>." << next << endl;
#endif
  if(color && colored != NULL)
    delete colored;

  os << "</td>" << endl;

  os << "<td>" << endl;
#ifdef USE_STD_MODS
  os << "<a target=\"protseq\" href=\"" << getCgiUrl() << "madcaps.pl?Ref=" << result_->protein_ << "&amp;Db=" << (*dbs)[database_index_] << "&amp;Pep=" << result_->peptide_ << "&amp;MassType=" << masstype_ << "\">" << result_->protein_ << "</a>" << endl;
#endif
#ifndef USE_STD_MODS
  os << "<a target=\"protseq\" href=\"" << getCgiUrl() << "madcaps.pl?Ref=" << result_->protein_ << "&amp;Db=" << (*dbs)[database_index_] << "&amp;Pep=" << stripped << "&amp;MassType=" << masstype_ << "\">" << result_->protein_ << "</a>" << endl;
#endif

  if (result_->proteins_->size() > 1) {
    os << " <a title='view alignment of sequences in MADCAPS' target='madcaps' href='" << getCgiUrl() << "madcaps.pl?;Pep=" << result_->peptide_ <<  "&amp;Db=" << (*dbs)[database_index_] << "&amp;Ref=";
    for (int t = 0; t < result_->proteins_->size(); t++)
      os << (*result_->proteins_)[t] << " ";
    os << "'>+" << (result_->proteins_->size()-1) << " &cong;</a>" << endl;
  }

  os << "</td>" << endl;

  // XPRESS
  os << "<td class='value'>" << endl;
  if(strlen(xpress_) > 0.0) {
    size_t len;
    char *mzname = new char[len = (strlen((*basenames)[basename_index_])+12)];
    rampConstructInputFileName(mzname,(int)len,(*basenames)[basename_index_]);
    const char *mzExt = rampValidFileType(mzname);
    if (!mzExt) {
       mzExt = ".?";
    }
    os << "<a title='XPRESS quant' target=\"xpresspep\" href=\"" << getCgiUrl() << "XPressPeptideUpdateParser.cgi?LightFirstScan=" << LightFirstScan_ << "&amp;LightLastScan=" << LightLastScan_ << "&amp;HeavyFirstScan=" << HeavyFirstScan_ << "&amp;HeavyLastScan=" << HeavyLastScan_ << "&amp;XMLFile=" << (*basenames)[basename_index_] << mzExt << "&amp;ChargeState=" << result_->charge_ << "&amp;LightMass=" << LightMass_ << "&amp;HeavyMass=" << HeavyMass_ << "&amp;MassTol=" << MassTol_ << "&amp;index=" << data_index_ << "&amp;xmlfile=" << (*inputfiles)[xml_index_] << "&amp;bXpressLight1=" << bXpressLight_ << "&amp;OutFile=" << result_->spectrum_ << "\">";
    delete [] mzname;
    os << xpress_;
    //sprintf(text, "%0.2f", xpress_);
    //os << text << endl;
    //os << xpress_ << endl;
    os << "</a>" << endl;
  }
  else
    os << "n/a";
  os << "</td>\n";

  // ASAP
  // borrowed from PepXMLViewer:
  int intratio = int(50+33*log(0.00000001 + asap_mean_)); // log scale for -4.5 < log(r) < +4.5
  string qclass =
    (asap_mean_ < 0) ? "quantNC"  :
    (intratio   < 0) ? "quant0"   :
    (intratio > 100) ? "quant100" :
    "quant" + std::to_string(intratio);

  os << "<td class='value " << qclass << "'>" << endl;

  if(asap_mean_ > -3.0) {
    if(heavy2light_) {
      sprintf(text, "<nobr>%0.2f %s %0.2f</nobr>", asap_inv_mean_, "&plusmn;", asap_inv_error_);
      os << "<a title='inspect/curate ASAPRatio quantification' class='quant' target='asappep' href='" << getCgiUrl() << "ASAPRatioPeptideCGIDisplayParser.cgi?Xmlfile=" << (*inputfiles)[xml_index_] << "&amp;Basename=" << (*basenames)[basename_index_] << "&amp;Indx=" << asap_index_ << "&amp;Timestamp=" << (*asap_times)[asap_timestamp_index_] << "&amp;Spectrum=" << result_->spectrum_ << "&amp;ratioType=1";
    }
    else {
      sprintf(text, "<nobr>%0.2f %s %0.2f</nobr>", asap_mean_, "&plusmn;", asap_error_);
      os << "<a title='inspect/curate ASAPRatio quantification' class='quant' target='asapep' href='" << getCgiUrl() << "ASAPRatioPeptideCGIDisplayParser.cgi?Xmlfile=" << (*inputfiles)[xml_index_] << "&amp;Basename=" << (*basenames)[basename_index_] << "&amp;Indx=" << asap_index_ << "&amp;Timestamp=" << (*asap_times)[asap_timestamp_index_] << "&amp;Spectrum=" << result_->spectrum_ ;
    }

    if (asap_quantHighBGs != NULL && asap_timestamp_index_ < asap_quantHighBGs->length())
      os << "&amp;quantHighBG=" << (int)(*asap_quantHighBGs)[asap_timestamp_index_];

    if (asap_zeroBGs != NULL && asap_timestamp_index_ < asap_zeroBGs->length())
      os << "&amp;zeroBG=" << (int)(*asap_zeroBGs)[asap_timestamp_index_];

    if (asap_mzBounds != NULL && asap_timestamp_index_ < asap_mzBounds->length())
      os << "&amp;mzBound=" << (*asap_mzBounds)[asap_timestamp_index_];

    if (asap_wavelets != NULL && asap_timestamp_index_ < asap_wavelets->length())
      os << "&amp;wavelet=" << (int)(*asap_wavelets)[asap_timestamp_index_];

    os << "'>" << text << "</a>";

    double error;
    if(asap_mean_ == 0.0)
      error = 0.0;
    else
      error = asap_error_ * 100 / asap_mean_;
    if(error < 10.0)
      sprintf(text, "%0.1f%%", error);
    else
      sprintf(text, "%0.0f%%", error);

    os << " (" << text << ")</td>" << endl;

    //os << "<td>" << radio << endl;
  }

  os << "</td>" << endl;


#ifndef USE_STD_MODS
  if(peptide != NULL)
    delete peptide;
  if(stripped != NULL)
    delete stripped;
#endif
}


Boolean ProDataComponent::setModifiedPeptide(char* modpep, char* buffer, int buf_len) {
  if(buffer == NULL || ! buf_len || !modpep)
    return False;
  buffer[0] = 0;
  char subst1[] = "<font size=\"-2\">";
  char subst2[] = "</font>";
  size_t index = 0;
  for(int k = 0;modpep[k]; k++) 
    if(modpep[k] == '[') {
      if((int) (strlen(buffer) + strlen(subst1)) < buf_len)
	strcat(buffer, subst1);
      else
	return False;
      index = strlen(buffer);
    }
    else if(modpep[k] == ']') {
      if((int) (strlen(buffer) + strlen(subst2)) < buf_len)
	strcat(buffer, subst2);
      else
	return False;
      index = strlen(buffer);
    }
    else {
      if((int)(strlen(buffer) + 1) < buf_len) {
	buffer[index++] = modpep[k];
	buffer[index] = 0;
      }
      else
	return False;
    }
  buffer[index] = 0;
  return True;
}
#ifdef USE_STD_MODS
void ProDataComponent::writeStandardFormat(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times, Boolean include_header, Array<char*>* misc_run_conds, char* colored_aas) 
#endif
#ifndef USE_STD_MODS
void ProDataComponent::writeStandardFormat(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times, Boolean include_header, Array<char*>* aa_mods, Array<char*>* term_mods, Array<char*>* misc_run_conds) 
#endif
  {

#ifdef USE_STD_MODS
   writeStandardFormat(os, inputfiles, basenames, pepproph_times, iproph_times, dbs,
		       asap_times, NULL, NULL, NULL, NULL, include_header,  
		       misc_run_conds,  colored_aas);
#endif
#ifndef USE_STD_MODS
   writeStandardFormat(os, inputfiles, basenames, pepproph_times, iproph_times, dbs, 
		       asap_times, NULL, NULL, NULL, NULL, include_header, 
		       aa_mods, term_mods, misc_run_conds);
#endif

  }


#ifdef USE_STD_MODS
void ProDataComponent::writeStandardFormat(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times, Array<Boolean>* asap_quantHighBGs, Array<Boolean>* asap_zeroBGs, Array<double>* asap_mzBounds,  Array<bool>* asap_wavelets, Boolean include_header, Array<char*>* misc_run_conds, char* colored_aas) {
  Boolean color = colored_aas != NULL;
#endif
#ifndef USE_STD_MODS
void ProDataComponent::writeStandardFormat(ostream& os, Array<char*>* inputfiles, Array<char*>* basenames, Array<char*>* pepproph_times, Array<char*>* iproph_times, Array<char*>* dbs, Array<char*>* asap_times,  Array<Boolean>* asap_quantHighBGs, Array<Boolean>* asap_zeroBGs, Array<double>* asap_mzBounds,  Array<bool>* asap_wavelets, Boolean include_header, Array<char*>* aa_mods, Array<char*>* term_mods, Array<char*>* misc_run_conds) {
#endif

  if(basename_index_ >= basenames->length()) {
    cout << "error with basenames" << endl;
    exit(1);
  }
  if(pepproph_timestamp_index_ >= pepproph_times->length()) {
    cout << "error with pepproph_times: " << pepproph_timestamp_index_ << " vs " << pepproph_times->length() << endl;
    exit(1);
  }
  if(iproph_timestamp_index_ >= iproph_times->length()) {
    cout << "error with iproph_times: " << iproph_timestamp_index_ << " vs " << iproph_times->length() << endl;
    exit(1);
  }
  if(asap_timestamp_index_ >= asap_times->length()) {
    cout << "error with asap_times" << endl;
    exit(1);
  }
  if(database_index_ >= dbs->length()) {
    cout << "error with databases" << endl;
    exit(1);
  }

  os << "<table cellpadding=\"2\" bgcolor=\"white\">" << endl;

  char text[500];
  char font[] = ""; //<font face=\"Courier New\" size=\"-1\">";
  char endfont[] = ""; //</font>";

  if(include_header) {
    os << "<tr>";

    // pepproph
    os << "<td>" << endl;
    if((pepproph_timestamp_index_ >= 0 || iproph_timestamp_index_ >= 0) && result_->probability_ > -4.0) 
      os << "<font color=\"brown\"><b>probability</b></font>" << endl;
    os << "</td>" << endl;
    // spec
    os << "<td>" << endl;
    if(basename_index_ >= 0)
      os << "<font color=\"brown\"><b>spectrum</b></font>" << endl;
    os << "</td>" << endl;

    // scores
  os << "<td>" << endl;
  if(! strcasecmp(result_->getName(), "Sequest")) {
    os << "<table cellpadding=\"2\" bgcolor=\"white\"><TR><td width=\"40\" align=\"right\">";
    os << "<font color=\"brown\"><b>xcorr</b></font>" << endl;
    os << "</td><td width=\"40\" align=\"right\">";
    os << "<font color=\"brown\"><b>deltacn</b></font>" << endl;
    os << "</td><td width=\"40\" align=\"right\">" << endl;
    os << "<font color=\"brown\"><b>sprank</b></font>" << endl;
    os << "</td></TR></table>" << endl;
  } // if SEQUEST
  else if(! strcasecmp(result_->getName(), "Mascot")) {
    //MascotResult* mas = (MascotResult*)result_;
    os << "<table cellpadding=\"2\" bgcolor=\"white\"><TR><td width=\"40\" align=\"right\">";
    os << "<font color=\"brown\"><b>ionscore</b></font>" << endl;
    os << "</td><td width=\"40\" align=\"right\">" << endl;
    os << "<font color=\"brown\"><b>identity</b></font>" << endl;

    os << "</td><td width=\"40\" align=\"right\">";
    os << "<font color=\"brown\"><b>homology</b></font>" << endl;
    os << "</td></TR></table>" << endl;
  } // mascot
  else if(! strcasecmp(result_->getName(), "Comet")) {
    os << "<table cellpadding=\"2\" bgcolor=\"white\"><TR><td width=\"40\" align=\"right\">";
    os << "<font color=\"brown\"><b>dot product</b></font>" << endl;
    os << "</td><td width=\"40\" align=\"right\">" << endl;
    os << "<font color=\"brown\"><b>delta</b></font>" << endl;

    os << "</td><td width=\"40\" align=\"right\">";
    os << "<font color=\"brown\"><b>zscore</b></font>" << endl;
    os << "</td></TR></table>" << endl;
  } // comet
  os << "</td>" << endl;

  // matched ions
  os << "<td>" << endl;
  if(basename_index_ >= 0) 
    os << "<font color=\"brown\"><b>matched ions</b></font>" << endl;
  os << "</td>" << endl;
  // peptide
  os << "<td>" << endl;
  os << "<font color=\"brown\"><b>peptide</b></font>" << endl;
  os << "</td>" << endl;
  // protein
  os << "<td>" << endl;
  os << "<font color=\"brown\"><b>protein</b></font>" << endl;
  os << "</td>" << endl;
  // xpress
  os << "<td>" << endl;
  if(basename_index_ >= 0 && strlen(xpress_) > 0.0) 
    os << "<font color=\"brown\"><b>XPRESS</b></font>" << endl;
  os << "</td>" << endl;
  // asap
  os << "<td colspan=\"16\">" << endl;
  if(basename_index_ >= 0 && asap_timestamp_index_ >= 0 && asap_mean_ > -3.0) 
    os << "<font color=\"brown\"><b>ASAPRatio</b></font>" << endl;
  os << "</td>" << endl;

  os << "</tr>";
  }

  os << "<tr>";

  // pepproph
  os << "<td>" << endl;
  if((pepproph_timestamp_index_ >= 0 || iproph_timestamp_index_ >= 0 ) && result_->probability_ > -4.0) {
    if(result_->incomplete_prob_)
      sprintf(text, "%f", result_->probability_);
    else 
      sprintf(text, "%0.4f", result_->probability_);

    if (iproph_timestamp_index_ >= 0) {
      os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "ModelParser.cgi?Xmlfile=" << (*inputfiles)[xml_index_] << "&amp;Scores=" << score_summary_ << "&amp;Timestamp=" << (*iproph_times)[iproph_timestamp_index_] << "&amp;Spectrum=" << result_->spectrum_ << "&amp;Prob=" << text << "\">" << font;

    }
    else {
      os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "ModelParser.cgi?Xmlfile=" << (*inputfiles)[xml_index_] << "&amp;Scores=" << score_summary_ << "&amp;Timestamp=" << (*pepproph_times)[pepproph_timestamp_index_] << "&amp;Spectrum=" << result_->spectrum_ << "&amp;Prob=" << text << "\">" << font;
    }

    if(result_->adjusted_prob_)
      os << "<font color=\"#FF00FF\">";
    os << text;
    if(result_->adjusted_prob_)
      os << "</font>";
    os << endfont << "</A>" << endl;

  }
  os << "</td>" << endl;

  // spectrum
  os << "<td>" << endl;
  if(basename_index_ >= 0) {

    if(! strcasecmp(result_->getName(), "Sequest")) 
      os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "sequest-tgz-out.cgi?OutFile=" << (*basenames)[basename_index_] << "/" << result_->spectrum_ << ".out\">";
    else if(! strcasecmp(result_->getName(), "Mascot")) 
      os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "mascotout.pl?OutFile=" << (*basenames)[basename_index_] << "/" << result_->spectrum_ << ".out\">";
    else if(! strcasecmp(result_->getName(), "Comet")) 
      os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "cometresult.cgi?TarFile=" << (*basenames)[basename_index_] << ".cmt.tar.gz&amp;File=./" << result_->spectrum_ << ".cmt\"/>";
    else
      os << "<A TARGET=\"Win3\" HREF=\"\"/>";
    os << font << result_->spectrum_ << endfont << "</A>" << endl;
  }

  //os << "<A TARGET=\"Win3\" HREF=\"/cgi-bin/akeller/sequest-tgz-out.cgi?OutFile=" << (*basenames)[basename_index_] << "/" << result_->spectrum_ << ".out\">" << font << result_->spectrum_ << endfont << "</A>" << endl;
  os << "</td>" << endl;


  // scores
  os << "<td>" << endl;
  if(! strcasecmp(result_->getName(), "Sequest")) {
    SequestResult* seq = (SequestResult*)result_;
    os << "<table cellpadding=\"2\" bgcolor=\"white\"><TR><td width=\"40\" align=\"right\">";
    //os << "<table cellpadding=\"2\" bgcolor=\"white\" style=\"font-family: 'Courier New', Courier, mono; font-size: 10pt;\"><TR><td width=\"50\" align=\"right\">";
    sprintf(text, "%0.3f",seq->xcorr_);
    os << font << text << endfont << "</td><td width=\"40\" align=\"right\">";
    sprintf(text, "%0.3f",seq->delta_);
    os << font << text;

      //os << seq->xcorr_ << "</td><td width=\"50\" align=\"right\">" << seq->delta_;
    if(seq->deltastar_)
      os << "*";
    os << endfont << "</td><td width=\"40\" align=\"right\">" << font << seq->rank_ << endfont << "</td></TR></table>" << endl;
  } // if SEQUEST
  else if(! strcasecmp(result_->getName(), "Mascot")) {
    MascotResult* mas = (MascotResult*)result_;
    os << "<table cellpadding=\"2\" bgcolor=\"white\"><TR><td width=\"40\" align=\"right\">";
    //    os << "<table cellpadding=\"2\" bgcolor=\"white\" style=\"font-family: 'Courier New', Courier, mono; font-size: 10pt;\"><TR><td width=\"50\" align=\"right\">";
    sprintf(text, "%0.3f",mas->ionscore_);
    os << text;
    if(mas->star_)
      os << "*";
    sprintf(text, "%0.3f",mas->identity_);
    os << "</td><td width=\"40\" align=\"right\">" << text << "</td><td width=\"40\" align=\"right\">";

    sprintf(text, "%0.3f",mas->homology_);
    os << text;
    os << "</td></TR></table>" << endl;
  } // mascot
  else if(! strcasecmp(result_->getName(), "Comet")) {
    CometResult* com = (CometResult*)result_;
    os << "<table cellpadding=\"2\" bgcolor=\"white\"><TR><td width=\"40\" align=\"right\">";
    //    os << "<table cellpadding=\"2\" bgcolor=\"white\" style=\"font-family: 'Courier New', Courier, mono; font-size: 10pt;\"><TR><td width=\"50\" align=\"right\">";
    sprintf(text, "%0.3f",com->xcorr_);
    os << "<td width=\"40\" align=\"right\">" << text << "</td><td width=\"40\" align=\"right\">";

    sprintf(text, "%0.3f",com->delta_);
    if(com->deltastar_)
      os << "*";
    os << "</td><td width=\"40\" align=\"right\">" << text << "</td><td width=\"40\" align=\"right\">";

    sprintf(text, "%0.3f",com->expect_);
    os << text;
    os << "</td></TR></table>" << endl;
  } // comet

  os << "</td>" << endl;

#ifdef USE_STD_MODS
  char prev = result_->prev_aa_;
  char next = result_->next_aa_;
  // get the modification encoding to pass to spectrum view CGI program
  const int mod_len = 500;
  char mod_encoding[mod_len];
  result_->setModificationEncoding(mod_encoding, mod_len);
  char mod_peptide[mod_len];
  const int modquerylen = 1000;
  char modquery[modquerylen];

  if(result_->mod_info_ != NULL) {
    setModifiedPeptide(result_->mod_info_->getModifiedPeptide(), mod_peptide, mod_len); // modify for html outpu
    modquery[0] = 0;
    if(! result_->mod_info_->setQueryString(modquery, modquerylen)) {
      cout << "error: not enough length (" << modquerylen << ") for modification query string" << endl;
      exit(1);
    }
  }
  if(result_->mod_info_ == NULL || strlen(mod_peptide) == 0)
    strcpy(mod_peptide, result_->peptide_); // default

  char* colored = NULL;
  if(color)
    colored = getColoredPeptide(mod_peptide, colored_aas, "<font color=\"red\">", "</font>");

  os << "<td align=\"RIGHT\">" << endl;
  if(basename_index_ >= 0) {

    os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "plot-msms-js.cgi?" << modquery << "MassType=" << masstype_ << "&amp;NumAxis=1&amp;Pep=" << result_->peptide_ << "&amp;Dta=" << (*basenames)[basename_index_] << "/" << result_->spectrum_ << ".dta";
    if(! strcasecmp(result_->getName(), "Comet")) 
      os << "&amp;COMET=1";

    os << "\"><nobr>";


    /*
    if(! strcmp(result_->getName(), "Comet")) {
      if(database_index_ >= misc_run_conds->length()) {
	cout << "error with COMET md5 checksum" << endl;
	exit(1);
      }
      os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "cometplot.cgi?" << modquery << "TarFile=" << (*basenames)[basename_index_] << ".cmt.tar.gz&amp;File=./" << result_->spectrum_ << ".dta&amp;Xmin=0&amp;Xmax=0&amp;Ymin=2&amp;Ymax=3&amp;LabelType=0&amp;NumAxis=1&amp;Pep=" << result_->peptide_ << "&amp;ConfigFile=comet.def&amp;MD5=" << (*misc_run_conds)[basename_index_] << "&amp;PepMass=" << result_->neutral_mass_ << "&amp;ShowB=1&amp;ShowY=1&amp;PeptideMods=" << mod_encoding << "\"><nobr>";
    }
    else
      os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "sequest-tgz-plot.cgi?" << modquery << "MassType=" << masstype_ << "&amp;NumAxis=1&amp;Pep=" << result_->peptide_ << "&amp;Dta=" << (*basenames)[basename_index_] << "/" << result_->spectrum_ << ".dta\"><nobr>";
    */

#endif

#ifndef USE_STD_MODS

  char* peptide = new char[strlen(result_->peptide_)+1];
  char* stripped = new char[strlen(result_->peptide_)+1];
  int start = 0;
  char prev = '-';
  char next = '-';

  int stop = strlen(result_->peptide_)-1;
  if(stop > start && result_->peptide_[1] == '.') {
    start = 2;
    prev = result_->peptide_[0];
  }
  if(stop > start && result_->peptide_[strlen(result_->peptide_)-2] == '.') {
    stop = strlen(result_->peptide_)-3;
    next = result_->peptide_[strlen(result_->peptide_)-1];
  }
  int index = 0;
  int stripped_index = 0;
  for(int k = start; k <= stop; k++) {
    peptide[index++] = result_->peptide_[k];
    if(result_->peptide_[k] >= 'A' && result_->peptide_[k] <= 'Z')
      stripped[stripped_index++] = result_->peptide_[k];
  }
  peptide[index] = 0;
  stripped[stripped_index] = 0;


  // matched ions
  os << "<td align=\"RIGHT\">" << endl;
  if(basename_index_ >= 0) {
    if(! strcasecmp(result_->getName(), "Comet")) {
      if(database_index_ >= misc_run_conds->length()) {
	cout << "error with COMET md5 checksum" << endl;
	exit(1);
      }
      os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "cometplot.cgi?TarFile=" << (*basenames)[basename_index_] << ".cmt.tar.gz&amp;File=./" << result_->spectrum_ << ".dta&amp;Xmin=0&amp;Xmax=0&amp;Ymin=2&amp;Ymax=3&amp;LabelType=0&amp;NumAxis=1&amp;Pep=" << stripped << "&amp;ConfigFile=comet.def&amp;MD5=" << (*misc_run_conds)[basename_index_] << "&amp;PepMass=" << result_->neutral_mass_ << "&amp;ShowB=1&amp;ShowY=1&amp;AAMods=" << (*aa_mods)[basename_index_] << "&amp;TerminalMods=" << (*term_mods)[basename_index_] << "\"><nobr>";
    }
    else
      os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "sequest-tgz-plot.cgi?MassType=" << masstype_ << "&amp;NumAxis=1&amp;Pep=" << stripped << "&amp;Dta=" << (*basenames)[basename_index_] << "/" << result_->spectrum_ << ".dta\"><nobr>";
#endif



    //os << "<A TARGET=\"Win3\" HREF=\"/cgi-bin/akeller/sequest-tgz-plot.cgi?MassType=" << fragmasstype_ << "&amp;NumAxis=1&amp;Pep=" << stripped << "&amp;Dta=" << (*basenames)[basename_index_] << "/" << result_->spectrum_ << ".dta\"><nobr>";
    os << font;
    if(result_->num_matched_ions_ < 100)
      os << " ";
    if(result_->num_matched_ions_ < 10)
      os << " ";
  
    os << result_->num_matched_ions_ << "/";

    if(result_->tot_num_ions_ < 100)
      os << " ";
    if(result_->tot_num_ions_ < 10)
      os << " ";
    os << result_->tot_num_ions_ << endfont << "</nobr></A>" << endl;
  }
  os << "</td>" << endl;

  // peptide
  os << "<td>" << font << endl;

#ifdef USE_STD_MODS
  if(color)
    os << prev << ".<A TARGET=\"Win3\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=" << result_->peptide_ << "\">" << colored << "</A>." << next << endl;
  else
    os << prev << ".<A TARGET=\"Win3\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=" << result_->peptide_ << "\">" << mod_peptide << "</A>." << next << endl;


#endif
#ifndef USE_STD_MODS
  os << prev << ".<A TARGET=\"Win3\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=" << stripped << "\">" << peptide << "</A>." << next << endl;
#endif
  os << endfont << "</td>" << endl;

  if(color && colored != NULL)
    delete colored;

  // protein
  os << "<td>" << endl;

#ifdef USE_STD_MODS
  os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "madcaps.pl?Ref=" << result_->protein_ << "&amp;Db=" << (*dbs)[database_index_] << "&amp;Pep=" << result_->peptide_ << "&amp;MassType=" << masstype_ << "\">" << font << result_->protein_ << endfont << "</A>" << endl;
  os << "</td>" << endl;
#endif
#ifndef USE_STD_MODS
  os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "madcaps.pl?Ref=" << result_->protein_ << "&amp;Db=" << (*dbs)[database_index_] << "&amp;Pep=" << stripped << "&amp;MassType=" << masstype_ << "\">" << font << result_->protein_ << endfont << "</A>" << endl;
  os << "</td>" << endl;
#endif

  // xpress ratio
  os << "<td>" << endl;
  if(basename_index_ >= 0 && strlen(xpress_) > 0.0) {
    size_t len;
    char *mzname = new char[len=(strlen((*basenames)[basename_index_])+12)];
    rampConstructInputFileName(mzname,(int)len,(*basenames)[basename_index_]);
    const char *mzExt = rampValidFileType(mzname);
    if (!mzExt) {
       mzExt = ".?";
    }
    os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "XPressPeptideUpdateParser.cgi?LightFirstScan=" << LightFirstScan_ << "&amp;LightLastScan=" << LightLastScan_ << "&amp;HeavyFirstScan=" << HeavyFirstScan_ << "&amp;HeavyLastScan=" << HeavyLastScan_ << "&amp;XMLFile=" << (*basenames)[basename_index_] <<  mzExt << "&amp;ChargeState=" << result_->charge_ << "&amp;LightMass=" << LightMass_ << "&amp;HeavyMass=" << HeavyMass_ << "&amp;MassTol=" << MassTol_ << "&amp;index=" << data_index_ << "&amp;xmlfile=" << (*inputfiles)[xml_index_] << "&amp;bXpressLight1=" << bXpressLight_ << "&amp;OutFile=" << result_->spectrum_ << "\">";
    delete [] mzname;
    os << font << xpress_ << endfont << endl;
    //sprintf(text, "%0.2f", xpress_);
    //os << text << endl;
    //os << xpress_ << endl;
    os << "</A>" << endl;
  }
  os << "</td>";

  //asapratio
  os << "<td colspan=\"16\">" << endl;
  if(basename_index_ >= 0 && asap_timestamp_index_ >= 0 && asap_mean_ > -3.0) {

    //os << "here with asap" << endl;
    sprintf(text, "<nobr>%0.2f %s %0.2f</nobr>", asap_mean_, "&plusmn;", asap_error_);

    if(heavy2light_) {
      //os << "here!" << endl;
      os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "ASAPRatioPeptideCGIDisplayParser.cgi?Xmlfile=" << (*inputfiles)[xml_index_] << "&amp;Basename=" << (*basenames)[basename_index_] << "&amp;Indx=" << asap_index_ << "&amp;Timestamp=" << (*asap_times)[asap_timestamp_index_] << "&amp;Spectrum=" << result_->spectrum_ << "&amp;ratioType=1";
    }
    else {
      //os << "there" << endl;
      os << "<A TARGET=\"Win3\" HREF=\"" << getCgiUrl() << "ASAPRatioPeptideCGIDisplayParser.cgi?Xmlfile=" << (*inputfiles)[xml_index_] << "&amp;Basename=" << (*basenames)[basename_index_] << "&amp;Indx=" << asap_index_ << "&amp;Timestamp=" << (*asap_times)[asap_timestamp_index_] << "&amp;Spectrum=" << result_->spectrum_ ;
    }
    if (asap_quantHighBGs != NULL && asap_timestamp_index_ < asap_quantHighBGs->length()) {
      os << "&amp;quantHighBG=" << (int)(*asap_quantHighBGs)[asap_timestamp_index_];  
    }
    if (asap_wavelets != NULL && asap_timestamp_index_ < asap_wavelets->length()) {
      os << "&amp;wavelet=" << (int)(*asap_wavelets)[asap_timestamp_index_];  
    }
    if (asap_zeroBGs != NULL && asap_timestamp_index_ < asap_zeroBGs->length()) {
      os << "&amp;zeroBG=" << (int)(*asap_zeroBGs)[asap_timestamp_index_];
    }
    if (asap_mzBounds != NULL && asap_timestamp_index_ < asap_mzBounds->length()) {
      os << "&amp;mzBound=" << (*asap_mzBounds)[asap_timestamp_index_];
    }
    os << "\">" << font << text << endfont << "</A>";



    double error;
    if(asap_mean_ == 0.0)
      error = 0.0;
    else 
      error = asap_error_ * 100 / asap_mean_;
    if(error < 10.0)
      sprintf(text, "%0.1f%%", error);
    else
      sprintf(text, "%0.0f%%", error);
  }
  os << "</td>" << endl;

  os << "</tr>" << endl;

  os << "</table>" << endl;
#ifndef USE_STD_MODS
  if(peptide != NULL)
    delete peptide;
  if(stripped != NULL)
    delete stripped;
#endif
}
