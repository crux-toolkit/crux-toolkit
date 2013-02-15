
/**
 * \file PinXMLWriter.cpp
 * \brief Writes search results in the TPP-compliant .pin.xml format.
 */
#include <iostream>
#include "PinXMLWriter.h"
#include "MatchCollection.h"
#include "crux-utils.h"
#include "parameter.h"
#include "MatchCollectionParser.h"
#include <cctype>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <vector>
#include "Spectrum.h"
#include "DelimitedFile.h"
#include "SQTReader.h"
#include "SpectrumZState.h"
#include <fstream>
#include <limits>

using namespace std;
using namespace Crux; 

void PinXMLWriter::init() {
  output_file_ = NULL;
  mass_precision_ = get_int_parameter("mass-precision");
  ENZYME_T enzyme= get_enzyme_type_parameter("enzyme");
  enzyme_= enzyme_type_to_string(enzyme);
  precision_ = get_int_parameter("precision");
  is_sp_=get_boolean_parameter("compute-sp");
  scan_number_=-1;
  decoy_prefix_=get_string_parameter("decoy-prefix");
  is_decoy_=false;
  decoy_file_path_="";
  target_file_path_=""; 
  output_file_path_="";
  directory_="";
}

PinXMLWriter::PinXMLWriter(){
  init();
}

PinXMLWriter::~PinXMLWriter(){ 
  closeFile(); 
}

/**
 * Open a file of the given name.  Replace an existing file if
 * overwrite is true, else exit if an existing file is found.
 */
void PinXMLWriter::openFile(const char* filename, const char* output_dir,bool overwrite){
  output_file_=create_file_in_path(filename,output_dir,overwrite);
  if(output_file_==NULL){
    carp(CARP_FATAL,"Cann't open %s file ",filename);
    return;
  }
}

/**
 * Close the file, if open.
 */

void PinXMLWriter::closeFile(){
  if(output_file_!=NULL){
    fclose(output_file_);
    output_file_ = NULL;
  }
}

//
void PinXMLWriter::write(
  MatchIterator* iterator,
  Spectrum* spectrum,
  int top_rank) {
  vector<Match*> matches; 
  
  while (iterator -> hasNext()) {
    matches.push_back(iterator->next());
  }

  calculateDeltaCN(matches);
  for(unsigned idx=0; idx<matches.size(); idx++){
   
    if(matches[idx]->getRank(XCORR)<=top_rank)
      printPSM(matches[idx],spectrum,matches[idx]->getNullPeptide());
    else 
     return;
  }

}


void PinXMLWriter::calculateDeltaCN(map<pair<int, int>, vector<Match*> >& scan_charge_to_matches) {
  map<pair<int, int>, vector<Match*> >::iterator iter;

  for (iter = scan_charge_to_matches.begin();
    iter != scan_charge_to_matches.end();
    ++iter) {
    calculateDeltaCN(iter->second);
  }
}

void PinXMLWriter::calculateDeltaCN(vector<Match*>& collection) {

  MatchCollection tmp_matches = MatchCollection();
  tmp_matches.setScoredType(XCORR, true);
  for (vector<Match*>::iterator iter = collection.begin();
    iter != collection.end();
    ++iter) {
    tmp_matches.addMatch(*iter);
  }
  tmp_matches.calculateDeltaCn();

}

void PinXMLWriter::calculateDeltaCN(
  MatchCollection* target_collection, 
  vector<MatchCollection*>& decoys
) {

 
  map<pair<int,int>, vector<Match*> > scan_charge_to_target_matches;
  pair<int,int> scan_charge;
  MatchIterator target_match_iterator(target_collection, XCORR, true);  
  while (target_match_iterator.hasNext()){
    Match* match = target_match_iterator.next();
    int scan = match->getSpectrum()->getFirstScan();
    int charge = match->getZState().getCharge();
    scan_charge.first=scan; 
    scan_charge.second=charge; 
    if (
      scan_charge_to_target_matches.find(scan_charge) 
        == scan_charge_to_target_matches.end()
    ) {
      vector<Match*> matches;
      scan_charge_to_target_matches[scan_charge] = matches;
    }
    scan_charge_to_target_matches[scan_charge].push_back(match);
  }
  calculateDeltaCN(scan_charge_to_target_matches);

  // Allow for multiple sets of decoys
  for(
    vector<MatchCollection*>::iterator decoy_set = decoys.begin(); 
    decoy_set != decoys.end(); 
    decoy_set++
  ) {
    map<pair<int,int>, vector<Match*> > scan_charge_to_decoy_matches;
    MatchIterator decoy_match_iterator(*decoy_set, XCORR, true);  
    while (decoy_match_iterator.hasNext()){
      Match* match = decoy_match_iterator.next();
      int scan = match->getSpectrum()->getFirstScan();  
      int charge= match->getZState().getCharge();
      scan_charge.first=scan;
      scan_charge.second=charge;
      if (scan_charge_to_decoy_matches.find(scan_charge) == scan_charge_to_decoy_matches.end()) {
        vector<Match*> matches;
        scan_charge_to_decoy_matches[scan_charge] = matches;
      }
      scan_charge_to_decoy_matches[scan_charge].push_back(match);
    }
    calculateDeltaCN(scan_charge_to_decoy_matches);
  }

}

void PinXMLWriter::write( 
  MatchCollection* target_collection,
  vector<MatchCollection*>& decoys,
  Spectrum* spectrum,
  int top_rank
 ){
  carp(CARP_DEBUG,"Start writing PIN xml file!");
 //set SP 
 if(target_collection->getScoredType(SP))
   is_sp_=true; 

  calculateDeltaCN(target_collection, decoys);

  MatchIterator target_match_iterator(target_collection, XCORR, true);  

  if (scan_number_ != spectrum->getFirstScan()) { 
    if (scan_number_ > 0) {
      printFragSpecFooter();
    }    
    scan_number_ = spectrum->getFirstScan();
    printHeaderFragSpectrumMatch(scan_number_);
  }
  
  write(&target_match_iterator, spectrum, top_rank);

  // Allow for multiple sets of decoys
  for(
    vector<MatchCollection*>::iterator decoy_set = decoys.begin();
    decoy_set != decoys.end(); 
    decoy_set++
  ) {
    MatchIterator decoy_match_iterator(*decoy_set, XCORR, true);  
    write(&decoy_match_iterator, spectrum, top_rank);
  }

}


/*creates a pin.xml file from two sqt, txt ,or pep.xml files*/
void PinXMLWriter::write( 
  MatchCollection* target_collection,
  vector<MatchCollection*>& decoys
 ){
  carp(CARP_DEBUG,"Start writing PIN xml file!");
  
  calculateDeltaCN(target_collection, decoys);
  

  map<int, vector<Match*> > scan_to_matches;
  
 //set SP 
 if(target_collection->getScoredType(SP))
   is_sp_=true;
  
  MatchIterator target_match_iterator = MatchIterator(target_collection);  
  while (target_match_iterator.hasNext()) {
 
    Match* match = target_match_iterator.next();
    //Spectrum* spectrum = match->getSpectrum();
    int scan = match->getSpectrum()->getFirstScan();

    if (scan_to_matches.find(scan) == scan_to_matches.end()) {
      vector<Match*> matches;
      scan_to_matches[scan] = matches;
    }
    scan_to_matches[scan].push_back(match);
    
  }
    
  // Allow for multiple sets of decoys
  for(
    vector<MatchCollection*>::iterator decoy_set = decoys.begin(); 
    decoy_set != decoys.end(); 
    decoy_set++
  ) {
    MatchIterator decoy_match_iterator = MatchIterator(*decoy_set);  
    while (decoy_match_iterator.hasNext()){
      
      Match* match = decoy_match_iterator.next();
      int scan = match->getSpectrum()->getFirstScan();
      
      if (scan_to_matches.find(scan) == scan_to_matches.end()) {
        vector<Match*> matches;
        scan_to_matches[scan] = matches;
      }
      scan_to_matches[scan].push_back(match); 
      
    }
  }

  //print header
  printHeader();
  for (map<int, vector<Match*> >::iterator iter = scan_to_matches.begin();
    iter != scan_to_matches.end();
    ++iter) {

    vector<Match*> matches = iter->second;
    
    if(scan_number_!=iter->first){
      if(scan_number_>0)
        printFragSpecFooter();

      scan_number_=iter->first;
      printHeaderFragSpectrumMatch(scan_number_);
     
    }
    for (unsigned idx = 0;idx < matches.size();idx++) {   
      bool is_decoy=isDecoy(matches[idx]);
      if (matches[idx]->getRank(XCORR) <= 1){
        printPSM(matches[idx],matches[idx]->getSpectrum(),is_decoy);
      }
    }
  } 
}

bool PinXMLWriter:: isDecoy(Match* match){
  string protein_id=match->getPeptide()->getProteinIds(); 
  if(protein_id.find(decoy_prefix_)!=string::npos)
    return true;
  else if(match->getNullPeptide())
    return true;
  return false; 
}

void PinXMLWriter::printHeader(
){
  fprintf(
    output_file_, 
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" 
    "<experiment xmlns=\"http://per-colator.com/percolator_in/12\""
    " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\""
    " xsi:schemaLocation=\"http://per-colator.com/percolator_in/12 "
    "https://github.com/percolator/percolator/raw/pin-1-2/src/xml/percolator_in.xsd\">\n\n"
    "<enzyme>%s</enzyme>\n",
    enzyme_.c_str()
  );
  
  printProcessInfo();
  printFeatureDescription();
}

void PinXMLWriter::printProcessInfo(){
   
  fprintf(
    output_file_, 
    "<process_info>\n"
    "  <command_line>/scratch/percolator-convertersBuild/sqt2pin -o %s %s %s"
    "</command_line>\n</process_info>",
    output_file_path_.c_str(), 
    target_file_path_.c_str(),
    decoy_file_path_.c_str()
  );  
  
}
void PinXMLWriter::printFooter(){

  //print FragScpctrumScan for last scan
  printFragSpecFooter();
  fprintf(
    output_file_,
    "</experiment>\n"
  );

}
void PinXMLWriter:: printFeatureDescription(){
  fprintf(
    output_file_, 
    "\n<featureDescriptions xmlns=\"http://per-colator.com/percolator_in/12\">"
  );
  
  if (is_sp_) {
    fprintf(output_file_,"\n  <featureDescription name=\"lnrSp\"/>");
  }

  fprintf(output_file_,"\n  <featureDescription name=\"deltLCn\"/>");
  fprintf(output_file_,"\n  <featureDescription name=\"deltCn\"/>");
  fprintf(output_file_,"\n  <featureDescription name=\"Xcorr\"/>");
  if (is_sp_) {
    fprintf(output_file_,"\n  <featureDescription name=\"Sp\"/>");
    fprintf(output_file_,"\n  <featureDescription name=\"IonFrac\"/>");
  }
  fprintf(output_file_,"\n  <featureDescription name=\"Mass\"/>");
  fprintf(output_file_,"\n  <featureDescription name=\"PepLen\"/>");
  fprintf(output_file_,"\n  <featureDescription name=\"Charge1\"/>");
  fprintf(output_file_,"\n  <featureDescription name=\"Charge2\"/>");
  fprintf(output_file_,"\n  <featureDescription name=\"Charge3\"/>");
  fprintf(output_file_,"\n  <featureDescription name=\"enzN\"/>");
  fprintf(output_file_,"\n  <featureDescription name=\"enzC\"/>");
  fprintf(output_file_,"\n  <featureDescription name=\"enzInt\"/>");
  fprintf(output_file_,"\n  <featureDescription name=\"lnNumSP\"/>");
  fprintf(output_file_,"\n  <featureDescription name=\"dM\"/>");
  fprintf(output_file_,"\n  <featureDescription name=\"absdM\"/>");
  fprintf(output_file_,"\n</featureDescriptions>\n");

}

void PinXMLWriter:: printHeaderFragSpectrumMatch(int scan_number){
  
  fprintf(
    output_file_, 
    "\n<fragSpectrumScan xmlns=\"http://per-colator.com/percolator_in/12\""
    " scanNumber=\"%u\">",
    scan_number
  );
}

void PinXMLWriter:: printPSM(
  Match* match,
  Spectrum* spectrum, 
  bool is_decoy 
){ 
   
  int scan_number= spectrum->getFirstScan();
  Peptide* peptide= match->getPeptide();
  string decoy; 

  //finding flanking_aas
  char *flanking_aas= peptide->getFlankingAAs();
  char flankN=flanking_aas[0];
  char flankC=flanking_aas[1];
  SpectrumZState& zstate = match->getZState();
  int charge=zstate.getCharge();
   
  //TODO - Figure out what is apppropiate here.  The xml asks for m/z, but we
  //are providing singly charge mass for experimental m/z and mass for 
  //calculated m/z.
  //calculating singly charged mass
  //FLOAT_T exp_mz= zstate.getMZ();
  FLOAT_T exp_mass = zstate.getSinglyChargedMass();
  //isDecoy(match);
  if(is_decoy)
    decoy="true";
  else 
    decoy="false";

  string id = getId(charge,is_decoy,scan_number,1); 
    
  //FLOAT_T calculated_mz = (peptide->calcMass(get_mass_type_parameter("isotopic-mass")) + (FLOAT_T)charge * MASS_PROTON) / (FLOAT_T)charge;

  FLOAT_T calculated_mass = peptide->calcMass(get_mass_type_parameter("isotopic-mass"));

  fprintf(
    output_file_,
    "\n  <peptideSpectrumMatch calculatedMassToCharge=\"%.*f\" "
    "chargeState=\"%i\""
    " experimentalMassToCharge=\"%.*f\""
    " id=\"%s\" isDecoy=\"%s\">\n",
    mass_precision_, calculated_mass,
    charge, 
    mass_precision_,exp_mass,
    id.c_str(),
    decoy.c_str()
  );
   
  printFeatures(match);
  printPeptideSequence(peptide);
  printOccurence(flankC, flankN, peptide); 
  printPSMsFooter();
  
  myfree(flanking_aas);
}

void PinXMLWriter::printFeatures(
  Match* match
){
  
  //Spectrum* spectrum= match->getSpectrum();
  int charge_state = match->getCharge();
  bool charge1, charge2, charge3;
  SpectrumZState zstate= match->getZState ();
  charge1=charge2=charge3=false;   
  bool enz_c=false; 
  bool enz_n=false;  
  Peptide* peptide=match->getPeptide(); 
  FLOAT_T obs_mass=match->getZState().getSinglyChargedMass();
  FLOAT_T calc_mass=peptide->calcMass(get_mass_type_parameter("isotopic-mass"));
  FLOAT_T dM=(obs_mass - calc_mass)/charge_state;
  
  char c_flank = peptide->getCTermFlankingAA();
  char n_flank = peptide->getNTermFlankingAA();
  char* sequence = peptide->getSequence();
  char sequence_first = sequence[0];
  char sequence_last = sequence[strlen(sequence) - 1];
  delete[] sequence;
  if (n_flank == '-' ||
      ((n_flank == 'K' || n_flank == 'R') && sequence_first != 'P'))
    enz_n =true;
  if (c_flank == '-' ||
      ((sequence_last == 'K' || sequence_last == 'R') && c_flank != 'P'))
    enz_c =true; 
 
  switch(charge_state){
    case 1:
      charge1=true;
      break; 
    case 2: 
      charge2=true; 
      break; 
    case 3:
      charge3=true; 
      break; 
  }
  
  FLOAT_T ln_num_sp=match->getLnExperimentSize();
  FLOAT_T lnrSp=0.0; 
  if(match->getRank(SP)>0)  
    lnrSp=log(match->getRank(SP));
  FLOAT_T delta_cn = match->getDeltaCn() ; 
  FLOAT_T delta_lcn= match->getDeltaLCn();
 
  
  if(fabs(delta_cn)==numeric_limits<FLOAT_T>::infinity()){
    delta_cn=0;
  }else if(fabs(delta_lcn)==numeric_limits<FLOAT_T>::infinity()){
    delta_lcn=0;     
  }
  
  if(is_sp_){
    fprintf(//< print 17 features 
      output_file_, 
      "    <features>\n"
      "      <feature>%.*f</feature>\n"/*1.lnrSp*/
      "      <feature>%.*f</feature>\n"/*2.deltLCN*/
      "      <feature>%.*f</feature>\n"/*3.deltCN*/
      "      <feature>%.*f</feature>\n"/*4.Xcorr*/
      "      <feature>%.*f</feature>\n"/*5.SP*/
      "      <feature>%.*f</feature>\n"/*6.IonFrac */
      "      <feature>%.*f</feature>\n"/*7.Mass */
      "      <feature>%u</feature>\n"/*8.PepLen */
      "      <feature>%u</feature>\n"/*9.Chrge1 */
      "      <feature>%u</feature>\n"/*10.charge2*/
      "      <feature>%u</feature>\n"/*11.charge3*/
      "      <feature>%u</feature>\n"/*12.enzN*/
      "      <feature>%u</feature>\n"/*13.enzC*/
      "      <feature>%u</feature>\n"/*14.enzInt*/
      "      <feature>%.*f</feature>\n"/*15.LnNumSP*/   
      "      <feature>%.*f</feature>\n"/*16.dM*/
      "      <feature>%.*f</feature>\n"/*17.absdM*/
      "    </features>\n",
      precision_,lnrSp,/*1*/
      precision_,delta_lcn,/*2. delta_l_Cn*/ 
      precision_,delta_cn,/*3. delta_cn*/
      precision_,match->getScore(XCORR),/*4*/
      precision_,match->getScore(SP),/*5*/
      precision_,match-> getBYIonFractionMatched(),/*6*/
      mass_precision_,obs_mass,/*7*/
      peptide->getLength(),/*8*/
      charge1,/*9*/ 
      charge2,/*10*/ 
      charge3,/*11*/ 
      enz_n,/*12*/
      enz_c,/*13*/
      peptide->getMissedCleavageSites(), /*14. enzInt */
      precision_,ln_num_sp,/*15*/
      mass_precision_,dM,/*16*/  
      mass_precision_,fabs(dM)/*17*/
    );
  }else{ ///< print 14 features when Sp score is not valid 
    fprintf(
      output_file_, 
      "    <features>\n"
      "      <feature>%.*f</feature>\n"/*1.deltLCn*/
      "      <feature>%.*f</feature>\n"/*2.deltCn*/ 
      "      <feature>%.*f</feature>\n"/*3.Xcorr*/
      "      <feature>%.*f</feature>\n"/*4.Mass*/
      "      <feature>%u</feature>\n"/*5.PepLen */
      "      <feature>%u</feature>\n"/*6.Charge1*/  
      "      <feature>%u</feature>\n"/*7.charge2*/
      "      <feature>%u</feature>\n"/*8.charge3*/
      "      <feature>%u</feature>\n"/*9.enzN */
      "      <feature>%u</feature>\n"/*10.enzC */
      "      <feature>%u</feature>\n"/*11.enzInt*/
      "      <feature>%.*f</feature>\n"/*12.LnNumSp */ 
      "      <feature>%.*f</feature>\n"/*13.dM */
      "      <feature>%.*f</feature>\n"/*14.absdM */
      "     </features>\n",
      precision_,match->getDeltaLCn(), /*1.deltLCn*/
      precision_,match->getDeltaCn(),/*2.deltCn*/
      precision_,match->getScore(XCORR),/*3.Xcorr*/
      mass_precision_,obs_mass,/*4.Mass*/
      peptide->getLength(),/*5.PepLen*/
      charge1,/*6.charge1*/
      charge2,/*7.charge2*/
      charge3,/*8.charge3*/
      enz_n,/*9.enzN*/
      enz_c,/*10.EnzC*/
      peptide->getMissedCleavageSites(), /*11.enzInt */
      precision_,ln_num_sp,/*12.LnNumSp*/
      mass_precision_,dM,/*13*/  
      mass_precision_,fabs(dM)/*14*/
    );
  }
 
}

void PinXMLWriter:: printPeptideSequence(Peptide* peptide){
  
  // Get the modifications and convert them to XML for output
  char* modified_sequence =
    peptide->getModifiedSequenceWithMasses(MOD_MASSES_SEPARATE);
  stringstream mod_output;
  int aa_read = 0;  // stores the number of AA read for location attribute
  string mass_delta = ""; // stores the mass shift of the current mod
  // Iterate over each character in the modified sequence
  for (int i = 0; i < strlen(modified_sequence); ++i) {
    char cur = modified_sequence[i];
    if (cur == ',' || cur == ']') {  // Found end of modification
      mod_output << "      <modification location=\"" << aa_read << "\""
                    " monoisotopicMassDelta=\"" << mass_delta << "\">\n"
                    "        <uniMod accession=\"0\"/>\n"
                    "      </modification>\n";
      mass_delta.clear();
    } else if (cur == '-' && mass_delta == "") {
      mass_delta = '-';
    } else if (isdigit(cur) ||
              (cur == '.' && mass_delta.find('.') == string::npos)) {
      mass_delta += cur;
    } else if (isalpha(cur)) {
      ++aa_read;
    }
  }
  delete modified_sequence;

  // Output peptide element with sequence and modifications
  fprintf(
    output_file_,
    "    <peptide>\n"
    "      <peptideSequence>%s</peptideSequence>\n"
    "%s"
    "    </peptide>\n",
    peptide->getSequence(),
    mod_output.str().c_str()
  );
}

//wite Occurance 
void PinXMLWriter::printOccurence(
  char flankC, 
  char flankN, 
  Peptide* peptide
){
  
  vector<string> protein_ids;
  vector<string> protein_description; 
  unsigned num_protein=peptide->getProteinInfo(
    protein_ids,
    protein_description
  );
  
   
  //string decoy_prefix="rand_";
 
  for(unsigned prot_idx=0;prot_idx<num_protein;prot_idx++){
    string protein_id=protein_ids[prot_idx];
    fprintf(
      output_file_,
      "    \n<occurence flankC=\"%c\" flankN=\"%c\" proteinId=\"%s\"/>",
      flankC,
      flankN,
      protein_id.c_str()
   );
  }
 
}


void PinXMLWriter::printPSMsFooter(){

  fprintf(
    output_file_,
    "\n  </peptideSpectrumMatch>\n"
  );
}

void PinXMLWriter::printFragSpecFooter(){
  fprintf(
    output_file_,
    "</fragSpectrumScan>\n\n"
  );
}

string PinXMLWriter::getId(
  int charge,
  bool is_decoy, 
  int scan_number, 
  int rank
  
){
  ostringstream psm_id; 
  string fileId; 
  if(is_decoy)
    fileId="decoy";
  else 
   fileId="target";
  psm_id<<fileId<<"_"<<scan_number<<"_";
  psm_id<<charge<<"_"<<rank;
  return psm_id.str();   
}

bool PinXMLWriter::setProcessInfo(
  const char* target_file, 
  const char* decoy_file
){
  if (decoy_file ==NULL || target_file==NULL)
   return false; 
  decoy_file_path_=absolutPath(decoy_file);
  target_file_path_=absolutPath(target_file);
  //output_file_path_=absolutPath(output_file);
  return true; 
}


string PinXMLWriter::absolutPath(const char* file){
 string file_path;
 
 #if DARWIN
   char path_buffer[PATH_MAX];
   file_path =  string(realpath(file, path_buffer));
  #else
    file_path =  string(realpath(file, NULL));
  #endif
   
 return file_path;
}
 //TODO:reomve this code after some time of debugging 

#ifdef MAIN

int main(int argc, char** argv) {
  //argv[1] - target path
  //argv[2] - decoy path
  //argv[3] - database index path

  initialize_parameters();
 
  char* target_path = argv[1];
  char* decoy_path = argv[2];
  char* database_path = argv[3];
  char* output_filename=argv[4];
   
  vector<MatchCollection*> decoys;
  MatchCollection *target_collection= MatchCollectionParser::create(target_path, database_path);
  MatchCollection *decoy_collection = MatchCollectionParser::create(decoy_path, database_path);
  PinXMLWriter* writer=new PinXMLWriter();

  decoys.push_back(decoy_collection);
  
  string output_dir="crux-output";
  if(output_filename==NULL)
    output_filename="crux.pin.xml";
  writer->openFile(output_filename,output_dir.c_str(),true);
  if(output_filename!=NULL && target_path!=NULL&& decoy_path!=NULL )
    writer->setProcessInfo(target_path, decoy_path);
  writer->write(target_collection, decoys);
  
  writer->printFooter();
  writer->closeFile();
  
  //delete writer;
  //delete decoy_collection;
  //delete target_collection; 
 
  
 return 0;
}


#endif


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

