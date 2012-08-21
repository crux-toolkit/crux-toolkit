
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
#include <cstdio>
#include <cstring>
#include <vector>
#include "Spectrum.h"
#include "DelimitedFile.h"
#include "SQTReader.h"
#include "SpectrumZState.h"
#include <fstream>

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
 
}

/**
 * Close the file, if open.
 */

void PinXMLWriter::closeFile(){
  if(output_file_!=NULL){
      
  }
}


void PinXMLWriter::write(
  MatchIterator* iterator,
  Spectrum* spectrum,
  int top_rank) {

  while (iterator -> hasNext()) {
    Match* match = iterator->next();
    if (match->getRank(XCORR) <= top_rank) {
      printPSM(match, spectrum, match->getNullPeptide());
    }
  }
}

void PinXMLWriter::write( 
  MatchCollection* target_collection,
  vector<MatchCollection*>& decoys,
  Spectrum* spectrum,
  int top_rank
 ){
  carp(CARP_DEBUG,"Start writing PIN xml file!");
  MatchIterator* target_match_iterator=new MatchIterator(target_collection);  
  MatchIterator* decoy_match_iterator=new MatchIterator(decoys[0]);  

  if (scan_number_ != spectrum->getFirstScan()) {
    if (scan_number_ > 0) {
      printFragSpecFooter();
    }
    scan_number_ = spectrum->getFirstScan();
    printHeaderFragSpectrumMatch(scan_number_);
  }
  
  write(target_match_iterator, spectrum, top_rank);
  write(decoy_match_iterator, spectrum, top_rank);

  delete target_match_iterator; 
  delete decoy_match_iterator;
}

/*creates a pin.xml file from two sqt, txt ,or pep.xml files*/
void PinXMLWriter::write( 
  MatchCollection* target_collection,
  vector<MatchCollection*>& decoys
 ){
  carp(CARP_DEBUG,"Start writing PIN xml file!");
  
  MatchIterator* target_match_iterator=new MatchIterator(target_collection);  
  MatchIterator* decoy_match_iterator=new MatchIterator(decoys[0]);  

  map<int, vector<Match*> > scan_to_matches;
  
 //set SP 
 if(target_collection->getScoredType(SP))
   is_sp_=true;
  
  while (target_match_iterator->hasNext()) {
 
    Match* match = target_match_iterator->next();
    //Spectrum* spectrum = match->getSpectrum();
    int scan = match->getSpectrum()->getFirstScan();

    if (scan_to_matches.find(scan) == scan_to_matches.end()) {
      vector<Match*> matches;
      scan_to_matches[scan] = matches;
    }
    scan_to_matches[scan].push_back(match);
    delete match;
  }
    
  while (decoy_match_iterator->hasNext()){
    
    Match* match = decoy_match_iterator->next();
    int scan = match->getSpectrum()->getFirstScan();
    
    if (scan_to_matches.find(scan) == scan_to_matches.end()) {
      vector<Match*> matches;
      scan_to_matches[scan] = matches;
    }
    scan_to_matches[scan].push_back(match); 
    delete match;
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
      if (matches[idx]->getRank(XCORR) <= 1){
        cerr<<"is Decoy? "<<matches[idx]->getNullPeptide()<<endl;
        printPSM(matches[idx],matches[idx]->getSpectrum(),matches[idx]->getNullPeptide());
      }
    }
  } 
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
  printProcessInfo("","","");
  printFeatureDescription();
}

void PinXMLWriter::printProcessInfo(
  string output_file_path, 
  string target_file_path,
  string decoy_file_path
){
   
  fprintf(
    output_file_, 
    "<process_info>"
    "  <command_line>/scratch/percolator-convertersBuild/sqt2pin -o "
    "</command_line></process_info>"
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
    fprintf(output_file_,"\n  <featureDescription name=\"lnrSP\"/>");
  }

  fprintf(output_file_,"\n  <featureDescription name=\"deltLCn\"/>");
  fprintf(output_file_,"\n  <featureDescription name=\"deltCn\"/>");
  fprintf(output_file_,"\n  <featureDescription name=\"Xcorr\"/>");
  if (is_sp_) {
    fprintf(output_file_,"\n  <featureDescription name=\"SP\"/>");
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
  fprintf(output_file_,"\n  <featureDescription name=\"LnNumSp\"/>");
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
  int is_decoy 
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
    
    //calculating experimental mz
   // double neutral_mass=zstate.getNeutralMass();
    double exp_mz= zstate.getMZ();
    
    if(is_decoy==1)
      decoy="true";
    else 
      decoy="false";
    string id = getId(charge,is_decoy,scan_number,1); 
    
    fprintf(
      output_file_,
      "\n  <peptideSpectrumMatch calculatedMassToCharge=\"%.*f\" "
      "chargeState=\"%i\""
      " experimentalMassToCharge=\"%.*f\""
      " id=\"%s\" isDecoy=\"%s\">\n",
      mass_precision_,zstate.getSinglyChargedMass(),
      charge, 
      mass_precision_,exp_mz,
      id.c_str(),
      decoy.c_str()
    );
   
    printFeatures(match);
    printPeptideSequence(peptide);
    printOccurence(flankC, flankN, peptide); 
    printPSMsFooter();
  

  delete []flanking_aas;
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
  double dM=(match->getNeutralMass())-(peptide->getPeptideMass());  
  if(peptide->getCTermFlankingAA()!='-')
     enz_c =true; 
  if(peptide->getNTermFlankingAA()!='-')
    enz_n =true;
 
  setCharges(charge_state, charge1, charge2, charge3);
  
  
  if(is_sp_){
    fprintf(//< print 17 features 
      output_file_, 
      "    <features>\n"
      "      <feature name=\"lnrSp\">%.*f</feature>\n"/*1.lnrSp*/
      "      <feature name=\"deltLCn\">%.*f</feature>\n"/*2.deltLCN*/
      "      <feature name=\"deltCn\">%.*f</feature>\n"/*3.deltCN*/
      "      <feature name=\"Xcorr\">%.*f</feature>\n"/*4.Xcorr*/
      "      <feature name=\"SP\">%.*f</feature>\n"/*5.SP*/
      "      <feature name=\"IonFrac\">%.*f</feature>\n"/*6.IonFrac */
      "      <feature name=\"Mass\">%.*f</feature>\n"/*7.Mass */
      "      <feature name=\"PepLen\">%u</feature>\n"/*8.PepLen */
      "      <feature name=\"Charge1\">%u</feature>\n"/*9.Chrge1 */
      "      <feature name=\"Charge2\">%u</feature>\n"/*10.charge2*/
      "      <feature name=\"Charge3\">%u</feature>\n"/*11.charge3*/
      "      <feature name=\"enzN\">%u</feature>\n"/*12.enzN*/
      "      <feature name=\"enzC\">%u</feature>\n"/*13.enzC*/
      "      <feature name=\"enzInt\">%u</feature>\n"/*14.enzInt*/
      "      <feature name=\"LnNumSP\">%f</feature>\n"/*15.LnNumSP*/   
      "      <feature name=\"dM\">%.*f</feature>\n"/*16.dM*/
      "      <feature name=\"absdM\">%.*f</feature>\n"/*17.absdM*/
      "    </features>\n",
      precision_,log(match->getRank(SP)),/*1*/
      precision_,0.0,/*2*/ 
      precision_,match->getDeltaCn(),/*3*/
      precision_,match->getScore(XCORR),/*4*/
      precision_,match->getScore(SP),/*5*/
      precision_,match-> getBYIonFractionMatched(),/*6*/
      mass_precision_,match->getNeutralMass(),/*7*/
      peptide->getLength(),/*8*/
      charge1,/*9*/ 
      charge2,/*10*/ 
      charge3,/*11*/ 
      enz_n,/*12*/
      enz_c,/*13*/
      peptide->getMissedCleavageSites(), /*14. enzInt */
      match->getLnExperimentSize(),/*15*/
      mass_precision_,dM,/*16*/  
      mass_precision_,fabs(dM)/*17*/
    );
  }else{ ///< print 14 features when Sp score is not valid 
    fprintf(
      output_file_, 
      "    <features>\n"
      "      <feature name=\"deltLCn\">%.*f</feature>\n"/*1.deltLCn*/
      "      <feature name=\"deltCn\">%.*f</feature>\n"/*2.deltCn*/ 
      "      <feature name=\"Xcorr\">%.*f</feature>\n"/*3.Xcorr*/
      "      <feature name=\"Mass\">%.*f</feature>\n"/*4.Mass*/
      "      <feature name=\"PepLen\">%u</feature>\n"/*5.PepLen */
      "      <feature name=\"Charge1\" >%u</feature>\n"/*6.Charge1*/  
      "      <feature name=\"Charge2\">%u</feature>\n"/*7.charge2*/
      "      <feature name=\"Charge3\">%u</feature>\n"/*8.charge3*/
      "      <feature name=\"enzC\">%u</feature>\n"/*9.enzC */
      "      <feature name=\"enzN\">%u</feature>\n"/*10.enzN */
      "      <feature name=\"enzInt\">%u</feature>\n"/*11.enzInt*/
      "      <feature name=\"LnNumSp\">%f</feature>\n"/*12.LnNumSp */ 
      "      <feature name=\"dM\">%.*f</feature>\n"/*13.dM */
      "      <feature name=\"absdM\">%.*f</feature>\n"/*14.absdM */
      "     </features>\n",
      precision_,0.0, /*deltLCn*/
      precision_,match->getDeltaCn(),/*deltCn*/
      precision_,match->getScore(XCORR),/*Xcorr*/
      mass_precision_,match->getNeutralMass(),/*Mass*/
      peptide->getLength(),/*PepLen*/
      charge1,/*charge1*/
      charge2,/*charge2*/
      charge3,/*charge3*/
      enz_n,/*enzN*/
      enz_c,/*EnzC*/
      peptide->getMissedCleavageSites(), /*enzInt */
      match->getLnExperimentSize(),/*LnNumSp*/
      mass_precision_,dM,  
      mass_precision_,fabs(dM)
    );
  }
 
}

void PinXMLWriter:: printPeptideSequence(Peptide* peptide){
  
  fprintf(
    output_file_,
    "    <peptide>\n      <peptideSequence>%s</peptideSequence>\n"
    "    </peptide>\n",
    peptide->getSequence() 
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
  for(unsigned prot_idx=0;prot_idx<num_protein;prot_idx++){
    string protein_id=protein_ids[prot_idx];
    fprintf(
      output_file_,
      "    <occurence flankC=\"%c\" flankN=\"%c\" proteinId=\"%s\"/>",
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
  string decoy; 
  if(is_decoy)
    decoy="reverce";
  else 
   decoy="target";
  psm_id<<decoy<<"_"<<setfill('0')<<setw(6)<<scan_number<<"_";
  psm_id<<charge<<"_"<<rank;
  return psm_id.str();   
}

//set chargers for features 
void PinXMLWriter::setCharges(
  int charge, 
  bool ch1, 
  bool ch2, 
  bool ch3){
  if (charge==1)
    ch1=true;
  else if (charge==2)
    ch2=false; 
  else if (charge==3)
    ch3=true;
  else 
    carp(CARP_DEBUG, "invalid charge!"); 
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
 
  vector<MatchCollection*> decoys;
  MatchCollection *target_collection= MatchCollectionParser::create(target_path, database_path);
  MatchCollection *decoy_collection = MatchCollectionParser::create(decoy_path, database_path);
 
  
  target_collection->spectrumSort(PERCOLATOR_SCAN);
  decoy_collection->spectrumSort(PERCOLATOR_SCAN);
  
  
  PinXMLWriter* writer=new PinXMLWriter();

  decoys.push_back(decoy_collection);

  string filename="sqt.pin.xml";
  string output_dir="crux-output";
  
  writer->openFile(filename.c_str(),output_dir.c_str(),true);
  
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

