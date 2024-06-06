/*
Copyright 2017-2023, Michael R. Hoopmann, Institute for Systems Biology
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include "Parsers/mzIMLTools/mzIMLTools.h"
#include "Parsers/NeoPepXMLParser/NeoPepXMLParser.h"
#include "Parsers/NeoProtXMLParser/NeoProtXMLParser.h"
#include <cstdio>
#include <iostream>

using namespace std;

enum eTPPFType{
  tppUnknown = 0,
  tppPepXML,
  tppProtXML,
  tppMzID
};

enum eTPPSoftware{
  swUnknown = 0,
  swComet,
  swIProphet,
  swPeptideProphet,
  swProteinProphet,
  swPTMProphet
};

typedef struct sParams{
  int version;
  string inFile;
  string outFile;
  string decoy;
  eTPPFType fileType;
  sParams(){
    version = 0;
    inFile.clear();
    outFile.clear();
    decoy = "DECOY";
    fileType = tppUnknown;
  }
} sParams;

typedef struct sDB{
  string name;
  string seq;
} sDB;

typedef struct sPE{
  string dbRef;
  int start;
  int stop;
  char pre;
  char post;
} sPE;

typedef struct sAltPep{
  string pep;
  string pRef;
  vector<sPE> pe;
} sAltPep;

//Add various elements to the mzID file
void addAnalysisSoftware(CMzIdentML& m, NeoPepXMLParser& p);
string addDBSequence(CMzIdentML& m, string protSeq, string protDesc, string dbRef);
void addPTMProphet(CMzIdentML& m, CSpectrumIdentificationItem& m_sii, CnpxSearchHit* sh, size_t ptmIndex);
void addPTMProphetToPeptide(CPeptide& pep, CnpxSearchHit* sHit, size_t ptmIndex);
void addPTMProphetToSIP(CSpectrumIdentificationProtocol* ptm_sip, CSpectrumIdentificationProtocol* m_sip);

//Create various mzID elements (but not necessarily add them yet).
CPeptide createPeptide(CnpxSearchHit* sHit, CSpectrumIdentificationProtocol* sip);

//sigh...
double guessMassDiff(double mass, char aa);

size_t checkAnalysis(CnpxSearchHit* sHit, eTPPSoftware sw);
bool compareDB(const sDB& a, const sDB& b);
bool convertPep(CMzIdentML& m, const char* in, sParams& params);
bool convertProt(CMzIdentML& m, const char* in, sParams& params);
bool findPepPos(string pep, string prot, vector<sDB>& db, int& start, int& end, string& alt);
//bool makePeptideReferences(CMzIdentML& m, CSpectrumIdentificationProtocol* sip, string& peptide_ref, string& peptide_ref2, string& cvValue, NeoPepXMLParser& p, size_t i, size_t r);
bool processCMD(int argc, char* argv[], sParams& par);
void processProtXMLMods(CPeptide& m_pep, CnprIndistinguishablePeptide* ip, string& sKey);
void processProtXMLMods(CPeptide& m_pep, CnprPeptide* ip, string& sKey);
bool readDB(string fName, vector<sDB>& db);
bool summary(const char* in);
void swapIL(vector<string>& v, string pep, size_t pos);
void usage();

int main(int argc, char* argv[]){
  CMzIdentML m;
  NeoPepXMLParser pep;
  NeoProtXMLParser prot;

  cout << "tpp2mzid v0.9.15 (August 22 2023), copyright Mike Hoopmann, Institute for Systems Biology." << endl;
  cout << "Built using mzIMLTools: " << m.getMzIMLToolsVersion() << endl;
  cout << "Built using NeoPepXMLParser: " << pep.versionNeo() << endl;
  cout << "Built using NeoProtXMLParser: " << prot.versionNeo() << endl;

  sParams params;
  if(!processCMD(argc, argv, params)){
    usage();
    return 1;
  }

  m.setVersion(params.version);
  if (params.fileType == tppPepXML){
    if (!convertPep(m, params.inFile.c_str(),params)) return -2;
    cout << "1 Exporting: " << params.outFile << endl;
    m.writeFile(params.outFile.c_str());
    cout << "1 Done" << endl;
  } else if (params.fileType == tppProtXML){
    if (!convertProt(m, params.inFile.c_str(),params)) return -2;
    cout << "2 Exporting: " << params.outFile << endl;
    m.writeFile(params.outFile.c_str());
    cout << "2 Done" << endl;
  } else if (params.fileType == tppMzID){
    cout << "3 Summary: " << endl;
    summary(params.inFile.c_str());
  }

  return 0;
}

void addAnalysisSoftware(CMzIdentML& m, NeoPepXMLParser& p){
  eTPPSoftware sw;

  for (size_t i = 0; i < p.msms_pipeline_analysis[0].analysis_summary.size(); i++){
    sw = swUnknown;
    CnpxAnalysisSummary* pAnalysis = &p.msms_pipeline_analysis[0].analysis_summary[i];
    string sAnal = pAnalysis->analysis;
    string sVer = pAnalysis->version;
    if (pAnalysis->analysis.compare("ptmprophet") == 0){ // special handling for PTMProphet
      sw = swPTMProphet;
      sAnal = "PTMProphet";
      CnpxPTMProphetSummary* ptmS = &pAnalysis->ptmprophet_summary[0];
      sVer = ptmS->version;
    }
    string analysisSoftware_ref = m.addAnalysisSoftware(sAnal, sVer);

    //Add software protocol for specific software analyses
    sCvParam cv;
    switch (sw){
    case swPTMProphet:
    {
      CnpxPTMProphetSummary* ptmS = &pAnalysis->ptmprophet_summary[0];
      CSpectrumIdentificationProtocol* m_sip = m.analysisProtocolCollection.addSpectrumIdentificationProtocol(analysisSoftware_ref);
      CAdditionalSearchParams sp;
      cv.cvRef = "PSI-MS"; cv.accession = "MS:1001494"; cv.name = "no threshold"; m_sip->threshold.cvParam.push_back(cv);
      cv.accession = "MS:1001083"; cv.name = "ms-ms search"; m_sip->searchType.cvParam = cv;
      cv.accession = "MS:1002491"; cv.name = "modification localization scoring"; sp.cvParam.push_back(cv);
      cv.accession = "MS:1001055"; cv.name = "modification parameters"; cv.value = ptmS->mod_string; sp.cvParam.push_back(cv);
      cv.accession = "MS:1001413"; cv.name = "search tolerance minus value"; cv.value = ptmS->frag_ppm_tol; sp.cvParam.push_back(cv);
      cv.accession = "MS:1001412"; cv.name = "search tolerance plus value"; sp.cvParam.push_back(cv);
      m_sip->additionalSearchParams.push_back(sp);
      break;
    }
    default:
      break;
    }

  }
}

string addDBSequence(CMzIdentML& m, string protSeq, string protDesc, string dbRef){
  CDBSequence m_dbs;
  sCvParam cv;
  m_dbs.accession = protSeq;
  m_dbs.searchDatabaseRef = dbRef;
  cv.cvRef = "PSI-MS";
  cv.accession = "MS:1001088";
  cv.name = "protein description";
  cv.value = protDesc;
  m_dbs.cvParam.push_back(cv);
  return m.sequenceCollection.addDBSequence(m_dbs);
}

void addPTMProphet(CMzIdentML& m, CSpectrumIdentificationItem& m_sii, CnpxSearchHit* sh, size_t ptmIndex){
  CPeptide* ptmPep = m.sequenceCollection.getPeptide(m_sii.peptideRef);
  sCvParam pcv;
  pcv.cvRef = "PSI-MS";
  pcv.accession = "MS:1003147";
  pcv.name = "PTMProphet probability";
  int mi = 1;

  //iterate over PTMProphet results
  for (size_t e = 0; e<sh->analysis_result[ptmIndex].ptmprophet_result.size(); e++){
    CnpxPTMProphetResult* ppr = &sh->analysis_result[ptmIndex].ptmprophet_result[e];
    string sites = ppr->ptm.substr(0, ppr->ptm.find(':'));
    string sMass = ppr->ptm.substr(ppr->ptm.find(':') + 1, ppr->ptm.size());
    double dMass = atof(sMass.c_str());

    //Add modification scores
    vector<sCvParam> vcv;
    for (size_t g = 0; g<ppr->parameter.size(); g++){
      CnpxParameter* par = &ppr->parameter[g];
      sCvParam scv;
      scv.cvRef = "PSI-MS";
      if (par->name.compare("mean_best_prob") == 0) {
        scv.accession = "MS:1003148";
        scv.name = "PTMProphet mean best probability";
      } else if (par->name.compare("norm_info_gain") == 0) {
        scv.accession = "MS:1003149";
        scv.name = "PTMProphet normalized information content";
      } else if (par->name.compare("localized_mods") == 0) {
        scv.accession = "MS:1003150";
        scv.name = "PTMProphet information content";
      }
      if (scv.accession.empty()) continue;
      scv.value = par->value;
      vcv.push_back(scv);
    }

    string mIndex;

    //iterate over the peptide mod sites to find matching ones
    for (size_t g = 0; g<ptmPep->modification.size(); g++){
      CModification* ptmM = &ptmPep->modification[g];
      if (sites.find(ptmM->residues) == string::npos) continue;

      //iterate over cvparams to see if we've got a modification index
      size_t h;
      for (h = 0; h<ptmM->cvParam.size(); h++){
        if (ptmM->cvParam[h].accession.compare("MS:1002504") == 0) break;
      }
      if (h == ptmM->cvParam.size()) { //need to add index
        sCvParam mcv;
        mcv.cvRef = "PSI-MS"; mcv.accession = "MS:1002504"; mcv.name = "modification index";
        mcv.value = to_string(mi++);
        ptmM->cvParam.push_back(mcv);
        //cout << "Bad PTMProphet modification: " << ptmPep->id << "\t" << ptmPep->peptideSequence.text << "\t" << m_sii.peptideRef << "\t" << ptmM->residues << endl;
        //exit(666);
      }
      mIndex = ptmM->cvParam[h].value;

      //Add our PTMProphet results to the SII
      for (h = 0; h<ppr->mod_amino_acid_probability.size(); h++){
        CnpxModAminoAcidProbability* maap = &ppr->mod_amino_acid_probability[h];
        if (maap->position == ptmM->location){
          //need to add all the other scores
          for (size_t i = 0; i<vcv.size(); i++){
            sCvParam scv = vcv[i];
            scv.value = mIndex + ":" + scv.value + ":" + to_string(maap->position) + ":true";
            m_sii.addCvParam(scv);
          }
          pcv.value = mIndex + ":" + to_string(maap->probability).substr(0, 6) + ":" + to_string(maap->position) + ":true";
          maap->position = -maap->position; //mark this position as dealt-with
          m_sii.addCvParam(pcv);
          break;
        }
      }
    }

    //iterate over any remaining ambiguous positions in PTMProphet. Add them with the final mIndex;
    for (size_t g = 0; g<ppr->mod_amino_acid_probability.size(); g++){
      CnpxModAminoAcidProbability* maap = &ppr->mod_amino_acid_probability[g];
      if (maap->position <0) continue;
      //need to add all the other scores
      for (size_t i = 0; i<vcv.size(); i++){
        sCvParam scv = vcv[i];
        scv.value = mIndex + ":" + scv.value + ":" + to_string(maap->position) + ":true";
        m_sii.addCvParam(scv);
      }
      pcv.value = mIndex + ":" + to_string(maap->probability).substr(0, 6) + ":" + to_string(maap->position) + ":true";
      m_sii.addCvParam(pcv);
    }
  }
}

void addPTMProphetToPeptide(CPeptide& pep, CnpxSearchHit* sHit, size_t ptmIndex){
  //iterate over PTMProphet results
  int mi = 1;
  for (size_t e = 0; e<sHit->analysis_result[ptmIndex].ptmprophet_result.size(); e++){
    CnpxPTMProphetResult* ppr = &sHit->analysis_result[ptmIndex].ptmprophet_result[e];
    string sites = ppr->ptm.substr(0, ppr->ptm.find(':'));
    string sMass = ppr->ptm.substr(ppr->ptm.find(':') + 1, ppr->ptm.size());
    double dMass = atof(sMass.c_str());

    //iterate over the peptide mod sites to find matching ones
    for (size_t g = 0; g<pep.modification.size(); g++){
      CModification* ptmM = &pep.modification[g];
      if (sites.find(ptmM->residues) == string::npos) continue;

      //iterate over cvparams to see if we've got a modification index
      size_t h;
      for (h = 0; h<ptmM->cvParam.size(); h++){
        if (ptmM->cvParam[h].accession.compare("MS:1002504") == 0) break;
      }
      if (h == ptmM->cvParam.size()) { //need to add index
        sCvParam mcv;
        mcv.cvRef = "PSI-MS"; mcv.accession = "MS:1002504"; mcv.name = "modification index";
        mcv.value = to_string(mi++);
        ptmM->cvParam.push_back(mcv);
      }
    }
  }
}

void addPTMProphetToSIP(CSpectrumIdentificationProtocol* ptm_sip, CSpectrumIdentificationProtocol* m_sip){
  //Get PTMProphet Search Parameters
  size_t a;
  string sPTM;
  for(a=0;a<ptm_sip->additionalSearchParams[0].cvParam.size();a++){
    sCvParam* cv = &ptm_sip->additionalSearchParams[0].cvParam[a];
    if (cv->accession.compare("MS:1001055")==0){
      sPTM=cv->value;
      break;
    }
  }
 
  //Parse the PTM string
  string aa;
  string mass;
  double dMass;
  bool bMass=false;
  for(a=0;a<sPTM.size();a++){
    if(sPTM[a]==','){ //process PTM
      dMass=atof(mass.c_str());
      for(size_t b=0;b<aa.size();b++){
        string sAA;
        sAA+=aa[b];
        if (aa[b] == 'n' || aa[b] == 'c'){
          m_sip->modificationParams[0].addSearchModification(false, dMass, sAA, true);
        } else {
          m_sip->modificationParams[0].addSearchModification(false,dMass,sAA,false);
        }
      }

      aa.clear();
      mass.clear();
      bMass=false;
      dMass=0;
      continue;
    }
    if(sPTM[a]==':'){
      bMass=true;
      continue;
    }
    if(bMass) mass+=sPTM[a];
    else aa+=sPTM[a];
  }
  if(aa.size()>0){ //process last one
    dMass = atof(mass.c_str());
    for (size_t b = 0; b<aa.size(); b++){
      string sAA;
      sAA += aa[b];
      if (aa[b] == 'n' || aa[b] == 'c'){
        m_sip->modificationParams[0].addSearchModification(false, dMass, sAA, true);
      } else {
        m_sip->modificationParams[0].addSearchModification(false, dMass, sAA, false);
      }
    }
  }

}

size_t checkAnalysis(CnpxSearchHit* sHit, eTPPSoftware sw){
  for (size_t a = 0; a<sHit->analysis_result.size(); a++){
    switch (sw){
    case swPeptideProphet:
      if (sHit->analysis_result[a].analysis.compare("peptideprophet") == 0) return a;
      break;
    case swIProphet:
      if (sHit->analysis_result[a].analysis.compare("interprophet") == 0) return a;
      break;
    case swPTMProphet:
      if (sHit->analysis_result[a].analysis.compare("ptmprophet") == 0) return a;
      break;
    default:
      return SIZE_MAX;
    }
  }
  return SIZE_MAX;
}

bool compareDB(const sDB& a, const sDB& b){
  return (a.name.compare(b.name)<0);
}

bool convertPep(CMzIdentML& m, const char* in, sParams& params){

  cout << "Reading: " << in << endl;

  NeoPepXMLParser p;
  if (!p.read(in)){
    cout << "Failed to read pepXML: " << in << endl;
    //strip path and try locally
    string fname = in;
    string str;
    str = fname.substr(fname.rfind('/') + 1, fname.size());
    if (!p.read(&str[0])) {
      str = fname.substr(fname.rfind('\\') + 1, fname.size());
      if (!p.read(&str[0])) {
        cout << "Failed to find pepXML in current working directory." << endl;
        return false;
      } else {
        cout << "Opened " << str << " from current working directory." << endl;
      }
    } else {
      cout << "Opened " << str << " from current working directory." << endl;
    }
  }

  cout << p.size() << " psms read." << endl;

  //Iterate through all psms
  string analysisSoftware_ref;
  vector<sDB> db;
  string lastDB;
  int start, end;
  string altPep;

  //populate software used from the pepxml - eventually push this to its own function block
  addAnalysisSoftware(m, p);

  char dbid[32];
  size_t a;
  for (a = 0; a<p.msms_pipeline_analysis[0].msms_run_summary.size(); a++){ //iterate over all run summaries

    cout << "Processing run " << a + 1 << " of " << p.msms_pipeline_analysis[0].msms_run_summary.size() << ": ";
    int iPercent;
    int iTmp;
    iPercent = 0;
    printf("%2d%%", iPercent);
    fflush(stdout);

    CnpxMSMSRunSummary* mrs = &p.msms_pipeline_analysis[0].msms_run_summary[a]; //equivalent to <SpectrumIdentification>?
    if (mrs->search_summary.size()>1){//what do to if there is more than one search summary? How is this possible?
      cerr << "WARNING: MULTIPLE search_summary" << endl;
    }
    CnpxSearchSummary* ss = &mrs->search_summary[0];
    analysisSoftware_ref = m.addAnalysisSoftware(ss->search_engine, ss->search_engine_version);

    CSpectraData* m_sd = m.dataCollection.inputs.addSpectraData(mrs->base_name + mrs->raw_data);

    CSearchDatabase* m_db = m.dataCollection.inputs.addSearchDatabase(ss->search_database[0].local_path);

    //Read in database for this search (needed to find peptide positions).
    if (ss->search_database[0].local_path.compare(lastDB) != 0){
      lastDB = ss->search_database[0].local_path;
      if (!readDB(lastDB, db)){
        cout << "Cannot read pepXML search database: " << lastDB << endl;
        exit(-75);
      }
      sort(db.begin(), db.end(), compareDB); //Sort database by name (for faster lookup)
    }

    CSpectrumIdentificationProtocol* m_sip = m.analysisProtocolCollection.addSpectrumIdentificationProtocol(analysisSoftware_ref);
    sCvParam cv2;
    cv2.accession = "MS:1001494";
    cv2.cvRef = "PSI-MS";
    cv2.name = "no threshold";
    m_sip->threshold.cvParam.push_back(cv2);
    cv2.accession = "MS:1001083";
    cv2.name = "ms-ms search";
    m_sip->searchType.cvParam = cv2;

    CModificationParams m_mp;
    size_t b;
    for (b = 0; b<ss->aminoacid_modification.size(); b++){
      bool bFixed = true;
      if (ss->aminoacid_modification[b].variable[0] == 'Y') bFixed = false;
      bool bProtTerm = false;
      if (!ss->aminoacid_modification[b].protein_terminus.empty() && ss->aminoacid_modification[b].protein_terminus[0] == 'Y') bProtTerm = true;
      m_mp.addSearchModification(bFixed, ss->aminoacid_modification[b].massdiff, ss->aminoacid_modification[b].aminoacid, bProtTerm);
    }
    for (b = 0; b<ss->terminal_modification.size(); b++){
      bool bFixed = true;
      if (ss->terminal_modification[b].variable[0] == 'Y') bFixed = false;
      bool bProtTerm = true;
      string term = "c";
      if (ss->terminal_modification[b].terminus[0] == 'N') term = "n";
      m_mp.addSearchModification(bFixed, ss->terminal_modification[b].massdiff, term, bProtTerm);
    }
    m_sip->modificationParams.push_back(m_mp);

    //Special case: if PTMProphet was used, add any additional modifications as defined in PTMProphet
    for(b=0;b<m.analysisProtocolCollection.spectrumIdentificationProtocol.size();b++){ //This block checks for PTMProphet use.
      CSpectrumIdentificationProtocol* t_sip = &m.analysisProtocolCollection.spectrumIdentificationProtocol[b];
      size_t c;
      for(c=0;c<m.analysisSoftwareList.analysisSoftware.size();c++){
        CAnalysisSoftware* t_as = &m.analysisSoftwareList.analysisSoftware[c];
        if(t_as->id.compare(t_sip->analysisSoftwareRef)==0 && t_as->name.compare("PTMProphet")==0) break;
      }
      if (c == m.analysisSoftwareList.analysisSoftware.size()) continue;
      break;
    }
    if (b<m.analysisProtocolCollection.spectrumIdentificationProtocol.size()) { //if in this block, then PTMProphet was used.
      addPTMProphetToSIP(&m.analysisProtocolCollection.spectrumIdentificationProtocol[b],m_sip);
    }

    CSpectrumIdentificationList* m_sil = NULL;
    CSpectrumIdentification* m_si = m.addSpectrumIdentification(m_sd->id, m_db->id, m_sip->id, m_sil);


    for (b = 0; b<mrs->spectrum_query.size(); b++){

      //Update progress meter
      iTmp = (int)((double)b / mrs->spectrum_query.size() * 100);
      if (iTmp>iPercent){
        iPercent = iTmp;
        printf("\b\b\b%2d%%", iPercent);
        fflush(stdout);
      }

      CnpxSpectrumQuery* sq = &mrs->spectrum_query[b]; //equivalent to <SpectrumIdentificationResult>

      CSpectrumIdentificationResult m_sir;
      if (!sq->spectrumNativeID.empty())m_sir.spectrumID = sq->spectrumNativeID;
      else m_sir.spectrumID = sq->spectrum; //this might break systems that rely on the nativeID. However nativeID is OPTIONAL in pepXML.
      m_sir.name = sq->spectrum;
      m_sir.spectraDataRef = m_sd->id;
      sprintf(dbid, "%s_%d", m_sd->id.c_str(), (int)m_sil->spectrumIdentificationResult.size());
      m_sir.id = dbid;

      //Add RT
      m_sir.addCvParam("MS:1000894", (double)sq->retention_time_sec);
      size_t c;
      for (c = 0; c<sq->search_result.size(); c++){
        CnpxSearchResult* sr = &sq->search_result[c];

        size_t d;
        for (d = 0; d<sr->search_hit.size(); d++){
          CnpxSearchHit* sh = &sr->search_hit[d];  //equivalent to <SpectrumIdentificationItem>
          CPeptide m_p = createPeptide(sh, m_sip); //everything but the id (which will be assigned automatically

          CSpectrumIdentificationItem m_sii;
          vector<sAltPep> vap;
          char dbid[32];
          sprintf(dbid, "%s_%d", m_sir.id.c_str(), (int)m_sir.spectrumIdentificationItem.size());
          m_sii.id = dbid;
          m_sii.calculatedMassToCharge = (sh->calc_neutral_pep_mass + sq->assumed_charge*1.007276466) / sq->assumed_charge;
          m_sii.chargeState = sq->assumed_charge;
          m_sii.experimentalMassToCharge = (sq->precursor_neutral_mass + sq->assumed_charge*1.007276466) / sq->assumed_charge;
          m_sii.rank = sh->hit_rank;
          m_sii.peptideRef = m.sequenceCollection.addPeptide(m_p);

          //Add protein information
          CPeptideEvidence m_pe;
          m_pe.peptideRef = m_sii.peptideRef;
          m_pe.dbSequenceRef = addDBSequence(m, sh->protein, sh->protein_descr, m_db->id);
          if (!findPepPos(sh->peptide, sh->protein, db, start, end, altPep)){
            cout << "Cannot find peptide: " << sh->peptide << " in FASTA DB for " << sh->protein << endl;
            exit(-76);
          }

          m_pe.start = start;
          m_pe.end = end;
          m_pe.pre = sh->peptide_prev_aa[0];
          m_pe.post = sh->peptide_next_aa[0]; //need a way to determine if this is a c-term peptide?
          if(sh->protein.find(params.decoy)==0) m_pe.isDecoy=true;
          m_sii.peptideEvidenceRef.push_back(m.sequenceCollection.addPeptideEvidence(m_pe));
          for (size_t e = 0; e<sh->alternative_protein.size(); e++){ //Additional proteins are problematic because they might be to a different, yet equivalant, peptide sequence, such as STUPID instead of STUPLD.
            m_pe.dbSequenceRef = addDBSequence(m, sh->alternative_protein[e].protein, sh->alternative_protein[e].protein_descr, m_db->id);
            if (!findPepPos(sh->peptide, sh->alternative_protein[e].protein, db, start, end, altPep)){
              cout << "Cannot find peptide: " << sh->peptide << " in FASTA DB for " << sh->alternative_protein[e].protein << endl;
              exit(-76);
            }
            if (!altPep.empty()) { //we have an alternative peptide sequence, must be recorded as separate SpectrumIdentificationItem
              size_t f;
              for (f = 0; f<vap.size(); f++){
                if (vap[f].pep.compare(altPep) == 0) break;
              }
              if (f == vap.size()){
                sAltPep ap;
                ap.pep = altPep;
                m_p.peptideSequence.text = altPep;
                ap.pRef = m.sequenceCollection.addPeptide(m_p);
                vap.push_back(ap);
              }
              sPE altPE;
              altPE.dbRef = m_pe.dbSequenceRef;
              altPE.start = start;
              altPE.stop = end;
              altPE.pre = sh->alternative_protein[e].peptide_prev_aa[0];
              altPE.post = sh->alternative_protein[e].peptide_next_aa[0];
              vap[f].pe.push_back(altPE);
            } else {
              m_pe.start = start;
              m_pe.end = end;
              //m_pe.pre=m_pe.post='?';
              m_pe.pre = sh->alternative_protein[e].peptide_prev_aa[0];
              m_pe.post = sh->alternative_protein[e].peptide_next_aa[0];
              m_sii.peptideEvidenceRef.push_back(m.sequenceCollection.addPeptideEvidence(m_pe));
            }
          }

          //add scores
          for (size_t e = 0; e<sh->search_score.size(); e++){
            m_sii.addPSMValue(ss->search_engine, sh->search_score[e].name, sh->search_score[e].value);
          }
          //for XL
          ///*for (j = 0; j<p[i].psms->at(r).peptide->xlScores->size(); j++){
          //  sii->addPSMValue(pSearch.alg, p.getScoreLabel(p[i].psms->at(r).peptide->xlScores->at(j).index), p[i].psms->at(r).peptide->xlScores->at(j).value,"pep");
          //}*/  

          //add probabilities
          for (size_t e = 0; e<sh->analysis_result.size(); e++){
            if (sh->analysis_result[e].analysis.compare("peptideprophet") == 0){
              m_sii.addPSMValue("PeptideProphet", "Probability", sh->analysis_result[e].peptide_prophet_result.probability);
            } else if (sh->analysis_result[e].analysis.compare("interprophet") == 0){
              m_sii.addPSMValue("iProphet", "Probability", sh->analysis_result[e].interprophet_result.probability);
            }
          }

          //check for PTMProphet
          size_t ptmIndex = checkAnalysis(sh, swPTMProphet);
          if (ptmIndex != SIZE_MAX) addPTMProphet(m, m_sii, sh, ptmIndex); //this mod was analyzed by PTMProphet

          m_sir.spectrumIdentificationItem.push_back(m_sii);
          
          //if there are alternative, yet equivalent, sequences (i.e. STUPID vs STUPLD), add new SpectrumIdentificationItems
          for (size_t e = 0; e<vap.size(); e++){
            m_sii.peptideRef = vap[e].pRef;
            sprintf(dbid, "%s_%d", m_sir.id.c_str(), (int)m_sir.spectrumIdentificationItem.size());
            m_sii.id = dbid;
            m_sii.peptideEvidenceRef.clear();
            for (size_t f = 0; f<vap[e].pe.size(); f++){
              m_pe.peptideRef = vap[e].pRef;
              m_pe.dbSequenceRef = vap[e].pe[f].dbRef;
              m_pe.start = vap[e].pe[f].start;
              m_pe.end = vap[e].pe[f].stop;
              m_pe.pre = vap[e].pe[f].pre;
              m_pe.post = vap[e].pe[f].post;
              m_sii.peptideEvidenceRef.push_back(m.sequenceCollection.addPeptideEvidence(m_pe));
            }
            m_sir.spectrumIdentificationItem.push_back(m_sii);
          }

        }

      }

      m_sil->spectrumIdentificationResult.push_back(m_sir);

    }

    printf("\b\b\b100%%");
    cout << endl;

    mrs->clear(); //maybe recover some memory for very large files?

  }

  return true;
}

bool convertProt(CMzIdentML& m, const char* in, sParams& params){

  cout << "Reading: " << in << endl;

  size_t i;

  //read protXML file
  NeoProtXMLParser p;
  if (!p.read(in)){
    cout << "Failed to read protXML: " << in << endl;
    return false;
  }

  //read pepXML source files
  vector<string> source_files;
  string str;
  char cstr[10000];
  char* tok;
  strcpy(cstr, p.protein_summary.protein_summary_header.source_files.c_str());
  tok = strtok(cstr, " \t");
  while (tok != NULL){
    str = tok;
    source_files.push_back(str);
    tok = strtok(NULL, " \t");
  }
  for (i = 0; i<source_files.size(); i++){
    if (!convertPep(m, source_files[i].c_str(),params)) return false;
  }

  //sorted shortcut table?
  m.dataCollection.analysisData.buildPeptideEvidenceTable();

  string analysisSoftware_ref;
  CProteinDetectionProtocol* m_pdp;
  sCvParam cv;

  //Add protein inference software
  CnprProgramDetails* pd = &p.protein_summary.protein_summary_header.program_details;
  analysisSoftware_ref = m.addAnalysisSoftware(pd->analysis, pd->version);

  //ProteinDetectionProtocol
  cv.accession = "MS:1001494";
  cv.cvRef = "PSI-MS";
  cv.name = "no threshold";
  m_pdp = m.analysisProtocolCollection.addProteinDetectionProtocol(analysisSoftware_ref);
  m_pdp->threshold.cvParam.push_back(cv);

  //Add ProteinDetection to the AnalysisCollection (This generates a pointer to use for the ProteinDetectionList)
  vector<string> vISI;
  for (i = 0; i<m.dataCollection.analysisData.spectrumIdentificationList.size(); i++) vISI.push_back(m.dataCollection.analysisData.spectrumIdentificationList[i].id);
  CProteinDetectionList* m_pdl = NULL;
  CProteinDetection* m_pd = m.addProteinDetection(vISI, m_pdp->id, m_pdl);

  //ProteinDetectionList
  cv.cvRef = "PSI-MS";
  cv.accession = "MS:1002404";
  cv.name = "count of identified proteins";
  sprintf(cstr, "%d", (int)p.protein_summary.protein_group.size());
  cv.value = cstr;
  m_pdl->cvParam.push_back(cv);

  //Iterate over all protein groups
  size_t a;
  for (a = 0; a<p.protein_summary.protein_group.size(); a++){
    CnprProteinGroup* pg = &p.protein_summary.protein_group[a];
    CProteinAmbiguityGroup m_pag;

    //set id
    sprintf(cstr, "PAG_%d", (int)m_pdl->proteinAmbiguityGroup.size());
    m_pag.id = cstr;

    //add scores
    m_pag.addParamValue(pd->analysis, "probability", pg->probability);
    m_pag.addParamValue(pd->analysis, "protein_group_passes_threshold", "true");

    //iterate over all proteins in the group
    size_t b;
    for (b = 0; b<pg->protein.size(); b++){
      CnprProtein* prot = &pg->protein[b];
      CDBSequence m_dbs = m.getDBSequenceByAcc(prot->protein_name);
      CProteinDetectionHypothesis* m_pdh = m_pag.addProteinDetectionHypothesis(m_pag.id, m_dbs.id);

      //protein scores
      m_pdh->addParamValue(pd->analysis, "probability", prot->probability);
      m_pdh->addParamValue(pd->analysis, "coverage", prot->percent_coverage);
      if (b == 0) m_pdh->addParamValue(pd->analysis, "leading_protein", "");
      if (b == 0) m_pdh->addParamValue(pd->analysis, "group_representative", "");
      if (b>0 && prot->probability == pg->protein[0].probability) m_pdh->addParamValue(pd->analysis, "leading_protein", "");
      else m_pdh->addParamValue(pd->analysis, "non_leading_protein", "");

      size_t c;
      for (c = 0; c<prot->peptide.size(); c++){
        CnprPeptide* pep = &prot->peptide[c];

        typedef struct sKeyRef{
          string key;
          vector<string> m_peRef;
        } sKeyRef;
        vector<sKeyRef> vKeyRef;

        CPeptide m_pep;  //Create a peptide object which we can use to lookup and match to existing peptides in the mzID.
        m_pep.peptideSequence.text = pep->peptide_sequence;

        size_t d;
        for (d = 0; d<pep->indistinguishable_peptide.size(); d++){ //iterate over all indistinguishable peptides																																	 
          //cout << "Next up: " << pep->indistinguishable_peptide[d].peptide_sequence << " " << pep->indistinguishable_peptide[d].charge << endl;
          m_pep.peptideSequence.text = pep->indistinguishable_peptide[d].peptide_sequence;  //update the peptide text, which might differ in cases of I/L.
          int charge = pep->indistinguishable_peptide[d].charge;
          string sKey = pep->indistinguishable_peptide[d].peptide_sequence; //this is the modified peptide sequence as a string for tracking vKeyRef identifiers.
          processProtXMLMods(m_pep, &pep->indistinguishable_peptide[d], sKey); //process any modifications

          size_t e;
          for (e = 0; e<vKeyRef.size(); e++){ //key is used to identify if peptide was already seen (for example at a different charge, which is processed elsewhere).
            if (sKey.compare(vKeyRef[e].key) == 0) break;
          }
          if (e == vKeyRef.size()){ //never seen this peptide variant before
            sKeyRef kr;
            bool bRet;
            kr.key = sKey;
            bRet = m.sequenceCollection.getPeptideEvidenceFromPeptideAndProtein(m_pep, m_dbs.id,kr.m_peRef);
            if (!bRet){
			  //cout << "Swap necessary" << endl;								 
              //Maybe I/L need to be swapped (STUPID==STUPLD)
              vector<string> vs;
              swapIL(vs, m_pep.peptideSequence.text, 0);
              size_t f;
              for (f = 0; f<vs.size(); f++){
                m_pep.peptideSequence.text = vs[f];
                bRet = m.sequenceCollection.getPeptideEvidenceFromPeptideAndProtein(m_pep, m_dbs.id, kr.m_peRef);
                if (bRet) break;
              }
              if (f == vs.size()){
                cout << "Warning (A), no peptide evidence found for (mod): " << pep->peptide_sequence << " " << str << " from " << prot->protein_name << endl;
                for (size_t tr = 0; tr<m_pep.modification.size(); tr++) cout << m_pep.modification[tr].monoisotopicMassDelta << endl;
                exit(1);
              }
            }
            vKeyRef.push_back(kr);
          }

          size_t pa;
          vector<string> vSII;
          for(pa=0;pa<vKeyRef[e].m_peRef.size();pa++){
            m.dataCollection.analysisData.getSpectrumIdentificationItems(vKeyRef[e].m_peRef[pa], charge, vSII);
            if(vSII.size()>0) break;
          }
          if (pa == vKeyRef[e].m_peRef.size()){
            cout << "Big fail (A) on :" << pep->peptide_sequence << " (" << c << "," << d << ")\t" << pg->group_number << "\t" << pg->probability << "\t" << prot->probability << "\t" << prot->percent_coverage << "\t" << prot->protein_name << "(" << b << ")" << "\t" << m_dbs.id << endl;
            cout << "Looking for: " << pep->indistinguishable_peptide[d].peptide_sequence<< "\t" << sKey << "," << charge << "\t" << m_pep.peptideSequence.text << endl;
            cout << "KeyRef: " << vKeyRef.size() << "\t" << e << endl;
            cout << "m_pep: " << m_pep.peptideSequence.text << "\t";	  
            for(size_t x=0;x<m_pep.modification.size(); x++) {
              cout << m_pep.modification[x].location << "," << m_pep.modification[x].monoisotopicMassDelta << "  ";
            }
            cout << endl;
            /*for(size_t y=0;y<vKeyRef.size();y++){
              cout << vKeyRef[y].key << " - " << vKeyRef[y].m_peRef.size() << endl;
              for (pa = 0; pa<vKeyRef[y].m_peRef.size(); pa++){
                m.dataCollection.analysisData.getSpectrumIdentificationItems(vKeyRef[e].m_peRef[pa], charge, vSII);
                cout << vKeyRef[e].m_peRef[pa] << endl;
                for (size_t x = 0; x<vSII.size(); x++) cout << " " << x << "\t" << vSII[x] << endl;
              }
            }*/
            exit(666);
          }
          CPeptideHypothesis* m_ph = m_pdh->addPeptideHypothesis(vKeyRef[e].m_peRef[pa]);
          size_t f;
          for (f = 0; f < vSII.size(); f++){
            sSpectrumIdentificationItemRef s;
            s.text = vSII[f];
            m_ph->spectrumIdentificationItemRef.push_back(s);
          }

        } //d

        //for case where there are no indistinguishable peptides
        if (pep->indistinguishable_peptide.size() == 0){
          int charge = pep->charge;
          string sKey;
          processProtXMLMods(m_pep, pep, sKey); //process any modifications

          //redundant code from here down...perhaps find a better way to deal with this.
          vector<string> peptideEvidence_ref;
          bool bRet;
          bRet = m.sequenceCollection.getPeptideEvidenceFromPeptideAndProtein(m_pep, m_dbs.id,peptideEvidence_ref);
          if (!bRet){
            cout << "Warning (B), no peptide evidence found for (mod): " << m_pep.peptideSequence.text << " " << str << " from " << prot->protein_name << endl;
            exit(1);
          }
          
          size_t pa;
          vector<string> vSII;
          for(pa=0;pa<peptideEvidence_ref.size();pa++){
            m.dataCollection.analysisData.getSpectrumIdentificationItems(peptideEvidence_ref[pa], charge, vSII);
            if(vSII.size()>0) break;
          }
          if (pa == peptideEvidence_ref.size()) {
            cout << "Big fail (B) on " << peptideEvidence_ref[0] << "\t" << sKey << "\t" << m_pep.peptideSequence.text << endl;
            for (size_t x = 0; x<m_pep.modification.size(); x++){
              cout << m_pep.modification[x].location << "\t" << m_pep.modification[x].monoisotopicMassDelta << endl;
            }
            exit(666);
          }
          
          CPeptideHypothesis* m_ph = m_pdh->addPeptideHypothesis(peptideEvidence_ref[pa]);
          size_t f;
          for (f = 0; f < vSII.size(); f++){
            sSpectrumIdentificationItemRef s;
            s.text = vSII[f];
            m_ph->spectrumIdentificationItemRef.push_back(s);
          }
        } // done for all peptides

      }

    }

    m_pdl->proteinAmbiguityGroup.push_back(m_pag);
  }

  return true;
}

CPeptide createPeptide(CnpxSearchHit* sHit, CSpectrumIdentificationProtocol* sip){
  CPeptide m_p; //everything but the id (which will be assigned automatically)
  sCvParam cv;
  m_p.peptideSequence.text = sHit->peptide;
  size_t e;
  int mi = 1;
  if (sHit->modification_info.size()>0){
    if (sHit->modification_info[0].mod_nterm_mass != 0){ //not handling ptm prophet here!!
      CModification m_m;
      m_m.location = 0;
      m_m.monoisotopicMassDelta = sHit->modification_info[0].mod_nterm_mass - 1.007825;
      cv = sip->modificationParams.back().getModificationCvParam(m_m.monoisotopicMassDelta, string("."), true, false);
      if (!cv.accession.empty()) m_m.cvParam.push_back(cv);
      m_p.modification.push_back(m_m);
    }
    if (sHit->modification_info[0].mod_cterm_mass != 0){ //not handling ptm prophet here!!
      CModification m_m;
      m_m.location = (int)m_p.peptideSequence.text.size() + 1;
      m_m.monoisotopicMassDelta = sHit->modification_info[0].mod_cterm_mass;
      cv = sip->modificationParams.back().getModificationCvParam(m_m.monoisotopicMassDelta, string("."), false, true);
      if (!cv.accession.empty()) m_m.cvParam.push_back(cv);
      m_p.modification.push_back(m_m);
    }
    for (e = 0; e<sHit->modification_info[0].mod_aminoacid_mass.size(); e++){
      CnpxModAminoAcidMass* mm = &sHit->modification_info[0].mod_aminoacid_mass[e];
      CModification m_m;
      m_m.location = mm->position;
      m_m.residues = sHit->peptide[mm->position - 1];
      m_m.monoisotopicMassDelta = 0;
      if (mm->staticMass != 0) {
        m_m.monoisotopicMassDelta += mm->staticMass;
        cv = sip->modificationParams.back().getModificationCvParam(mm->staticMass, m_m.residues, false, false);
        if (!cv.accession.empty()) m_m.cvParam.push_back(cv);
      }
      if (mm->variable != 0) { //can have both static and variable mods listed together on the same amino acid.
        m_m.monoisotopicMassDelta += mm->variable;
        cv = sip->modificationParams.back().getModificationCvParam(mm->variable, m_m.residues, false, false);
        if (!cv.accession.empty()) m_m.cvParam.push_back(cv);
      } else if (mm->variable == 0 && mm->staticMass == 0){ //neither variable or static is identified. Usually a result of David.
        m_m.monoisotopicMassDelta = guessMassDiff(mm->mass, m_m.residues[0]);
        cv = sip->modificationParams.back().getModificationCvParam(m_m.monoisotopicMassDelta, m_m.residues, false, false);
        if (!cv.accession.empty()) m_m.cvParam.push_back(cv);
      }
      m_p.modification.push_back(m_m);
    }
  }

  //check for PTMProphet
  size_t ptmIndex = checkAnalysis(sHit, swPTMProphet);
  if (ptmIndex != SIZE_MAX) addPTMProphetToPeptide(m_p, sHit, ptmIndex); //this mod was analyzed by PTMProphet

  return m_p;
}

bool findPepPos(string pep, string prot, vector<sDB>& db, int& start, int& end, string& alt){
  alt.clear();
  if (db.size() == 0) return false;

  size_t sz = db.size();
  size_t lower = 0;
  size_t mid = sz / 2;
  size_t upper = sz;
  int i;

  i = db[mid].name.compare(prot);
  while (i != 0){
    if (lower >= upper) return false;
    if (i>0){
      if (mid == 0) return false;
      upper = mid - 1;
      mid = (lower + upper) / 2;
    } else {
      lower = mid + 1;
      mid = (lower + upper) / 2;
    }
    if (mid == sz) return false;
    i = db[mid].name.compare(prot);
  }

  size_t pos;
  pos = db[mid].seq.find(pep);
  if (pos == string::npos) { //Maybe sequence has different I/L amino acids.
    vector<string> vs;
    swapIL(vs, pep, 0);
    size_t a;
    for (a = 0; a<vs.size(); a++){
      pos = db[mid].seq.find(vs[a]);
      if (pos != string::npos) break;
    }
    if (a == vs.size()) return false;
    alt = vs[a];
  }
  start = (int)pos + 1;
  end = start + (int)pep.size();
  return true;

}

double guessMassDiff(double mass, char aa){
  switch (aa){
  case 'G': return mass - 57.02146;
  case 'A': return mass - 71.03711;
  case 'S': return mass - 87.03203;
  case 'P': return mass - 97.05276;
  case 'V': return mass - 99.06841;
  case 'T': return mass - 101.04768;
  case 'C': return mass - 103.00918;
  case 'L': return mass - 113.08406;
  case 'I': return mass - 113.08406;
  case 'N': return mass - 114.04293;
  case 'D': return mass - 115.02694;
  case 'Q': return mass - 128.05858;
  case 'K': return mass - 128.09496;
  case 'E': return mass - 129.04259;
  case 'M': return mass - 131.04048;
  case 'O': return mass - 132.08988;
  case 'H': return mass - 137.05891;
  case 'F': return mass - 147.06841;
  case 'U': return mass - 150.95363;
  case 'R': return mass - 156.10111;
  case 'Y': return mass - 163.06333;
  case 'W': return mass - 186.07931;
  default: return mass;
  }
}

bool processCMD(int argc, char* argv[], sParams& par){
  size_t pos;
  int i;

  for (i = 1; i<argc; i++){

    //check parameters
    if (argv[i][0] == '-'){

      if(argv[i][1] == 'd'){ //set decoy identifier
        i++;
        if(i==argc) {
          cout << "Bad -d parameter." << endl;
          return false;
        }
        par.decoy=argv[i];
      } else if (argv[i][1] == 'v'){ //check version number
        if (par.version != 0){
          cout << "Duplicated parameter: " << argv[i] << endl;
          exit(-1);
        }
        switch (argv[i][2]){
        case '1': par.version = 1; break;
        case '2': par.version = 2; break;
        default:
          cout << "Unknown version number: " << argv[i] << endl;
          return false;
        }
      } else {
        cout << "Unknown parameter: " << argv[i] << endl;
        return false;
      }
      continue;
    }

    //check files
    if (par.inFile.size()>0) {
      cout << "Multiple file specified: " << argv[i] << endl;
      cout << "Please specify only a single file." << endl;
      exit(-1);
    }
    par.inFile = argv[i];

    pos = par.inFile.find(".pep.xml");
    if (pos != string::npos)par.fileType = tppPepXML;
    if (par.fileType == tppUnknown){
      pos = par.inFile.find(".prot.xml");
      if (pos != string::npos)par.fileType = tppProtXML;
    }
    if (par.fileType == tppUnknown){
      pos = par.inFile.find(".prot.xml");
      if (pos != string::npos)par.fileType = tppProtXML;
    }
    if (par.fileType == tppUnknown){
      pos = par.inFile.find(".mzid");
      if (pos != string::npos)par.fileType = tppMzID;
    }
    if (par.fileType == tppUnknown){
      cout << "Uknown or unsupported file type: " << argv[i] << endl;
      exit(-1);
    }

    par.outFile = par.inFile;
    switch (par.fileType){
    case tppPepXML: par.outFile.replace(pos + 4, 4, ".mzid"); break;
    case tppProtXML: par.outFile.replace(pos + 5, 4, ".mzid"); break;
    case tppMzID:
    default:
      break;
    }
  }

  if(par.inFile.empty()) return false;

  //always default the version to the latest if not specified.
  if (par.version == 0) par.version = 2;
  return true;
}

void processProtXMLMods(CPeptide& m_pep, CnprIndistinguishablePeptide* ip, string& sKey){
  size_t e;
  m_pep.modification.clear(); //clear any mods from the last pass
  string stripped_pep; //generate the stripped peptide
  for (e = 0; e<ip->modification_info.size(); e++){ //iterate over all modified forms; Is this always 0 or 1?
    CModification m_mod;
    string str = ip->modification_info[e].modified_peptide;
    sKey = str; //update key
    int z = 0;
    bool bMass = false;
    char aa;
    string sMass;
    for (size_t n = 0; n < str.size(); n++){
      if (str[n] == 'n') {
        z++;
        m_mod.location = 0;
      } else if (str[n] == '['){
        z++;
        sMass.clear();
        aa = str[n - 1];
        bMass = true;
      } else if (str[n] == ']'){
        z++;
        m_mod.monoisotopicMassDelta = atoi(sMass.c_str());
        switch (aa){
        case 'n': m_mod.monoisotopicMassDelta -= 1; break;
        case 'A': m_mod.monoisotopicMassDelta -= 71; break;
        case 'C': m_mod.monoisotopicMassDelta -= 103; break;
        case 'E': m_mod.monoisotopicMassDelta -= 129; break;
        case 'G': m_mod.monoisotopicMassDelta -= 57; break;
        case 'H': m_mod.monoisotopicMassDelta -= 137; break;
        case 'I': m_mod.monoisotopicMassDelta -= 113; break;
        case 'K': m_mod.monoisotopicMassDelta -= 128; break;
        case 'L': m_mod.monoisotopicMassDelta -= 113; break;
        case 'M': m_mod.monoisotopicMassDelta -= 131; break;
        case 'N': m_mod.monoisotopicMassDelta -= 114; break;
        case 'P': m_mod.monoisotopicMassDelta -= 97; break;
        case 'Q': m_mod.monoisotopicMassDelta -= 128; break;
        case 'S': m_mod.monoisotopicMassDelta -= 87; break;
        case 'T': m_mod.monoisotopicMassDelta -= 101; break;
        case 'V': m_mod.monoisotopicMassDelta -= 99; break;
        case 'W': m_mod.monoisotopicMassDelta -= 186; break;
        case 'Y': m_mod.monoisotopicMassDelta -= 163; break;
        default: break;
        }
        if (m_mod.location == -1) m_mod.location = (int)n - (int)z + 1;
        m_pep.modification.push_back(m_mod);
        m_mod.location = -1;
        bMass = false;
      } else {
        if (bMass) {
          sMass += str[n];
          z++;
        } else stripped_pep+=str[n];
      }
    } //n
	//update the m_pep structure to the stripped peptide;
    m_pep.peptideSequence.text=stripped_pep;													 
    //if additional modifications are listed as separate mod_mainoacid_mass elements, process them here
    for (size_t a = 0; a<ip->modification_info[e].mod_aminoacid_mass.size(); a++){
      int loc = ip->modification_info[e].mod_aminoacid_mass[a].position;
      size_t b;
      for (b = 0; b<m_pep.modification.size(); b++){
        if (m_pep.modification[b].location == loc) break;
      }
      if (b == m_pep.modification.size()){ //must add this modification
        int mMass = (int)(ip->modification_info[e].mod_aminoacid_mass[a].mass + 0.5);
        //correct amino acid mass
        switch (m_pep.peptideSequence.text[loc-1]) { //ip->peptide_sequence[loc - 1]){ //lookup peptide from unmodified sequence
        case 'C': mMass -= 103; break;
        default:
          cout << "Tell Mike to fix tpp2mzid::processProtXMLMods()" << endl;
		  cout << str << ": " << ip->peptide_sequence[loc - 1] << " " << mMass << " " << m_pep.peptideSequence.text << "\t" << ip->peptide_sequence << endl;																																				
          exit(1);
          break;
        }
        m_mod.location = loc;
        m_mod.monoisotopicMassDelta = mMass;
        m_pep.modification.push_back(m_mod);
      }
    }
  } //e

}

//This might need updating to reflect changes similar to the function above.																			
void processProtXMLMods(CPeptide& m_pep, CnprPeptide* ip, string& sKey){
  size_t e;
  m_pep.modification.clear(); //clear any mods from the last pass
  for (e = 0; e<ip->modification_info.size(); e++){ //iterate over all modified forms
    CModification m_mod;
    string str = ip->modification_info[e].modified_peptide;
    sKey = str; //update key
    int z = 0;
    bool bMass = false;
    char aa;
    string sMass;
    for (size_t n = 0; n < str.size(); n++){
      if (str[n] == 'n') {
        z++;
        m_mod.location = 0;
      } else if (str[n] == '['){
        z++;
        sMass.clear();
        aa = str[n - 1];
        bMass = true;
      } else if (str[n] == ']'){
        z++;
        m_mod.monoisotopicMassDelta = atoi(&sMass[0]);
        switch (aa){
        case 'n': m_mod.monoisotopicMassDelta -= 1; break;
        case 'A': m_mod.monoisotopicMassDelta -= 71; break;
        case 'C': m_mod.monoisotopicMassDelta -= 103; break;
        case 'E': m_mod.monoisotopicMassDelta -= 129; break;
        case 'G': m_mod.monoisotopicMassDelta -= 57; break;
        case 'H': m_mod.monoisotopicMassDelta -= 137; break;
        case 'I': m_mod.monoisotopicMassDelta -= 113; break;
        case 'K': m_mod.monoisotopicMassDelta -= 128; break;
        case 'L': m_mod.monoisotopicMassDelta -= 113; break;
        case 'M': m_mod.monoisotopicMassDelta -= 131; break;
        case 'N': m_mod.monoisotopicMassDelta -= 114; break;
        case 'P': m_mod.monoisotopicMassDelta -= 97; break;
        case 'Q': m_mod.monoisotopicMassDelta -= 128; break;
        case 'S': m_mod.monoisotopicMassDelta -= 87; break;
        case 'T': m_mod.monoisotopicMassDelta -= 101; break;
        case 'V': m_mod.monoisotopicMassDelta -= 99; break;
        case 'W': m_mod.monoisotopicMassDelta -= 186; break;
        case 'Y': m_mod.monoisotopicMassDelta -= 163; break;
        default: break;
        }
        if (m_mod.location == -1) m_mod.location = (int)n - (int)z + 1;
        m_pep.modification.push_back(m_mod);
        m_mod.location = -1;
        bMass = false;
      } else {
        if (bMass) {
          sMass += str[n];
          z++;
        }
      }
    } //n

    //if additional modifications are listed as separate mod_mainoacid_mass elements, process them here
    for (size_t a = 0; a<ip->modification_info[e].mod_aminoacid_mass.size(); a++){
      int loc = ip->modification_info[e].mod_aminoacid_mass[a].position;
      size_t b;
      for (b = 0; b<m_pep.modification.size(); b++){
        if (m_pep.modification[b].location == loc) break;
      }
      if (b == m_pep.modification.size()){ //must add this modification
        int mMass = (int)(ip->modification_info[e].mod_aminoacid_mass[a].mass + 0.5);
        //correct amino acid mass
        switch (ip->peptide_sequence[loc - 1]){
        case 'C': mMass -= 103; break;
        default:
          cout << "Tell Mike to fix tpp2mzid::processProtXMLMods()" << endl;
          exit(1);
          break;
        }
        m_mod.location = loc;
        m_mod.monoisotopicMassDelta = mMass;
        m_pep.modification.push_back(m_mod);
      }
    }
  } //e

}

bool readDB(string fName, vector<sDB>& db){
  char str[65000];
  char* tok;
  FILE* f;
  sDB d;
  char c;

  d.name.clear();

  db.clear();
  f = fopen(fName.c_str(), "rt");
  if (f == NULL) {
    string grr = fName.substr(fName.rfind('/') + 1, fName.size());
    f = fopen(grr.c_str(), "rt");
    if (f == NULL) return false;
  }

  while (!feof(f)){
    if (fgets(str, 65000, f) == NULL) continue;
    if (strlen(str)>0){
      tok = strtok(str, "\r\n"); //strip terminators
      if (tok == NULL) continue;
      strcpy(str, tok);
    }
    if (str[0] == '>') {
      if (d.name.size() != 0) {
        if (d.seq.length()>65000){
          cout << "  WARNING: " << &d.name[0] << " has a sequence that is too long. It will be skipped." << endl;
        } else {
          db.push_back(d);
        }
      }
      tok = strtok(str, "> \t");
      d.name = tok;
      d.seq = "";
    } else {
      for (size_t i = 0; i<strlen(str); i++){
        c = toupper(str[i]);
        if (c == ' ' || c == '\t') continue;
        d.seq += c;
      }
    }
  }
  fclose(f);
  if (d.seq.length()>65000){
    cout << "  WARNING: " << &d.name[0] << " has a sequence that is too long. It will be skipped." << endl;
  } else {
    db.push_back(d);
  }
  return true;

}

bool summary(const char* in){
  CMzIdentML m;
  CPSM psm;
  size_t i;
  int j;

  cout << "Reading: " << in << endl;

  if (!m.readFile(in)) return false;

  cout << in << endl;
  cout << "Analyzed with: " << endl;
  for (i = 0; i < m.analysisSoftwareList.analysisSoftware.size(); i++){
    cout << " " << m.analysisSoftwareList.analysisSoftware[i].name << " version: " << m.analysisSoftwareList.analysisSoftware[i].version << endl;
  }
  cout << "Search Databases: " << endl;
  for (i = 0; i < m.dataCollection.inputs.searchDatabase.size(); i++){
    cout << " " << m.dataCollection.inputs.searchDatabase[i].location << endl;
  }
  cout << "Total Unique Proteins in DB: " << m.sequenceCollection.dbSequence.size() << endl;
  if (m.dataCollection.analysisData.proteinDetectionList.size()>0){
    cout << "Total Protein Groups:        " << m.dataCollection.analysisData.proteinDetectionList[0].proteinAmbiguityGroup.size() << endl;
  }
  cout << "Total Unique Peptides:       " << m.sequenceCollection.peptide.size() << endl;
  cout << "Total PSMs:                  " << m.getPSMCount() << endl;

  cout << "First PSM that maps to at least 20 proteins: " << endl;
  for (i = 0; i < m.getPSMCount(); i++){
    psm = m.getPSM((int)i);
    if (psm.proteinCount<20) continue;
    cout << " " << psm.scanInfo.scanID
      << " " << psm.scanInfo.rTimeSec
      << " " << psm.sequenceMod
      << " " << psm.getScore(0).value << endl;
    for (j = 0; j < psm.proteinCount; j++){
      cout << "  " << psm.getProtein(j) << endl;
    }
    break;
  }
  cout << i << " psms were explored to reach this point." << endl;

  cout << "Exporting hopethisworks.mzid" << endl;
  m.writeFile("hopethisworks.mzid");

  return true;
}

void swapIL(vector<string>& v, string pep, size_t pos){
  size_t i;
  if (pos>0) v.push_back(pep);
  for (i = pos; i<pep.size(); i++){
    if (pep[i] == 'I') {
      string s = pep;
      s[i] = 'L';
      swapIL(v, s, i + 1);
    } else if (pep[i] == 'L'){
      string s = pep;
      s[i] = 'I';
      swapIL(v, s, i + 1);
    }
  }
}

void usage(){
  cout << "This program converts pepXML and protXML to mzIdentML (mzid)." << endl;
  cout << "When given a pepXML or protXML file, it will convert it to mzIdentML." << endl;
  cout << "When given an mzIdentML file, it will output a fairly useless summary, and make a copy to annoy you." << endl;
  cout << "USAGE:  tpp2mzid [options] <pepXML|protXML|mzid>" << endl;
  cout << "OPTIONS:" << endl;
  cout << "  -d <string>: decoy prefix to indicate FASTA sequences that are decoy sequences. Default value is 'DECOY'." << endl;
  cout << "  -v# : optional version specification. " << endl;
  cout << "        -v1 = 1.1" << endl;
  cout << "        -v2 = 1.2" << endl;
}
