/*
Copyright 2017, Michael R. Hoopmann, Institute for Systems Biology
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

#include "ProtXMLParser2.h"

using namespace std;

// Static callback handlers
static void ProtXMLParser2_startElementCallback(void *data, const XML_Char *el, const XML_Char **attr) {
  ((ProtXMLParser2*)data)->startElement(el, attr);
}

static void ProtXMLParser2_endElementCallback(void *data, const XML_Char *el){
  ((ProtXMLParser2*)data)->endElement(el);
}

static void ProtXMLParser2_charactersCallback(void *data, const XML_Char *s, int len){
  ((ProtXMLParser2*)data)->characters(s, len);
}

ProtXMLParser2::ProtXMLParser2(){
  algorithmName.clear();
  sourceFiles = new vector<string>;
  version.clear();
  parser = XML_ParserCreate(NULL);
  bFixMods=false;
  XML_SetUserData(parser, this);
  XML_SetElementHandler(parser, ProtXMLParser2_startElementCallback, ProtXMLParser2_endElementCallback);
  XML_SetCharacterDataHandler(parser, ProtXMLParser2_charactersCallback);
}

ProtXMLParser2::~ProtXMLParser2(){
  delete sourceFiles;
  XML_ParserFree(parser);
}

cpxProteinGroup& ProtXMLParser2::operator[](const size_t index){
  return proteinGroup[index];
}

void ProtXMLParser2::characters(const XML_Char *s, int len) {
}

void ProtXMLParser2::endElement(const XML_Char *el) {

  string s;

  if (isElement("indistinguishable_peptide", el)){
    if (elements.back() != indistinguishable_peptide) cout << "Error in elements" << endl;
    else elements.pop_back();

  } else if (isElement("modification_info",el)){
    if (bFixMods){
      bFixMods=false;
      char str[32];
      string s;
      string pep = proteinGroup.back().protein->back().peptide->back().peptideSequence;
      s.clear();

      cpxModInfo m;
      if (elements.back() == indistinguishable_peptide){
        m = proteinGroup.back().protein->back().peptide->back().indistinguishablePeptide->back().modificationInfo;
      } else if (elements.back() == peptide){
        m = proteinGroup.back().protein->back().peptide->back().modificationInfo;
      }

      if (m.nTermMod != 0){
        sprintf(str, "n[%.0lf]", m.nTermMod);
        s += str;
      }
      for (size_t i = 0; i < pep.size(); i++){
        s += pep[i];
        for (size_t j = 0; j < m.modAminoAcidMass->size(); j++){
          if (m.modAminoAcidMass->at(j).pos == i + 1){
            sprintf(str, "[%.0lf]", m.modAminoAcidMass->at(j).mass);
            s+=str;
            break;
          }
        }
      }
      if (m.nTermMod != 0){
        sprintf(str, "c[%.0lf]", m.cTermMod);
        s += str;
      }

      if (elements.back() == indistinguishable_peptide){
        proteinGroup.back().protein->back().peptide->back().indistinguishablePeptide->back().modificationInfo.modifiedPeptide=s;
      } else if (elements.back() == peptide){
        proteinGroup.back().protein->back().peptide->back().modificationInfo.modifiedPeptide=s;
      }

    }

  } else if (isElement("peptide", el)){
    if (elements.back() != peptide) cout << "Error in elements" << endl;
    else elements.pop_back();

  } else if (isElement("protein_group", el)){
  }

}

void ProtXMLParser2::startElement(const XML_Char *el, const XML_Char **attr){

  string s;

  if (isElement("indistinguishable_peptide", el)){
    elements.push_back(indistinguishable_peptide);
    cpxIndPeptide p;
    p.peptideSequence = getAttrValue("peptide_sequence", attr);
    s = getAttrValue("charge", attr);
    p.charge = atoi(&s[0]);
    s = getAttrValue("calc_neutral_pep_mass", attr);
    p.calcNeutralPepMass = atof(&s[0]);
    proteinGroup.back().protein->back().peptide->back().indistinguishablePeptide->push_back(p);

  } else if (isElement("modification_info", el)){
    if (elements.back() == indistinguishable_peptide){
      proteinGroup.back().protein->back().peptide->back().indistinguishablePeptide->back().modificationInfo.modifiedPeptide = getAttrValue("modified_peptide", attr);
    } else if (elements.back() == peptide){
      proteinGroup.back().protein->back().peptide->back().modificationInfo.modifiedPeptide = getAttrValue("modified_peptide", attr);
    }

  } else if (isElement("mod_aminoacid_mass", el)){
    bFixMods=true;
    spxMod m;
    s = getAttrValue("position",attr);
    m.pos = atoi(&s[0]); 
    s = getAttrValue("mass", attr);
    m.mass = atof(&s[0]);
    if (elements.back() == indistinguishable_peptide){
      proteinGroup.back().protein->back().peptide->back().indistinguishablePeptide->back().modificationInfo.modAminoAcidMass->push_back(m);
    } else if (elements.back() == peptide){
      proteinGroup.back().protein->back().peptide->back().modificationInfo.modAminoAcidMass->push_back(m);
    }

  } else if (isElement("peptide", el)){
    elements.push_back(peptide);
    cpxPeptide p;
    p.peptideSequence = getAttrValue("peptide_sequence", attr);
    s = getAttrValue("charge", attr);
    p.charge = atoi(&s[0]);
    s = getAttrValue("initial_probability", attr);
    p.initialProbability = atof(&s[0]);
    s = getAttrValue("nsp_adjusted_probability", attr);
    p.nspAdjustedProbability = atof(&s[0]);
    s = getAttrValue("fpkm_adjusted_probability", attr);
    p.fpkmAdjustedProbability = atof(&s[0]);
    s = getAttrValue("weight", attr);
    p.weight = atof(&s[0]);
    s = getAttrValue("is_nondegenerate_evidence", attr);
    if (s.size()>0 && s[0] == 'Y') p.isNondegenerateEvidence=true;
    else p.isNondegenerateEvidence=false;
    s = getAttrValue("n_enzymatic_termini", attr);
    p.nEnzymaticTermini = atoi(&s[0]);
    s = getAttrValue("n_sibling_peptides", attr);
    p.nSiblingPeptides = atof(&s[0]);
    s = getAttrValue("n_sibling_peptides_bin", attr);
    p.nSiblingPeptidesBin = atoi(&s[0]);
    s = getAttrValue("n_instances", attr);
    p.nInstances = atoi(&s[0]);
    s = getAttrValue("exp_tot_instances", attr);
    p.expTotInstances = atof(&s[0]);
    s = getAttrValue("is_contributing_evidence", attr);
    if (s.size()>0 && s[0] == 'Y') p.isContributingEvidence = true;
    else p.isContributingEvidence = false;
    proteinGroup.back().protein->back().peptide->push_back(p);

  } else if (isElement("program_details", el)){
    algorithmName = getAttrValue("analysis",attr);
    version = getAttrValue("version", attr);

  } else if (isElement("protein", el)){
    cpxProtein p;
    p.proteinName = getAttrValue("protein_name", attr);
    s = getAttrValue("n_indistinguishable_proteins", attr);
    p.nIndistinguishableProteins = atoi(&s[0]);
    s = getAttrValue("probability", attr);
    p.probability = atof(&s[0]);
    s = getAttrValue("percent_coverage", attr);
    p.percentCoverage = atof(&s[0]);
    p.groupSiblingID = getAttrValue("group_sibling_id",attr);
    s = getAttrValue("total_number_peptides", attr);
    p.totalNumberPeptides = atoi(&s[0]);
    s = getAttrValue("total_number_distinct_peptides", attr);
    p.totalNumberDistinctPeptides = atoi(&s[0]);
    s = getAttrValue("pct_spectrum_ids", attr);
    p.pctSpectrumIDs = atof(&s[0]);
    s = getAttrValue("confidence", attr);
    p.confidence = atof(&s[0]);

    s = getAttrValue("unique_stripped_peptides", attr);
    proteinGroup.back().protein->push_back(p);

  } else if (isElement("protein_group", el)){
    cpxProteinGroup pg;
    s = getAttrValue("group_number", attr);
    pg.groupNumber = atoi(&s[0]);
    s = getAttrValue("probability", attr);
    pg.probability = atof(&s[0]);
    proteinGroup.push_back(pg);
  
  } else if (isElement("protein_summary_header", el)){
    string s2;
    size_t pos,x;
    sourceFiles->clear();
    s = getAttrValue("source_files", attr);
    pos=0;
    x = s.find(' ',pos);
    while (x != string::npos){
      s2 = s.substr(pos,x-pos-1);
      sourceFiles->push_back(s2);
      pos=x+1;
      x = s.find(' ',pos);
    }
    s2 = s.substr(pos);
    sourceFiles->push_back(s2);

  } else if (isElement("StPeter_analysis_summary", el)){
    STPSummary.version = getAttrValue("version", attr);
    s = getAttrValue("probability", attr);
    STPSummary.probability=atof(s.c_str());
    s = getAttrValue("FDR", attr);
    STPSummary.FDR = atof(s.c_str());
    s = getAttrValue("sampleLoad", attr);
    STPSummary.sampleLoad = atof(s.c_str());
    s = getAttrValue("tolerance", attr);
    STPSummary.tolerance = atof(s.c_str());
    s = getAttrValue("FDR", attr);
    if(s.compare("yes")==0) STPSummary.bDegen=true;
    else STPSummary.bDegen=false;

  } else if (isElement("StPeterQuant", el)){
    s = getAttrValue("SI", attr);
    if(s.size()>0) proteinGroup.back().protein->back().stPeter.SI = atof(s.c_str());
    s = getAttrValue("SIn", attr);
    if (s.size()>0) proteinGroup.back().protein->back().stPeter.SIn = atof(s.c_str());
    s = getAttrValue("NSAF", attr);
    if (s.size()>0) proteinGroup.back().protein->back().stPeter.NSAF = atof(s.c_str());
    s = getAttrValue("ng", attr);
    if (s.size()>0) proteinGroup.back().protein->back().stPeter.ng = atof(s.c_str());
    s = getAttrValue("ngC", attr);
    if (s.size()>0) proteinGroup.back().protein->back().stPeter.ngC = atof(s.c_str());
    s = getAttrValue("counts", attr);
    if (s.size()>0) proteinGroup.back().protein->back().stPeter.counts = atoi(s.c_str());

  } else if (isElement("StPeterQuant_peptide", el)){
    spxSTPQPeptide stpPep;
    stpPep.sequence = getAttrValue("sequence", attr);
    s = getAttrValue("charge", attr);
    stpPep.charge = atoi(s.c_str());
    s = getAttrValue("SI", attr);
    stpPep.SI = atof(s.c_str());
    s = getAttrValue("SC", attr);
    stpPep.counts = atoi(s.c_str());
    proteinGroup.back().protein->back().stPeter.stpPeptide->push_back(stpPep);

  }

}

string ProtXMLParser2::getAlgorithm(){
  return algorithmName;
}

string ProtXMLParser2::getSourceFile(size_t index){
  if (index >= sourceFiles->size()) return "";
  return sourceFiles->at(index);
}

size_t ProtXMLParser2::getSourceFileCount(){
  return sourceFiles->size();
}

sSTPSummary ProtXMLParser2::getStPeterSummary(){
  if(hasStPeter())  return STPSummary;
  sSTPSummary blank;
  return blank;
}

string ProtXMLParser2::getVersion(){
  return version;
}

bool ProtXMLParser2::hasStPeter(){
  if(STPSummary.version.compare("0")!=0) return true;
  return false;
}

bool ProtXMLParser2::readFile(const char* fn) {

  XML_ParserFree(parser);
  parser = XML_ParserCreate(NULL);
  XML_SetUserData(parser, this);
  XML_SetElementHandler(parser, ProtXMLParser2_startElementCallback, ProtXMLParser2_endElementCallback);
  XML_SetCharacterDataHandler(parser, ProtXMLParser2_charactersCallback);

  // clear data
  algorithmName.clear();
  version.clear();
  proteinGroup.clear();
  STPSummary.version="0";

  FILE* fptr = fopen(fn, "rt");
  if (fptr == NULL){
    cerr << "Error parse(): No open file." << endl;
    return false;
  }

  char buffer[16384];
  int readBytes = 0;
  bool success = true;
  int chunk = 0;
  killRead = false;

  while (success && (readBytes = (int)fread(buffer, 1, sizeof(buffer), fptr)) != 0){
    success = (XML_Parse(parser, buffer, readBytes, false) != 0);
    if (killRead){
      fclose(fptr);
      return false;
    }
  }
  success = success && (XML_Parse(parser, buffer, 0, true) != 0);

  if (!success) {
    XML_Error error = XML_GetErrorCode(parser);

    cerr << fn << "(" << XML_GetCurrentLineNumber(parser) << ") : error " << (int)error << ": ";
    switch (error) {
    case XML_ERROR_SYNTAX:
    case XML_ERROR_INVALID_TOKEN:
    case XML_ERROR_UNCLOSED_TOKEN:
      cerr << "Syntax error parsing XML.";
      break;
    case XML_ERROR_TAG_MISMATCH:
      cerr << "XML tag mismatch.";
      break;
    case XML_ERROR_DUPLICATE_ATTRIBUTE:
      cerr << "XML duplicate attribute.";
      break;
    case XML_ERROR_JUNK_AFTER_DOC_ELEMENT:
      cerr << "XML junk after doc element.";
      break;
    default:
      cerr << "XML Parsing error.";
      break;
    }
    cerr << "\n";
    fclose(fptr);
    return false;
  }

  fclose(fptr);
  return true;
}

size_t ProtXMLParser2::size(){
  return proteinGroup.size();
}

bool ProtXMLParser2::writeFile(const char* fn){
  return true;
}