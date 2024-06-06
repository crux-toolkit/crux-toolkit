#include "ProtXMLParser.h"

// Static callback handlers
static void ProtXMLParser_startElementCallback(void *data, const XML_Char *el, const XML_Char **attr) {
	((ProtXMLParser*) data)->startElement(el, attr);
}

static void ProtXMLParser_endElementCallback(void *data, const XML_Char *el){
	((ProtXMLParser*) data)->endElement(el);
}

static void ProtXMLParser_charactersCallback(void *data, const XML_Char *s, int len){
	((ProtXMLParser*) data)->characters(s, len);
}

ProtXMLParser::ProtXMLParser(){
	parser = XML_ParserCreate(NULL);
	XML_SetUserData(parser, this);
	XML_SetElementHandler(parser, ProtXMLParser_startElementCallback, ProtXMLParser_endElementCallback);
	XML_SetCharacterDataHandler(parser, ProtXMLParser_charactersCallback);

  bIndistinguishablePeptide=false;
  bIndistinguishableProtein=false;

  //Build universal modification table
  //createModTable();
  pepXML="";
  database="";
}

ProtXMLParser::~ProtXMLParser(){
	XML_ParserFree(parser);
}

ProtXMLEntry& ProtXMLParser::operator[ ](const size_t& i){
  return vProteins[i];
}

void ProtXMLParser::characters(const XML_Char *s, int len) {
	//nothing is held here.
}

void ProtXMLParser::endElement(const XML_Char *el) {

	if(isElement("indistinguishable_peptide",el)){
		protein.peptides->push_back(peptide);
 		bIndistinguishablePeptide=true;

  } else if(isElement("indistinguishable_protein",el)){
    bIndistinguishableProtein=false;

  } else if(isElement("peptide",el)){

		if(!bIndistinguishablePeptide) protein.peptides->push_back(peptide);  //this is a redundant push if all peptide tags contain indistinguishable_peptide tags.
    bIndistinguishablePeptide=false;
		peptide.sequence="";
    peptide.modSequence="";
		peptide.nonDegenerate=false;
		peptide.evidence=false;
    peptide.charge=0;
    peptide.calcNeutralPepMass=0;
    peptide.initProb=0;
    peptide.nspProb=0;
    peptide.length=0;
    if(peptide.mods!=NULL) delete [] peptide.mods;
    peptide.mods=NULL;

	} else if(isElement("protein", el)) {
		/*if(protein.probability>0)*/ vProteins.push_back(protein);
		protein.probability=0;
    protein.stPeter.clear();
		protein.peptides->clear();
		protein.altProteinNames->clear();
		protein.proteinName="";
		probability=0;

	} else if(isElement("protein_group", el)) {
		protein.groupID=0;
		protein.subID=0;

  } else if(isElement("proteinprophet_details", el)){
    //new protein prophet in TPP5 reports probabilities/errors in reverse order of previous versions.
    //must reverse this to look like old versions so that all versions are on the same system and remain compatible.
    if(vError.size()>1 && vError[0].prob>0){
      size_t sz=vError.size();
      size_t i=sz-1;
      while(true){
        i--;
        vError.push_back(vError[i]);
        if(i==0) break;
      }
      vError.erase(vError.begin(),vError.begin()+sz-1);
    }
  }

}


void ProtXMLParser::startElement(const XML_Char *el, const XML_Char **attr){
	string s;

  if (isElement("annotation",el)){
    s = getAttrValue("protein_description", attr);
    if(s.size()>0 && !bIndistinguishableProtein) protein.proteinDescription=s;

  } else if (isElement("indistinguishable_peptide",el)){
		s = getAttrValue("calc_neutral_pep_mass", attr);
    peptide.calcNeutralPepMass=atof(&s[0]);
    s = getAttrValue("charge", attr);
    peptide.charge=atoi(&s[0]);
    peptide.modSequence=peptide.sequence;
    peptide.nMod=0;
    peptide.cMod=0;
    if(peptide.mods!=NULL) delete [] peptide.mods;
    peptide.mods = new int[peptide.length];
    for(int i=0;i<peptide.length;i++) peptide.mods[i]=0;

  } else if (isElement("indistinguishable_protein",el)){
    bIndistinguishableProtein=true;
    s = getAttrValue("protein_name",attr);
    //if(vProteinNames[protein.proteinNameIndex].find("RAND")==0 && s[0]!='R') {
    //	cout << "Not RAND: " << &s[0] << "\t" << &vProteinNames[protein.proteinNameIndex][0] << endl;
    //}
    protein.altProteinNames->push_back(s);
	
  } else if (isElement("modification_info",el)){
    s = getAttrValue("modified_peptide", attr);
    peptide.nMod=atof(getAttrValue("mod_nterm_mass",attr));
    peptide.cMod=atof(getAttrValue("mod_cterm_mass", attr)); 
    for(int i=0;i<peptide.length;i++) peptide.mods[i]=0;
    if(s.size()>0) parseModPeptide(s,peptide); //peptide.modSequence=s;

  } else if (isElement("mod_aminoacid_mass",el)){
    //double check mods. sometimes static mods are not recorded.
    s = getAttrValue("position", attr);
    int pos = atoi(&s[0]);
    if (peptide.mods[pos - 1] == 0){
      s = getAttrValue("mass", attr);
      double m=atof(&s[0]);
      int x = (int)m;
      peptide.mods[pos-1]=x;

      //rebuild mod sequence
      size_t i;
      char modStr[32];
      string sMod;

      if (peptide.nMod > 0) {
        sprintf(modStr, "n[%d]", peptide.nMod);
        peptide.modSequence = modStr;
      } else {
        peptide.modSequence = "";
      }
      for (i = 0; i<peptide.sequence.size(); i++){
        peptide.modSequence += peptide.sequence[i];
        if (peptide.mods[i]>0){
          sprintf(modStr, "[%d]", peptide.mods[i]);
          peptide.modSequence += modStr;
        }
      }

    }

  } else if (isElement("parameter",el)){
    s = getAttrValue("name", attr);
    if(s.compare("prot_length")==0) {
      s = getAttrValue("value", attr);
      protein.length = atoi(&s[0]);
    }

	} else if (isElement("peptide",el)){
		peptide.sequence = getAttrValue("peptide_sequence", attr);
    peptide.modSequence = peptide.sequence; //modified sequence is the same until a modification is defined.
    peptide.length = (int)peptide.sequence.size();
    peptide.nMod=0;
    peptide.cMod=0;
    if(peptide.mods!=NULL) delete [] peptide.mods;
    peptide.mods = new int[peptide.length];
    for(int i=0;i<peptide.length;i++) peptide.mods[i]=0;
		s = getAttrValue("is_nondegenerate_evidence", attr);
		if(s.size()>0 && s[0]=='Y') peptide.nonDegenerate=true;
		else peptide.nonDegenerate=false;
		s = getAttrValue("is_contributing_evidence", attr);
		if(s.size()>0 && s[0]=='Y') peptide.evidence=true;
		else peptide.evidence=false;
    s = getAttrValue("nsp_adjusted_probability", attr);
    peptide.nspProb=atof(&s[0]);
    s = getAttrValue("n_instances", attr);
    peptide.instances=atoi(&s[0]);
    s = getAttrValue("charge", attr);
    peptide.charge=atoi(&s[0]);
    s = getAttrValue("calc_neutral_pep_mass", attr);
    peptide.calcNeutralPepMass=atof(&s[0]);
    s = getAttrValue("initial_probability", attr);
    peptide.initProb=atof(&s[0]);

	} else if (isElement("protein",el)){
		protein.subID++;
		s = getAttrValue("probability", attr);
		protein.probability=atof(&s[0]);
		protein.peptides->clear();
		
		/*
		s = getAttrValue("unique_stripped_peptides", attr);
		strcpy(str,&s[0]);
		tok=strtok(str,"+\n");
		while(tok!=NULL){
			s=tok;
			vPeptides.push_back(s);
			protein.peptides->push_back(vPeptides.size()-1);
			tok=strtok(NULL,"+\n");
		}
		*/

		protein.proteinName = getAttrValue("protein_name",attr);
    s = getAttrValue("percent_coverage",attr);
    protein.coverage = (float)atof(&s[0]);

	} else if (isElement("protein_group",el)){
		s = getAttrValue("group_number",attr);
		protein.groupID=atoi(&s[0]);

	} else if (isElement("protein_summary_data_filter",el)){
    ProtXMLError e;
    s = getAttrValue("false_positive_error_rate",attr);
    e.error=atof(&s[0]);
    s = getAttrValue("min_probability",attr);
    e.prob=atof(&s[0]);
    vError.push_back(e);

  } else if (isElement("protein_summary_header",el)){
    pepXML = getAttrValue("source_files",attr);
    database = getAttrValue("reference_database",attr);

  } else if(isElement("StPeterQuant",el)){
	s = getAttrValue("dSIn", attr);
    protein.stPeter.dSIn = atof(&s[0]);							   								   
    s = getAttrValue("SIn",attr);
    protein.stPeter.SIn=atof(&s[0]);
    s = getAttrValue("SI", attr);
    protein.stPeter.SI = atof(&s[0]);
    s = getAttrValue("ng", attr);
    protein.stPeter.ng = atof(&s[0]);
    s = getAttrValue("ngC", attr);
    protein.stPeter.ngCounts = atof(&s[0]);
    s = getAttrValue("NSAF", attr);
    protein.stPeter.NSAF = atof(&s[0]);
    s = getAttrValue("counts", attr);
    protein.stPeter.counts = atoi(&s[0]);

  }
}

bool ProtXMLParser::parse(const char* fileName) {

  XML_ParserFree(parser);
  parser = XML_ParserCreate(NULL);
	XML_SetUserData(parser, this);
	XML_SetElementHandler(parser, ProtXMLParser_startElementCallback, ProtXMLParser_endElementCallback);
	XML_SetCharacterDataHandler(parser, ProtXMLParser_charactersCallback);

  vError.clear();
  vProteins.clear();
	
	FILE* fptr=fopen(fileName,"rt");
	if (fptr == NULL){
		cerr << "Error parse(): No open file. Cannot open " << fileName << endl;
		return false;
	}

	char buffer[16384];
	int readBytes = 0;
	bool success = true;
	int chunk=0;

	while (success && (readBytes = (int) fread(buffer, 1, sizeof(buffer), fptr)) != 0){
		success = (XML_Parse(parser, buffer, readBytes, false) != 0);
	}
	success = success && (XML_Parse(parser, buffer, 0, true) != 0);

	if (!success) {
		XML_Error error = XML_GetErrorCode(parser);

		cerr << fileName << "(" << XML_GetCurrentLineNumber(parser) << ") : error " << (int) error << ": ";
		switch (error) {
			case XML_ERROR_SYNTAX:
			case XML_ERROR_INVALID_TOKEN:
			case XML_ERROR_UNCLOSED_TOKEN:
				cerr << "Syntax error parsing XML.";
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

string ProtXMLParser::getDatabase(){
  return database;
}

string ProtXMLParser::getPepXML(){
  return pepXML;
}

double ProtXMLParser::getProbability(double err){
  size_t i;
  for(i=0;i<vError.size();i++){
    if(vError[i].error==err) return vError[i].prob;
    if(vError[i].error<err) break;
  }
  double slope=(vError[i-1].error-vError[i].error)/(vError[i-1].prob-vError[i].prob);
  double intercept=vError[i].error-(vError[i].prob*slope);
  return (err-intercept)/slope;
}

int ProtXMLParser::size(){
	return (int)vProteins.size();
}

void ProtXMLParser::sortProbabilityRev(){
	qsort(&vProteins[0],vProteins.size(),sizeof(ProtXMLEntry),compareProbRev);
}


/* save for later?
void ProtXMLParser::createModTable(){
  for(int i=0;i<128;i++) modTable[i]=0;
  modTable[0]=0;            //no modification
  modTable[1]=57.0214611;   //Carbamidomethyl
  modTable[2]=15.9949141;   //Oxidized
  modTable[3]=-18.0105633;  //Loss of water
  modTable[4]=-17.026547;   //Loss of ammonia
  modTable[5]=79.966331;    //Phosphorylated
  modTable[6]=0.984016;     //Deamidation
  modTable[7]=39.9949141;   //Loss of ammonia+Carbamidomethyl
  modTable[8]=42.0105633;   //Acetylation
  modTable[9]=79.9663289;   //Phophorylation
  modTable[10]=45.987721;   //Beta-methylthiolation
  modTable[11]=24.0;        //Acetylation+Loss of water
  modTable[12]=24.9840163;  //Acetylation+Loss of ammonia
  modTable[13]=8.014199;    //SILAC heavy mass (typically lysine)
  modTable[14]=10.008269;   //SILAC heavy mass (typically argenine)
  modTable[15]=145.019749;  //Biotin
}
*/

/* save for later?
int ProtXMLParser::modLookup(char aa, double mass){

  switch(aa){
    case 'C':
      if(fabs(mass-160)<0.001) return 1;
      if(fabs(mass-143)<0.001) return 7;
      if(fabs(mass-149)<0.001) return 10;
      if (fabs(mass - 160.030649)<0.001) return 1;
      break;
    case 'E':
      if(fabs(mass-111)<0.001) return 3;
      break;
    case 'K':
      if(fabs(mass-136)<0.001) return 13;
      if(fabs(mass-273)<0.001) return 15;
      break;
    case 'M':
      if(fabs(mass-147)<0.001) return 2;
      break;
    case 'n':
      if(fabs(mass-43)<0.001) return 8;  //Don't know why some people use 43 for acetylation
      if(fabs(mass-25)<0.001) return 11; //Ditto
      if(fabs(mass-26)<0.001) return 12; //Ditto
      break;
    case 'N':
      if(fabs(mass-115)<0.001) return 6;
      break;
    case 'Q':
      if(fabs(mass-129)<0.001) return 6;
      if(fabs(mass-111)<0.001) return 4;
      break;
    case 'R':
      if(fabs(mass-157)<0.001) return 6;
      if(fabs(mass-166)<0.001) return 14;
      break;
    case 'S':
      if(fabs(mass-167)<0.001) return 9;
      break;
    case 'T':
      if(fabs(mass-181)<0.001) return 9;
      break;
    case 'Y':
      if(fabs(mass-243)<0.001) return 9;
      break;
    case 'W':
      if (fabs(mass-202)<0.001) return 2;
      break;
    default:
      break;
  }
  return -1;

}
*/

void ProtXMLParser::parseModPeptide(string& s, ProtXMLPeptide& p){
  size_t i;
  int j;
  char aa;
  char modStr[32];
  string sMod;
  int modMass;
  bool bMod=false;

  
  j=0;
  for(i=0;i<s.size();i++){
    if(s[i]=='['){
      bMod=true;
      aa=s[i-1];
      sMod="";
    } else if(s[i]==']'){
      bMod=false;
      modMass=atoi(&sMod[0]);
      if(aa=='n') {
        p.nMod=modMass;
        j--; //n-terminus needs to be decremented so as not to throw off amino acid positions
      } else {
        p.mods[j-1]=modMass;
      }      
    } else if(bMod){
      sMod+=s[i];
    } else {
      //p.modSequence+=s[i];
      j++;
    }
  }

  p.modSequence="";
  if(p.nMod>0) {
    sprintf(modStr, "n[%d]", p.nMod);
    p.modSequence+=modStr;
  }
  for(i=0;i<p.sequence.size();i++){
    p.modSequence+=p.sequence[i];
    if(p.mods[i]>0){
      sprintf(modStr,"[%d]",p.mods[i]);
      p.modSequence+=modStr;
    }
  }

}


/* For the qsort */
int ProtXMLParser::compareProbRev(const void *p1, const void *p2){
  const ProtXMLEntry d1 = *(ProtXMLEntry *)p1;
  const ProtXMLEntry d2 = *(ProtXMLEntry *)p2;
	if(d1.probability>d2.probability) return -1;
	else if(d1.probability<d2.probability) return 1;
  else return 0;
}
