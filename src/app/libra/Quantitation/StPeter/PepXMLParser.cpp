#include "PepXMLParser.h"

// Static callback handlers
static void PepXMLParser_startElementCallback(void *data, const XML_Char *el, const XML_Char **attr) {
	((PepXMLParser*) data)->startElement(el, attr);
}

static void PepXMLParser_endElementCallback(void *data, const XML_Char *el){
	((PepXMLParser*) data)->endElement(el);
}

static void PepXMLParser_charactersCallback(void *data, const XML_Char *s, int len){
	((PepXMLParser*) data)->characters(s, len);
}

PepXMLParser::PepXMLParser(){
  bIprophet=false;
  bMassCalc=false;
  rank=0;
  bHit=false;

  for(int i=0;i<128;i++) aaMass[i]=0;
  aaMass['A']= 71.03711;
  aaMass['C'] = 103.00918;
  aaMass['D'] = 115.02694;
  aaMass['E'] = 129.04259;
  aaMass['F'] = 147.06841;
  aaMass['G'] = 57.02146; 
  aaMass['H'] = 137.05891;
  aaMass['I'] = 113.08406;
  aaMass['K'] = 128.09496;
  aaMass['L'] = 113.08406;
  aaMass['M'] = 131.04048;
  aaMass['N'] = 114.04293;
  aaMass['P'] = 97.05276; 
  aaMass['Q'] = 128.05858;
  aaMass['R'] = 156.10111;
  aaMass['S'] = 87.03203; 
  aaMass['T'] = 101.04768;
  aaMass['V'] = 99.06841; 
  aaMass['W'] = 186.07931; 
  aaMass['Y'] = 163.06333; 

	parser = XML_ParserCreate(NULL);
	XML_SetUserData(parser, this);
	XML_SetElementHandler(parser, PepXMLParser_startElementCallback, PepXMLParser_endElementCallback);
	XML_SetCharacterDataHandler(parser, PepXMLParser_charactersCallback);
}

PepXMLParser::~PepXMLParser(){
	XML_ParserFree(parser);
}

PepXMLEntry& PepXMLParser::operator[ ](const unsigned int& i){
  return vPeptides[i];
}

bool PepXMLParser::addMod(char aa, double massDiff, int pos){
  PepXMLMod mod;

  mod.mass = 0;
  switch (aa){
  case 'A': mod.mass = massDiff-71.03711; break;
  case 'C': mod.mass = massDiff-103.00918; break;
  case 'D': mod.mass = massDiff - 115.02694; break;
  case 'E': mod.mass = massDiff - 129.04259; break;
  case 'F': mod.mass = massDiff - 147.06841; break;
  case 'G': mod.mass = massDiff - 57.02146; break;
  case 'H': mod.mass = massDiff - 137.05891; break;
  case 'I': mod.mass = massDiff - 113.08406; break;
  case 'K': mod.mass = massDiff - 128.09496; break;
  case 'L': mod.mass = massDiff - 113.08406; break;
  case 'M': mod.mass = massDiff - 131.04048; break;
  case 'n': mod.mass = massDiff - 1.007825; break; //Note that MS-GF+ pepXMLs may have bad n-term values!!!
  case 'N': mod.mass = massDiff - 114.04293; break;
  case 'P': mod.mass = massDiff - 97.05276; break;
  case 'Q': mod.mass = massDiff - 128.05858; break;
  case 'R': mod.mass = massDiff - 156.10111; break;
  case 'S': mod.mass = massDiff - 87.03203; break;
  case 'T': mod.mass = massDiff - 101.04768; break;
  case 'V': mod.mass = massDiff - 99.06841; break;
  case 'W': mod.mass = massDiff - 186.07931; break;
  case 'Y': mod.mass = massDiff - 163.06333; break;
  default:
    break;
  }

  /* old way of doing things to reconcile different search params
  //Look up actual modification mass: these need to be adjusted!!
  mod.mass=0;
  switch(aa){
    case 'C':
      if(fabs(massDiff-160.0307)<0.001) mod.mass=57.021464;
      if(fabs(massDiff-143.0042)<0.001) mod.mass=39.9949141;
      if(fabs(massDiff-148.9969)<0.001) mod.mass=45.987721;
      break;
    case 'E':
      if(fabs(massDiff-111.0320)<0.001) mod.mass=-18.0105633;
      break;
    case 'K':
      if(fabs(massDiff-136.109162)<0.001) mod.mass=8.014199;
      if(fabs(massDiff-273.114712)<0.001) mod.mass=145.019749;
      break;
    case 'M':
      if(fabs(massDiff-147.0353937)<0.001) mod.mass=15.9949141;
      if(fabs(massDiff-147.039885)<0.001) mod.mass=15.9949141;
      break;
    case 'n':
      if(fabs(massDiff-43.0184)<0.001) mod.mass=42.0105633;
      if(fabs(massDiff-25.0078)<0.001) mod.mass=24.0;
      if(fabs(massDiff-25.9918)<0.001) mod.mass=24.9840163;
      if(fabs(massDiff+16.0187)<0.001) mod.mass=-17.026549;  //This fixes a bug in pepXML generated from MS-GF+
      if(fabs(massDiff+17.0027)<0.001) mod.mass=-18.0105633; //This fixes a bug in pepXML generated from MS-GF+
      break;
    case 'N':
      if(fabs(massDiff-115.0269)<0.001) mod.mass=0.984016;
      break;
    case 'Q':
      if(fabs(massDiff-129.0429)<0.001) mod.mass=0.984016;
      if(fabs(massDiff-111.0321)<0.001) mod.mass=-17.026549;
      break;
    case 'R':
      if(fabs(massDiff-157.0851)<0.001) mod.mass=0.984016;
      if(fabs(massDiff-166.109380)<0.001) mod.mass=10.008269;
      break;
    case 'S':
      if(fabs(massDiff-166.9984)<0.001) mod.mass=79.9663289;
      break;
    case 'T':
      if(fabs(massDiff-181.0140)<0.001) mod.mass=79.9663289;
      break;
    case 'Y':
      if(fabs(massDiff-243.0297)<0.001) mod.mass=79.9663289;
      break;
    case 'W':
      if (fabs(massDiff - 202.0742)<0.001) mod.mass = 15.9949141;
      break;
    default:
      break;
  }
  */

  if(mod.mass==0){
    cout << "Unknown amino acid, mass modification combination: " << aa << " = " << massDiff << endl;
    return false;
  }

  mod.pos=pos;
  aaModList.push_back(mod);
  
  return true;

}

void PepXMLParser::characters(const XML_Char *s, int len) {
	//nothing is held here.
}

void PepXMLParser::endElement(const XML_Char *el) {

	string s,s2;
  char str[32];
	unsigned int u,n;
  double modMass;

	if(isElement("search_hit", el))	{

    s="";
    s2="";
    for(n=0;n<aaModList.size();n++){ //look for n-terminal mods
      if(aaModList[n].pos==0){
        s+='n';
        s2+='n';
        sprintf(str, "[%d]", (int)(aaModList[n].mass+1.007825+0.5)); //add 0.5 for rounding up
        s += str;
        sprintf(str, "[%.4lf]", aaModList[n].mass);
        s2 += str;
        break;
      }
    }
    for(u=0;u<peptide.peptide.size();u++){
      s+=peptide.peptide[u];
      s2+=peptide.peptide[u];
      modMass=0;
      for(n=0;n<aaModList.size();n++){
        if(aaModList[n].pos==u+1) modMass+=aaModList[n].mass;
      }
      if(modMass!=0){
        sprintf(str, "[%d]", (int)(aaMass[peptide.peptide[u]]+modMass+0.5));//add 0.5 for rounding up
        s+=str;
      }
      if (modMass != 0){
        sprintf(str, "[%.4lf]", modMass);
        s2 += str;
      }
    }
		peptide.modifiedPeptidePlus=s;
    peptide.modifiedPeptide=s2;

    if(rank==1 && !bHit) {
      vPeptides.push_back(peptide);
      bHit=true;
    }
    rank=0;

	} else if(isElement("spectrum_query", el))	{
		peptide.charge=0;
		peptide.expect=0;
		peptide.iProphetProbability=0;
		peptide.monoMass=0;
		peptide.precursorMonoMass=0;
    peptide.probability=0;
		peptide.RTime=0;
		peptide.scanNum=0;
		peptide.peptide="";
    peptide.modifiedPeptide="";
    peptide.modifiedPeptidePlus="";
    peptide.label="";
    peptide.xLabel="";
    peptide.protein="";
    peptide.fileID=0;
    bHit=false;
	}
}


void PepXMLParser::startElement(const XML_Char *el, const XML_Char **attr){
  char          c;
  double        d1;
  string        s;
  unsigned int  u;
  
  if (isElement("error_point",el)){
    PepXMLError e;
    s = getAttrValue("error",attr);
    e.error=atof(&s[0]);
    s = getAttrValue("min_prob",attr);
    e.prob=atof(&s[0]);
    vError.push_back(e);
    
  } else if (isElement("interprophet_result",el)){
    s = getAttrValue("probability", attr);
		peptide.iProphetProbability=atof(&s[0]);
		
  } else if(isElement("interprophet_summary",el)){
    bIprophet=true;
    
  } else if( isElement("modification_info",el)){
    s = getAttrValue("mod_nterm_mass",attr);
    
    if(s.size()>0){
      d1 = atof(&s[0]);
      addMod('n',d1,0);
    }
    
  } else if (isElement("mod_aminoacid_mass",el)){
    s = getAttrValue("mass",attr);
    d1 = atof(&s[0]);
    s = getAttrValue("position",attr);
    c = peptide.peptide[atoi(&s[0])-1];
    addMod(c,d1,atoi(&s[0]));
    
  } else if (isElement("msms_run_summary",el)){
    s = getAttrValue("base_name",attr);
    s += getAttrValue("raw_data",attr);
    for(u=0;u<vFiles.size();u++){
      if(s.compare(vFiles[u])==0) break;
    }
    currentFileID=(int)u;
    if(u==vFiles.size()) vFiles.push_back(s);
    
  } else if (isElement("peptideprophet_result",el)){
    s = getAttrValue("probability", attr);
    peptide.probability=atof(&s[0]);
    
  } else if (isElement("search_hit",el)) {
    aaModList.clear();
    s = getAttrValue("hit_rank", attr);
    rank = atoi(&s[0]);
    s = getAttrValue("calc_neutral_pep_mass", attr);
    peptide.monoMass=atof(&s[0]);
    peptide.peptide = getAttrValue("peptide", attr);
    if(bMassCalc) peptide.monoMass=calcMonoMass(&peptide.peptide[0]);
    peptide.protein = getAttrValue("protein", attr);
    s = getAttrValue("peptide_prev_aa", attr);
    if(s.length()>0) peptide.prevAA=s[0];
    else peptide.prevAA='-';
    s = getAttrValue("peptide_next_aa", attr);
    if(s.length()>0) peptide.nextAA=s[0];
    else peptide.nextAA='-';
    
  } else if (isElement("search_score",el)) {
    s = getAttrValue("name", attr);
    if(s.compare("expect")==0){
      s = getAttrValue("value", attr);
      peptide.expect = atof(&s[0]);
    }
    if(s.compare("xcorr")==0){
      s = getAttrValue("value", attr);
      peptide.xcorr = atof(&s[0]);
    }
    
  } else if (isElement("search_summary",el)) {
    
  } else if (isElement("spectrum_query",el)) {
    s = getAttrValue("start_scan", attr);
    peptide.scanNum=atoi(&s[0]);
    s = getAttrValue("precursor_neutral_mass", attr);
    peptide.precursorMonoMass=atof(&s[0]);
    s = getAttrValue("assumed_charge", attr);
    peptide.charge=atoi(&s[0]);
    s = getAttrValue("retention_time_sec", attr);

    peptide.RTime=(float)(atof(&s[0])/60.0);
    peptide.label = getAttrValue("spectrum", attr);
    peptide.xLabel  = getAttrValue("experiment_label", attr);
    
    s = getAttrValue("spectrum", attr);
    s = s.substr(0,s.find("."));
    peptide.fileID=currentFileID;
    peptide.peptide="";
    peptide.protein="";
    peptide.expect=1000;
    peptide.iProphetProbability=0.0;
    peptide.monoMass=0.0;
    peptide.probability=0.0;
    
  } else if (isElement("terminal_modification",el)){
    
  }
}

double PepXMLParser::calcMonoMass(char *seq, bool water){

	double mass=0.0;

	int H=0;
	int C=0;
	int N=0;
	int O=0;
	int S=0;

	unsigned int i;

	for(i=0;i<strlen(seq);i++){
		switch(seq[i]){
		case 'A':	C+= 3; H+= 5; N+= 1; O+= 1;	break;
		case 'R':	C+= 6; H+= 12; N+= 4;	O+= 1; break;
		case 'N': C+= 4; H+= 6; N+= 2; O+= 2;	break;    
		case 'D': C+= 4; H+= 5; N+= 1; O+= 3; break;
		case 'C': C+= 3; H+= 5;	N+= 1; O+= 1; S+= 1; break; //IAA treated samples
		case 'Q': C+= 5; H+= 8; N+= 2; O+= 2; break;
		case 'E': C+= 5; H+= 7; N+= 1; O+= 3; break;
		case 'G': C+= 2; H+= 3; N+= 1; O+= 1; break;
		case 'H': C+= 6; H+= 7; N+= 3; O+= 1; break;
		case 'I':
    case 'L':	C+= 6; H+= 11; N+= 1; O+= 1; break;
		case 'K': C+= 6; H+= 12; N+= 2; O+= 1; break;
    case 'M': C+= 5; H+= 9; N+= 1; O+= 1; S+= 1; break;
		case 'F': C+= 9; H+= 9; N+= 1; O+= 1; break;
		case 'P': C+= 5; H+= 7; N+= 1; O+= 1; break;
		case 'S': C+= 3; H+= 5; N+= 1; O+= 2; break;
		case 'T': C+= 4; H+= 7; N+= 1; O+= 2; break;
		case 'W': C+= 11; H+= 10; N+= 2; O+= 1; break;
		case 'Y': C+= 9; H+= 9; N+= 1; O+= 2; break;
		case 'V': C+= 5; H+= 9; N+= 1; O+= 1; break;
    default: break;
		}
	}

  if(water){
	  H+=2;
	  O++;
  }

	mass+=1.0078246*H;
	mass+=12.0000000*C;
	mass+=14.0030732*N;
	mass+=15.9949141*O;
	mass+=31.972070*S;

	return mass;
}

void PepXMLParser::getFile(int index, char* str){
  strcpy(str,&vFiles[vPeptides[index].fileID][0]);
}

int PepXMLParser::getFileCount(){
  return (int)vFiles.size();
}

void PepXMLParser::getFileFromList(int index, char* str){
  strcpy(str,&vFiles[index][0]);
}

bool PepXMLParser::getIprophet(){
  return bIprophet;
}

double PepXMLParser::getProbability(double err){
  size_t i;
  for(i=0;i<vError.size();i++){
    if(vError[i].error==err) return vError[i].prob;
    if(vError[i].error>err) break;
  }
  //interpolate error
  double slope=(vError[i].error-vError[i-1].error)/(vError[i].prob-vError[i-1].prob);
  double intercept=vError[i].error-(vError[i].prob*slope);
  return (err-intercept)/slope;
}

int PepXMLParser::parse(const char* fileName) {

  XML_ParserFree(parser);
	parser = XML_ParserCreate(NULL);
	XML_SetUserData(parser, this);
	XML_SetElementHandler(parser, PepXMLParser_startElementCallback, PepXMLParser_endElementCallback);
	XML_SetCharacterDataHandler(parser, PepXMLParser_charactersCallback);

	vPeptides.clear();
  vFiles.clear();
  vError.clear();
  bIprophet=false;
	
	FILE* fptr=fopen(fileName,"rt");
	if (fptr == NULL){
		//cerr << "Error parse(): No open file." << endl;
		return -1;
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
		return -2;
	}
	
	fclose(fptr);
	return 0;
}

void PepXMLParser::setMassCalc(bool b){
  bMassCalc=b;
}

int PepXMLParser::size(){
	return (int)vPeptides.size();
}

void PepXMLParser::sortModPep(){
  qsort(&vPeptides[0],vPeptides.size(),sizeof(PepXMLEntry),compareModPep);
}

void PepXMLParser::sortScanNum(){
  qsort(&vPeptides[0],vPeptides.size(),sizeof(PepXMLEntry),compareScanNum);
}

int PepXMLParser::compareModPep(const void *p1, const void *p2){
  const PepXMLEntry d1 = *(PepXMLEntry *)p1;
  const PepXMLEntry d2 = *(PepXMLEntry *)p2;
  int i=d1.modifiedPeptide.compare(d2.modifiedPeptide);
  if(i<0) return -1;
  else if(i>0) return 1;
  else return 0;
}

int PepXMLParser::compareScanNum(const void *p1, const void *p2){
  const PepXMLEntry d1 = *(PepXMLEntry *)p1;
  const PepXMLEntry d2 = *(PepXMLEntry *)p2;
  if(d1.scanNum<d2.scanNum) return -1;
  else if(d1.scanNum>d2.scanNum) return 1;
  else return 0;
}

