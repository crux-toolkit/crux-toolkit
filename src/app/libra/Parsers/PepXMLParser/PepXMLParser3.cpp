#include "PepXMLParser3.h"

// *******************************
//    XML Parser Code
// *******************************

// Static callback handlers
static void PepXMLParser3_startElementCallback(void *data, const XML_Char *el, const XML_Char **attr) {
	((PepXMLParser3*) data)->startElement(el, attr);
}

static void PepXMLParser3_endElementCallback(void *data, const XML_Char *el){
	((PepXMLParser3*) data)->endElement(el);
}

static void PepXMLParser3_charactersCallback(void *data, const XML_Char *s, int len){
	((PepXMLParser3*) data)->characters(s, len);
}

void PepXMLParser3::characters(const XML_Char *s, int len) {
	//nothing is held here.
}

void PepXMLParser3::endElement(const XML_Char *el) {

	string s;
  char str[32];
	size_t i,j;

  if (isElement("analysis_result", el)){
    /*
    CPepXMLAnalysisResult ar;
    psm.analysisResults->push_back(ar);
    switch(analysisState){
    case asPeptideProphet:
      psm.analysisResults->back().setAnalysisResult(prophetResult,PeptideProphet);
      break;
    case asInterProphet:
      psm.analysisResults->back().setAnalysisResult(prophetResult, InterProphet);
      break;
    case asPTMProphet:
      psm.analysisResults->back().setAnalysisResult(ptmProphetResult, PTMProphet);
      break;
    default:
      cout << "WARNING: Unknown ending analysis_result" << endl;
      break;
    }
    */
    analysisState=asNone;

  } else if (isElement("analysis_summary", el)){
    vAnalysis.push_back(analysis);
    analysis.clear();

  } else if (isElement("interprophet_result", el)) {
    //CPepXMLAnalysisResult ar;
    //psm.analysisResults->push_back(ar);
    psm.analysisResults->back().setAnalysisResult(prophetResult, InterProphet);
    
  } else if (isElement("linked_peptide", el)){

    if(psm.peptide==NULL) {
      psm.peptide = new CPepXMLPeptide;
      *psm.peptide=peptide;
    } else {
      if(psm.xlPeptide==NULL) psm.xlPeptide = new CPepXMLPeptide;
      *psm.xlPeptide=peptide;
    }
    peptide.clear();

  } else if(isElement("modification_info",el)) {
    
    //Build out the modified peptide
    peptide.modifiedPeptide.clear();
    for (j = 0; j<peptide.sizeMod(); j++){
      if (peptide.mods->at(j).pos == -1){ //n-terminal mods
        sprintf(str, "n[%.0lf]", vMods[peptide.mods->at(j).index].massSearch);
        peptide.modifiedPeptide += str;
      }
    }
    for(i=0;i<peptide.peptide.size();i++){
      peptide.modifiedPeptide+=peptide.peptide[i];
      for(j=0;j<peptide.sizeMod();j++){
        if(peptide.mods->at(j).pos==i){
          sprintf(str,"[%.0lf]",vMods[peptide.mods->at(j).index].massSearch);
          peptide.modifiedPeptide+=str;
        }
      }
    }

  } else if (isElement("peptideprophet_result", el)) {
    //CPepXMLAnalysisResult ar;
    //psm.analysisResults->push_back(ar);
    psm.analysisResults->back().setAnalysisResult(prophetResult, PeptideProphet);

  } else if (isElement("ptmprophet_result", el)) {
    analysisState = asPTMProphet;
    psm.analysisResults->back().setAnalysisResult(ptmProphetResult, PTMProphet);

  } else if(isElement("search_hit", el))	{

    if(psm.xlType<2) {
      if(psm.peptide==NULL) psm.peptide = new CPepXMLPeptide;
      *psm.peptide=peptide;
    }
    if(psm.probability>=probabilityFilter && psm.iProphetProbability>=iProbabilityFilter){
      spec.psms->push_back(psm);
    }
    peptide.clear();
    psm.clear();

  } else if (isElement("search_summary", el)) {

    currentSearchSummary = (char)vSearch.size();
    vSearch.push_back(search);
    search.clear();

	} else if(isElement("spectrum_query", el))	{

    if(spec.psms->size()>0) vSpectra.push_back(spec);
    spec.clear();
    
	}

}


void PepXMLParser3::startElement(const XML_Char *el, const XML_Char **attr){

  char          c;
  double        d1,d2;
  string        s;
	size_t        u;

  if (isElement("alternative_protein", el)){
    s = getAttrValue("protein", attr);
    peptide.proteins->push_back(findProtein(s,currentDBID));

  } else if (isElement("aminoacid_modification",el)){
    s = getAttrValue("aminoacid",attr);
    c = s[0];
    s = getAttrValue("mass",attr);
    d1 = atof(&s[0]);
    s = getAttrValue("massdiff",attr);
    d2 = atof(&s[0]);
    addMod(c,d1,d2);

    PepXMLSearchMod m;
    m.aa=c;
    m.mass=d2;
    s = getAttrValue("variable", attr);
    if (s[0] == 'N') m.fixed=true;
    else m.fixed=false;
    m.proteinTerm=false;
    search.mods->push_back(m);

  } else if (isElement("analysis_result", el)){
    CPepXMLAnalysisResult ar;
    psm.analysisResults->push_back(ar);

    s = getAttrValue("analysis", attr);
    if(s.compare("peptideprophet")==0){
      analysisState=asPeptideProphet;
      prophetResult = new CPepXMLProphetResult;
    } else if (s.compare("interprophet")==0){
      analysisState=asInterProphet;
      prophetResult = new CPepXMLProphetResult;
    } else if (s.compare("ptmprophet")==0){
      analysisState=asPTMProphet;
      //ptmProphetResult = new CPepXMLPTMProphetResult;
    } else {
      cout << "WARNING: Unknown analysis_result: " << s << endl;
    }

  } else if (isElement("analysis_summary", el)){
    analysis.analysis = getAttrValue("analysis", attr);

  } else if (isElement("error_point",el)){
    PepXMLError e;
    s = getAttrValue("error",attr);
    e.error=atof(&s[0]);
    s = getAttrValue("min_prob",attr);
    e.prob=atof(&s[0]);
    vError.push_back(e);

	} else if (isElement("interprophet_result",el)){
		s = getAttrValue("probability", attr);
		psm.iProphetProbability=atof(&s[0]);

  } else if(isElement("interprophet_summary",el)){
    bPepProphet = true;  //because you can't iprophet without peptideprophet
    bIProphet=true;
    analysis.version = getAttrValue("version", attr);

  } else if (isElement("linked_peptide",el)){

    s = getAttrValue("calc_neutral_pep_mass", attr);
    peptide.calcPepNeutMass=atof(&s[0]);
		peptide.peptide = getAttrValue("peptide", attr);
    s = getAttrValue("protein", attr);
    peptide.proteins->push_back(findProtein(s,currentDBID));
    s = getAttrValue("peptide_prev_aa", attr);
    if(s.length()>0) peptide.prevAA=s[0];
    else peptide.prevAA='-';
    s = getAttrValue("peptide_next_aa", attr);
    if(s.length()>0) peptide.nextAA=s[0];
    else peptide.nextAA='-';
    s = getAttrValue("complement_mass", attr);
    peptide.complementMass=atof(&s[0]);
    s = getAttrValue("designation", attr);
    if(s.compare("alpha")==0) peptide.label=0;
    else peptide.label=1;

  } else if( isElement("modification_info",el)){

    PepXMLPepMod pm;
    s = getAttrValue("mod_nterm_mass",attr);  
    if(s.size()>0){
      d1 = atof(&s[0]);
      pm.index = findMod('n',d1);
      pm.pos = -1;
      peptide.addMod(pm);
    }
    s = getAttrValue("mod_cterm_mass",attr);  
    if(s.size()>0){
      d1 = atof(&s[0]);
      pm.index = findMod('c',d1);
      pm.pos = (char)peptide.peptide.size()-1; //this position might need to be marked after the sequence.
      peptide.addMod(pm);
    }

	} else if (isElement("mod_aminoacid_mass",el)){

    PepXMLPepMod pm;
    s = getAttrValue("mass",attr);
		d1 = atof(&s[0]);
    s = getAttrValue("position",attr);
    pm.pos = (char)atoi(&s[0])-1;
    c = peptide.peptide[pm.pos];
    pm.index=findMod(c,d1);
    peptide.addMod(pm);

  } else if (isElement("mod_aminoacid_probability", el)){

    PepXMLPTMMod pm;
    s = getAttrValue("position", attr);
    pm.position = (char)atoi(s.c_str());
    s = getAttrValue("pval", attr);
    pm.pVal=atof(s.c_str());
    s = getAttrValue("probability", attr);
    pm.probability = atof(s.c_str());
    s = getAttrValue("oscore", attr);
    pm.oscore = atof(s.c_str());
    s = getAttrValue("mscore", attr);
    pm.mscore = atof(s.c_str());
    s = getAttrValue("cterm_score", attr);
    pm.ctermscore = atof(s.c_str());
    s = getAttrValue("nterm_score", attr);
    pm.ntermscore = atof(s.c_str());
    ptmProphetResult->mods->push_back(pm);

  } else if (isElement("msms_run_summary",el)){

    s = getAttrValue("base_name",attr);
    s += getAttrValue("raw_data",attr);
    for(u=0;u<vFiles.size();u++){
      if(s.compare(vFiles[u])==0) break;
    }
    currentFileID=(int)u;
    if(u==vFiles.size()) vFiles.push_back(s);

  } else if (isElement("parameter", el)){

    PepXMLParam par;
    PepXMLPepScore ps;
    if(analysisState==asNone){
      par.name = getAttrValue("name", attr);
      par.value = getAttrValue("value", attr);
    } else {
      s = getAttrValue("name", attr);
      ps.index = findScore(s);
      s = getAttrValue("value", attr);
      ps.value = atof(&s[0]);
    }
    switch(analysisState){
    case asPeptideProphet:
    case asInterProphet:
      prophetResult->parameters->push_back(ps);
      break;
    case asPTMProphet:
      psm.analysisResults->back().parameters->push_back(ps);
      break;
    case asPTMProphetResult:
      ptmProphetResult->parameters->push_back(ps);
      break;
    default:
      search.params->push_back(par);
      break;
    }

	} else if (isElement("peptideprophet_result",el)){

		s = getAttrValue("probability", attr);
		psm.probability=atof(&s[0]);
    prophetResult->probability=psm.probability;

  } else if (isElement("peptideprophet_summary", el)){
    bPepProphet=true;

  } else if (isElement("ptmprophet_result", el)){
    if(psm.analysisResults->back().analysisResult!=NULL){ //if analysis result was already placed here, add a new one
      CPepXMLAnalysisResult ar;
      for(size_t xp=0;xp<psm.analysisResults->back().parameters->size();xp++){
        ar.parameters->push_back(psm.analysisResults->back().parameters->at(xp));
      }
      psm.analysisResults->push_back(ar);
    }
      
    analysisState = asPTMProphetResult;
    ptmProphetResult = new CPepXMLPTMProphetResult;
    s = getAttrValue("prior", attr);
    ptmProphetResult->prior=atof(s.c_str());
    ptmProphetResult->ptm = getAttrValue("ptm", attr);
    ptmProphetResult->ptm_peptide = getAttrValue("ptm_peptide", attr);

  } else if (isElement("search_database",el)){

    s = getAttrValue("local_path", attr);
    for (u = 0; u<vDBs.size(); u++){
      if (s.compare(vDBs[u]) == 0) break;
    }
    currentDBID = (char)u;
    if (u == vDBs.size()) vDBs.push_back(s);
    search.DBindex = currentDBID;

	} else if (isElement("search_hit",el)) {

    string desc;
		s = getAttrValue("hit_rank", attr);
		psm.rank = atoi(&s[0]);
    s = getAttrValue("xlink_type",attr);
    if(s.length()>0 && s.compare("xl")==0) {
      psm.xlType=2;
      s = getAttrValue("calc_neutral_pep_mass", attr);
      psm.calcPSMNeutMass=atof(&s[0]);
    } else { //only grab peptide info if non-linked or loop-linked
      if(s.length()>0 && s.compare("loop")==0) psm.xlType=1;
      else psm.xlType=0;
		  s = getAttrValue("calc_neutral_pep_mass", attr);
      peptide.calcPepNeutMass=atof(&s[0]);
      psm.calcPSMNeutMass=peptide.calcPepNeutMass;
		  peptide.peptide = getAttrValue("peptide", attr);
      desc = getAttrValue("protein_descr",attr);
      s = getAttrValue("protein", attr);
      peptide.proteins->push_back(findProtein(s,currentDBID,desc));
      s = getAttrValue("peptide_prev_aa", attr);
      if(s.length()>0) peptide.prevAA=s[0];
      else peptide.prevAA='-';
      s = getAttrValue("peptide_next_aa", attr);
      if(s.length()>0) peptide.nextAA=s[0];
      else peptide.nextAA='-';
    }

	} else if (isElement("search_score",el)) {

    PepXMLPepScore ps;
		s = getAttrValue("name", attr);
    ps.index = findScore(s);
    s = getAttrValue("value", attr);
    ps.value = atof(&s[0]);
    psm.psmScores->push_back(ps);

  } else if (isElement("search_summary",el)) {

    search.alg = getAttrValue("search_engine",attr);
    search.version = getAttrValue("search_engine_version",attr);

	} else if (isElement("spectrum_query",el)) {

		s = getAttrValue("start_scan", attr);
    spec.scanStart = atoi(&s[0]);
    spec.scanNumber = spec.scanStart;  //This might not always be true
    s = getAttrValue("end_scan", attr);
    spec.scanEnd = atoi(&s[0]);
    s = getAttrValue("precursor_neutral_mass", attr);
    spec.precursorNeutMass = atof(&s[0]);
		s = getAttrValue("assumed_charge", attr);
		spec.charge = atoi(&s[0]);
		s = getAttrValue("retention_time_sec", attr);
    spec.rTimeSec = (float)atof(&s[0]);
    spec.ID = getAttrValue("spectrum", attr);
    spec.nativeID = getAttrValue("spectrumNativeID",attr);
    spec.fileID = currentFileID;
    spec.searchID = currentSearchSummary;
    //s = getAttrValue("index", attr);
    //spec.index = atoi(&s[0]);

	} else if (isElement("terminal_modification",el)){

    s = getAttrValue("terminus",attr);
    if(s[0]=='N' || s[0]=='n') c='n';
    else c='c';
    s = getAttrValue("mass",attr);
    d1 = atof(&s[0]);
    s = getAttrValue("massdiff",attr);
    d2 = atof(&s[0]);
    addMod(c,d1,d2);

    PepXMLSearchMod m;
    m.aa = c;
    m.mass = d2;
    s = getAttrValue("variable", attr);
    if (s[0] == 'N') m.fixed = true;
    else m.fixed = false;
    s = getAttrValue("protein_terminus", attr);
    if (s[0] == 'N') m.proteinTerm = false;
    else m.proteinTerm = true;
    search.mods->push_back(m);

  } else if (isElement("xlink",el)){

    s = getAttrValue("mass", attr);
    d1 = atof(&s[0]);
    s = getAttrValue("identifier", attr);
    psm.xlIndex = findXL(s,d1);

  } else if (isElement("xlink_score",el)) {

    PepXMLPepScore ps;
		s = getAttrValue("name", attr);
    ps.index = findScore(s);
    s = getAttrValue("value", attr);
    ps.value = atof(&s[0]);
    peptide.xlScores->push_back(ps);

  }
}

PepXMLParser3::PepXMLParser3(){
  analysisState = asNone;
  bIProphet   = false;
  bPepProphet = false;
	parser      = XML_ParserCreate(NULL);
	XML_SetUserData(parser, this);
	XML_SetElementHandler(parser, PepXMLParser3_startElementCallback, PepXMLParser3_endElementCallback);
	XML_SetCharacterDataHandler(parser, PepXMLParser3_charactersCallback);
}

PepXMLParser3::~PepXMLParser3(){
	XML_ParserFree(parser);
}

CPepXMLSpectrum& PepXMLParser3::operator[ ](const size_t& i){
  return vSpectra[i];
}

bool PepXMLParser3::addMod(char aa, double mass, double massDiff){
  size_t i;
  PepXMLMod mod;
  
  mod.aa=aa;
  mod.massSearch=mass;
  mod.massDiff=massDiff;
  mod.massDiffStd=0;
  mod.massStd=0;
  mod.label.clear();

  //Look up actual modification mass: these need to be adjusted!!
  switch(aa){
    case 'C':
      if(fabs(mass-160.0307)<0.001) {
        mod.massDiffStd=57.021464;
        mod.massStd=calcMonoMass("C",false)+57.021464;
        mod.label="Carbamidomethyl_C";
      }
      if(fabs(mass-143.0042)<0.001) {
        mod.massDiffStd=39.9949141;
        mod.massStd=calcMonoMass("C",false)+39.9949141;
        mod.label="Carbamidomethyl_C-NH3";
      }
      break;
    case 'E':
      if(fabs(mass-111.0320)<0.001) {
        mod.massDiffStd=-18.0105633;
        mod.massStd=calcMonoMass("E",false)-18.0105633;
        mod.label="Pyroglutamate";
      }
      break;
    case 'M':
      if(fabs(mass-147.0353937)<0.001 || fabs(mass-147.039885)<0.001) {
        mod.massDiffStd=15.9949141;
        mod.massStd=calcMonoMass("M",false)+15.9949141;
        mod.label="Oxidized_M";
      }
      break;
    case 'K':
      if (fabs(mass - 272.197023)<0.001) {
        mod.massDiffStd = 144.102060;
        mod.massStd = calcMonoMass("K", false) + 144.102060;
        mod.label = "iTraq4";
      } else if (fabs(mass - 284.173556)<0.001) {
        mod.massDiffStd = 156.078644;
        mod.massStd = calcMonoMass("K", false) + 156.078644;
        mod.label = "Water_Quenched_DSS/BS3";
      } else if (fabs(mass - 283.189556)<0.001) {
        mod.massDiffStd = 155.094629;
        mod.massStd = calcMonoMass("K", false) + 155.094629;
        mod.label = "Ammonium_Quenched_DSS/BS3";
      }
      break;
    case 'n':
      if(fabs(mass-43.0184)<0.001) {
        mod.massDiffStd=42.0105633;
        mod.massStd=0;
        mod.label="Acetylation";
      } else if (fabs(mass - 145.109885)<0.001) {
        mod.massDiffStd = 144.102060;
        mod.massStd = 0;
        mod.label = "iTraq4";
      } else if (fabs(mass - 157.0864)<0.001) {
        mod.massDiffStd = 156.078644;
        mod.massStd = 0;
        mod.label = "Water_Quenched_DSS/BS3";
      } else if (fabs(mass - 156.1024)<0.001) {
        mod.massDiffStd = 155.094629;
        mod.massStd = 0;
        mod.label = "Ammonium_Quenched_DSS/BS3";
      }
      break;
    case 'N':
      if(fabs(mass-115.0269)<0.001) {
        mod.massDiffStd=-0.984016;
        mod.massStd=calcMonoMass("N",false)-0.984016;
        mod.label="Deamidation_N";
      }
      break;
    case 'Q':
      if(fabs(mass-129.0429)<0.001) {
        mod.massDiffStd=-0.984016;
        mod.massStd=calcMonoMass("Q",false)-0.984016;
        mod.label="Deamidation_Q";
      }
      if(fabs(mass-111.0321)<0.001) {
        mod.massDiffStd=-17.026549;
        mod.massStd=calcMonoMass("Q",false)-17.026549;
        mod.label="Pyrrolidone";
      }
      break;
    case 'R':
      if(fabs(mass-157.0851)<0.001) {
        mod.massDiffStd=-0.984016;
        mod.massStd=calcMonoMass("R",false)-0.984016;
        mod.label="Citrulline";
      }
      break;
    case 'S':
      if(fabs(mass-166.9984)<0.001) {
        mod.massDiffStd=79.9663289;
        mod.massStd=calcMonoMass("S",false)+79.9663289;
        mod.label="Phosphorylation_S";
      }
      break;
    case 'T':
      if(fabs(mass-181.0140)<0.001) {
        mod.massDiffStd=79.9663289;
        mod.massStd=calcMonoMass("T",false)+79.9663289;
        mod.label="Phosphorylation_T";
      }
      break;
    case 'W':
      if (fabs(mass - 202.074228)<0.001) {
        mod.massDiffStd = 15.9949141;
        mod.massStd = calcMonoMass("W", false) + 15.9949141;
        mod.label = "Oxidized_W";
      } else {
        mod.massDiffStd = massDiff;
        mod.massStd = calcMonoMass("W", false) + massDiff;
        mod.label = "Unknown_W";
      }
    case 'Y':
      if(fabs(mass-243.0297)<0.001) {
        mod.massDiffStd=79.9663289;
        mod.massStd=calcMonoMass("Y",false)+79.9663289;
        mod.label="Phosphorylation_Y";
      }
      break;
    default:
      break;
  }

  if(mod.label.size()==0){
    cout << "Unknown xx modification: " << aa << " = " << mass << "\t" << fabs(mass - 157.0864) << endl;
  }

  for(i=0;i<vMods.size();i++){
    if(vMods[i].aa==mod.aa && fabs(vMods[i].massSearch-mod.massSearch)<0.0000001) break;
  }
  if(i==vMods.size()){
    vMods.push_back(mod);
    return true;
  }
  return false;

}

double PepXMLParser3::calcMonoMass(const char *seq, bool water){

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

char PepXMLParser3::findMod(char aa, double mass){
  size_t i;
  for(i=0;i<vMods.size();i++){
    if(vMods[i].aa==aa && fabs(vMods[i].massSearch-mass)<0.0000001) return (char)i;
  }

  //didn't match mod, so add it
  if(!addMod(aa,mass)){
    cout << "WARNING: possible duplicate mod: " << aa << ", " << mass << endl;
  }
  return (char)i;
}

size_t PepXMLParser3::findProtein(string& s, char index, string desc){

  size_t sz = vProtTable.size() - vProtTable.size()%100;
  size_t lower = 0;
  size_t mid = sz / 2;
  size_t upper = sz;
  int i;

  if (sz>0){
    i = vProtTable[mid].name.compare(s);
    while (i != 0){
      if (lower >= upper) break;
      if (i>0){
        if (mid == 0) break;
        upper = mid - 1;
        mid = (lower + upper) / 2;
      } else {
        lower = mid + 1;
        mid = (lower + upper) / 2;
      }
      if (mid == sz) break;
      i = vProtTable[mid].name.compare(s);
    }

    //match by protein name, so check index
    if (i == 0){
      if (vProtTable[mid].DBIndex == index) return vProtTable[mid].protIndex;
      lower = mid-1;
      while (lower < mid && vProtTable[lower].name.compare(s) == 0){
        if (vProtTable[lower].DBIndex == index) return vProtTable[lower].protIndex;
        lower--;
      }
      upper=mid+1;
      while (upper < sz && vProtTable[upper].name.compare(s) == 0){
        if (vProtTable[upper].DBIndex == index) return vProtTable[upper].protIndex;
        upper++;
      }
    }
  }

  //check unsorted proteins
  for (mid = sz; mid < vProtTable.size(); mid++){
    if (vProtTable[mid].name.compare(s) == 0 && vProtTable[mid].DBIndex == index) return vProtTable[mid].protIndex;
  }

  //protein was not found, so add it;
  PepXMLProtein pr;
  pr.name=s;
  pr.DBindex=index;
  pr.description=desc;
  vProteins.push_back(pr);

  PepXMLProtTable pt;
  pt.name=s;
  pt.DBIndex=index;
  pt.protIndex = vProteins.size()-1;
  vProtTable.push_back(pt);

  //sort if buffer is full
  if (vProtTable.size() % 100 == 0) {
    sort(vProtTable.begin(),vProtTable.end(),compareProt);
  }

  return pt.protIndex;
}

char PepXMLParser3::findScore(string& s){
  size_t i;
  for(i=0;i<vScores.size();i++){
    if(vScores[i].compare(s)==0) return (char)i;
  }
  vScores.push_back(s);
  return (char)i;
}

char PepXMLParser3::findXL(string& s, double mass){
  size_t i;
  for(i=0;i<vXL.size();i++){
    if(vXL[i].ID.compare(s)==0 && fabs(vXL[i].mass-mass)<0.0000001) return (char)i;
  }
  PepXMLXL p;
  p.ID=s;
  p.mass=mass;
  vXL.push_back(p);
  return (char)i;
}

CPepXMLAnalysis PepXMLParser3::getAnalysis(size_t index){
  return vAnalysis[index];
}

size_t PepXMLParser3::getAnalysisCount(){
  return vAnalysis.size();
}

string PepXMLParser3::getFile(int index){
  return vFiles[vSpectra[index].fileID];
}

int PepXMLParser3::getFileCount(){
  return (int)vFiles.size();
}

void PepXMLParser3::getFileFromList(int index, char* str){
  strcpy(str,&vFiles[index][0]);
}

PepXMLXL PepXMLParser3::getLinker(size_t index, char rank){
  size_t i;
  PepXMLXL xl;
  xl.ID.clear();
  xl.mass=0;
  if(index>=vSpectra.size()) return xl;
  for(i=0;i<vSpectra[index].size();i++){
    if(vSpectra[index][i].rank==rank) break;
  }
  if(i==vSpectra[index].size()) return xl;
  if(vSpectra[index][i].xlType==0) return xl;
  return vXL[vSpectra[index][i].xlIndex];
}

PepXMLXL PepXMLParser3::getLinkerFromPos(size_t index, size_t nth){
  PepXMLXL xl;
  xl.ID.clear();
  xl.mass = 0;
  if (index >= vSpectra.size()) return xl;
  if (nth>=vSpectra[index].size()) return xl;
  if (vSpectra[index][nth].xlType == 0) return xl;
  return vXL[vSpectra[index][nth].xlIndex];
}

bool PepXMLParser3::getLinkSites(size_t index, char& a, char& b, char rank){
  size_t i,j,k;
  if(index>=vSpectra.size()) return false;
  for(i=0;i<vSpectra[index].size();i++){
    if(vSpectra[index][i].rank==rank) break;
  }
  if(i==vSpectra[index].size()) return false;
  switch(vSpectra[index][i].xlType) {
    case 1:
      k=0;
      for(j=0;j<vSpectra[index][i].peptide->xlScores->size();j++){
        if(getScoreLabel(vSpectra[index][i].peptide->xlScores->at(j).index).compare("link")==0){
          if(k==0) {
            a=(char)((int)(vSpectra[index][i].peptide->xlScores->at(j).value+0.5));
            k++;
          } else {
            b=(char)((int)(vSpectra[index][i].peptide->xlScores->at(j).value+0.5));
            break;
          }
        }
      }
      break;
    case 2:
      for(j=0;j<vSpectra[index][i].peptide->xlScores->size();j++){
        if(getScoreLabel(vSpectra[index][i].peptide->xlScores->at(j).index).compare("link")==0){
          a=(char)(vSpectra[index][i].peptide->xlScores->at(j).value+0.5);
          break;
        }
      }
      for(j=0;j<vSpectra[index][i].xlPeptide->xlScores->size();j++){
        if(getScoreLabel(vSpectra[index][i].xlPeptide->xlScores->at(j).index).compare("link")==0){
          b=(char)(vSpectra[index][i].xlPeptide->xlScores->at(j).value+0.5);
          break;
        }
      }
      break;
    default:
      return false;
  }
  return true;
}

bool PepXMLParser3::getLinkSitesFromPos(size_t index, char& a, char& b, size_t nth){
  size_t j, k;
  if (index >= vSpectra.size()) return false;
  if (nth>=vSpectra[index].size()) return false;
  switch (vSpectra[index][nth].xlType) {
  case 1:
    k = 0;
    for (j = 0; j<vSpectra[index][nth].peptide->xlScores->size(); j++){
      if (getScoreLabel(vSpectra[index][nth].peptide->xlScores->at(j).index).compare("link") == 0){
        if (k == 0) {
          a = (char)((int)(vSpectra[index][nth].peptide->xlScores->at(j).value + 0.5));
          k++;
        } else {
          b = (char)((int)(vSpectra[index][nth].peptide->xlScores->at(j).value + 0.5));
          break;
        }
      }
    }
    break;
  case 2:
    for (j = 0; j<vSpectra[index][nth].peptide->xlScores->size(); j++){
      if (getScoreLabel(vSpectra[index][nth].peptide->xlScores->at(j).index).compare("link") == 0){
        a = (char)(vSpectra[index][nth].peptide->xlScores->at(j).value + 0.5);
        break;
      }
    }
    for (j = 0; j<vSpectra[index][nth].xlPeptide->xlScores->size(); j++){
      if (getScoreLabel(vSpectra[index][nth].xlPeptide->xlScores->at(j).index).compare("link") == 0){
        b = (char)(vSpectra[index][nth].xlPeptide->xlScores->at(j).value + 0.5);
        break;
      }
    }
    break;
  default:
    return false;
  }
  return true;
}

char PepXMLParser3::getLinkType(size_t index, char rank){
  size_t i;
  if(index>=vSpectra.size()) return -1;
  for(i=0;i<vSpectra[index].size();i++){
    if(vSpectra[index][i].rank==rank) break;
  }
  if(i==vSpectra[index].size()) return -1;
  return vSpectra[index][i].xlType;
}

string PepXMLParser3::getPeptide(size_t index, bool mod, char rank, bool link){
  size_t i;
  if(index>=vSpectra.size()) return "";
  for(i=0;i<vSpectra[index].size();i++){
    if(vSpectra[index][i].rank==rank) break;
  }
  if(i==vSpectra[index].size()) return "";
  if(link){
    if(vSpectra[index][i].xlType<2) return "";
    if(mod){
      if(vSpectra[index][i].xlPeptide->modifiedPeptide.size()>0) return vSpectra[index][i].xlPeptide->modifiedPeptide;
    }
    return vSpectra[index][i].xlPeptide->peptide;
  }
  if(mod){
    if(vSpectra[index][i].peptide->modifiedPeptide.size()>0) return vSpectra[index][i].peptide->modifiedPeptide;
  }
  return vSpectra[index][i].peptide->peptide;
}

string PepXMLParser3::getPeptideFromPos(size_t index, bool mod, size_t nth, bool link){
  if (index >= vSpectra.size()) return "";
  if (nth>=vSpectra.size()) return "";
  if (link){
    if (vSpectra[index][nth].xlType<2) return "";
    if (mod){
      if (vSpectra[index][nth].xlPeptide->modifiedPeptide.size()>0) return vSpectra[index][nth].xlPeptide->modifiedPeptide;
    }
    return vSpectra[index][nth].xlPeptide->peptide;
  }
  if (mod){
    if (vSpectra[index][nth].peptide->modifiedPeptide.size()>0) return vSpectra[index][nth].peptide->modifiedPeptide;
  }
  return vSpectra[index][nth].peptide->peptide;
}

PepXMLMod PepXMLParser3::getPeptideMod(size_t pepIndex, size_t modIndex, char rank, bool link){
  size_t i;
  PepXMLMod p;
  p.aa=0;
  p.label.clear();
  p.massDiff=0;
  p.massDiffStd=0;
  p.massSearch=0;
  p.massStd=0;

  if(pepIndex>=vSpectra.size()) return p;
  for(i=0;i<vSpectra[pepIndex].size();i++){
    if(vSpectra[pepIndex][i].rank==rank) break;
  }
  if(i==vSpectra[pepIndex].size()) return p;
  if(link){
    if(vSpectra[pepIndex][i].xlType<2) return p;
    if(modIndex>=vSpectra[pepIndex][i].xlPeptide->mods->size()) return p;
    p=vMods[vSpectra[pepIndex][i].xlPeptide->mods->at(modIndex).index];
    p.aa=vSpectra[pepIndex][i].xlPeptide->mods->at(modIndex).pos;
    return p;
  }
  if(modIndex>=vSpectra[pepIndex][i].peptide->mods->size()) return p;
  p=vMods[vSpectra[pepIndex][i].peptide->mods->at(modIndex).index];
  p.aa=vSpectra[pepIndex][i].peptide->mods->at(modIndex).pos;
  return p;
}

PepXMLMod PepXMLParser3::getPeptideModFromPos(size_t pepIndex, size_t modIndex, size_t nth, bool link){
  PepXMLMod p;
  p.aa = 0;
  p.label.clear();
  p.massDiff = 0;
  p.massDiffStd = 0;
  p.massSearch = 0;
  p.massStd = 0;

  if (pepIndex >= vSpectra.size()) return p;
  if (nth>=vSpectra[pepIndex].size()) return p;
  if (link){
    if (vSpectra[pepIndex][nth].xlType<2) return p;
    if (modIndex >= vSpectra[pepIndex][nth].xlPeptide->mods->size()) return p;
    p = vMods[vSpectra[pepIndex][nth].xlPeptide->mods->at(modIndex).index];
    p.aa = vSpectra[pepIndex][nth].xlPeptide->mods->at(modIndex).pos;
    return p;
  }
  if (modIndex >= vSpectra[pepIndex][nth].peptide->mods->size()) return p;
  p = vMods[vSpectra[pepIndex][nth].peptide->mods->at(modIndex).index];
  p.aa = vSpectra[pepIndex][nth].peptide->mods->at(modIndex).pos;
  return p;
}

size_t PepXMLParser3::getPeptideModCount(size_t index, char rank, bool link){
  size_t i;
  if(index>=vSpectra.size()) return 0;
  for(i=0;i<vSpectra[index].size();i++){
    if(vSpectra[index][i].rank==rank) break;
  }
  if(i==vSpectra[index].size()) return 0;
  if(link){
    if(vSpectra[index][i].xlType<2) return 0;
    return vSpectra[index][i].xlPeptide->mods->size();
  }
  return vSpectra[index][i].peptide->mods->size();
}

size_t PepXMLParser3::getPeptideModCountFromPos(size_t index, size_t nth, bool link){
  if (index >= vSpectra.size()) return 0;
  if (nth >= vSpectra[index].size()) return 0;
  if (link){
    if (vSpectra[index][nth].xlType<2) return 0;
    if (vSpectra[index][nth].xlPeptide->mods==NULL) return 0;
    return vSpectra[index][nth].xlPeptide->mods->size();
  }
  if (vSpectra[index][nth].peptide->mods==NULL) return 0;
  return vSpectra[index][nth].peptide->mods->size();
}

//bool PepXMLParser::getIprophet(){
//  return bIprophet;
//}

double PepXMLParser3::getProbability(double err){
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

string PepXMLParser3::getProtein(size_t index, size_t protIndex, char rank, bool link){
  size_t i;
  if(index>=vSpectra.size()) return "";
  for(i=0;i<vSpectra[index].size();i++){
    if(vSpectra[index][i].rank==rank) break;
  }
  if(i==vSpectra[index].size()) return "";
  if(link){
    if(vSpectra[index][i].xlType<2) return "";
    if(protIndex>=vSpectra[index][i].xlPeptide->proteins->size()) return "";
    return vProteins[vSpectra[index][i].xlPeptide->proteins->at(protIndex)].name;
  }
  if(protIndex>=vSpectra[index][i].peptide->proteins->size()) return "";
  return vProteins[vSpectra[index][i].peptide->proteins->at(protIndex)].name;
}

string PepXMLParser3::getProteinFromList(int protIndex){
  if(protIndex>vProteins.size()) return "";
  if(protIndex<0) return "";
  return vProteins[(size_t)protIndex].name;
}

string PepXMLParser3::getProteinFromPos(size_t index, size_t protIndex, size_t nth, bool link){
  if (index >= vSpectra.size()) return "";
  if (nth >= vSpectra[index].size()) return "";
  if (link){
    if (vSpectra[index][nth].xlType<2) return "";
    if (protIndex >= vSpectra[index][nth].xlPeptide->proteins->size()) return "";
    return vProteins[vSpectra[index][nth].xlPeptide->proteins->at(protIndex)].name;
  }
  if (protIndex >= vSpectra[index][nth].peptide->proteins->size()) return "";
  return vProteins[vSpectra[index][nth].peptide->proteins->at(protIndex)].name;
}

string PepXMLParser3::getProteinDB(size_t index, size_t protIndex, char rank, bool link){
  size_t i;
  if (index >= vSpectra.size()) return "";
  for (i = 0; i<vSpectra[index].size(); i++){
    if (vSpectra[index][i].rank == rank) break;
  }
  if (i == vSpectra[index].size()) return "";
  if (link){
    if (vSpectra[index][i].xlType<2) return "";
    if (protIndex >= vSpectra[index][i].xlPeptide->proteins->size()) return "";
    return vDBs[(size_t)vProteins[vSpectra[index][i].xlPeptide->proteins->at(protIndex)].DBindex];
  }
  if (protIndex >= vSpectra[index][i].peptide->proteins->size()) return "";
  return vDBs[(size_t)vProteins[vSpectra[index][i].peptide->proteins->at(protIndex)].DBindex];
}

string PepXMLParser3::getProteinDesc(size_t index, size_t protIndex, char rank, bool link){
  size_t i;
  if (index >= vSpectra.size()) return "";
  for (i = 0; i<vSpectra[index].size(); i++){
    if (vSpectra[index][i].rank == rank) break;
  }
  if (i == vSpectra[index].size()) return "";
  if (link){
    if (vSpectra[index][i].xlType<2) return "";
    if (protIndex >= vSpectra[index][i].xlPeptide->proteins->size()) return "";
    return vProteins[vSpectra[index][i].xlPeptide->proteins->at(protIndex)].description;
  }
  if (protIndex >= vSpectra[index][i].peptide->proteins->size()) return "";
  return vProteins[vSpectra[index][i].peptide->proteins->at(protIndex)].description;
}

string PepXMLParser3::getProteinDescFromPos(size_t index, size_t protIndex, size_t nth, bool link){
  if (index >= vSpectra.size()) return "";
  if (nth >= vSpectra[index].size()) return "";
  if (link){
    if (vSpectra[index][nth].xlType<2) return "";
    if (protIndex >= vSpectra[index][nth].xlPeptide->proteins->size()) return "";
    return vProteins[vSpectra[index][nth].xlPeptide->proteins->at(protIndex)].description;
  }
  if (protIndex >= vSpectra[index][nth].peptide->proteins->size()) return "";
  return vProteins[vSpectra[index][nth].peptide->proteins->at(protIndex)].description;
}

CPepXMLPSM PepXMLParser3::getPSM(size_t index, char rank, bool link){
  size_t i;
  if (index >= vSpectra.size()) {
    CPepXMLPSM p;
    return p;
  }
  for (i = 0; i<vSpectra[index].size(); i++){
    if (vSpectra[index][i].rank == rank) break;
  }
  if (i == vSpectra[index].size()) {
    CPepXMLPSM p;
    return p;
  }
  if (link){
    if (vSpectra[index][i].xlType<2) {
      CPepXMLPSM p;
      return p;
    }
    return vSpectra[index][i];
  }
  return vSpectra[index][i];
}

CPepXMLPSM PepXMLParser3::getPSMFromPos(size_t index, size_t nth, bool link){
  if (index >= vSpectra.size()) {
    CPepXMLPSM p;
    return p;
  }
  if(nth>=vSpectra[index].size()){
    CPepXMLPSM p;
    return p;
  }
  if (link){
    if (vSpectra[index][nth].xlType<2) {
      CPepXMLPSM p;
      return p;
    }
    return vSpectra[index][nth];
  }
  return vSpectra[index][nth];
}


char PepXMLParser3::getScoreIndex(string s){
  size_t i;
  for (i = 0; i<vScores.size(); i++){
    if (vScores[i].compare(s) == 0) return (char)i;
  }
  return -1;
}

string PepXMLParser3::getScoreLabel(char scoreIndex){
  size_t i=(size_t)scoreIndex;
  if(i>=vScores.size()) return "";
  return vScores[i];
}

CPepXMLSearch PepXMLParser3::getSearchParams(size_t pepIndex){
  return vSearch[vSpectra[pepIndex].searchID];
}

bool PepXMLParser3::hasIProphet(){
  return bIProphet;
}

bool PepXMLParser3::hasPepProphet(){
  return bPepProphet;
}

bool PepXMLParser3::readFile(const char* fileName, double probFilter, double iProbFilter) {

  XML_ParserFree(parser);
	parser = XML_ParserCreate(NULL);
	XML_SetUserData(parser, this);
	XML_SetElementHandler(parser, PepXMLParser3_startElementCallback, PepXMLParser3_endElementCallback);
	XML_SetCharacterDataHandler(parser, PepXMLParser3_charactersCallback);

	vSpectra.clear();
  vScores.clear();
  vProteins.clear();
  vProtTable.clear();
  vMods.clear();
  vFiles.clear();
  vError.clear();
  vXL.clear();
  bIProphet=false;
  bPepProphet=false;
  analysisState = asNone;

  probabilityFilter=probFilter;
  iProbabilityFilter=iProbFilter;
	
	FILE* fptr=fopen(fileName,"rt");
	if (fptr == NULL){
		cerr << "Error parse(): No open file." << endl;
		return false;
	}

	char buffer[16384];
	int readBytes = 0;
	bool success = true;
	int chunk=0;

	while (success && (readBytes = (int) fread(buffer, 1, sizeof(buffer), fptr)) != 0){
		success = (XML_Parse(parser, buffer, readBytes, false) != 0);
    /*
    if(vSpectra.size()>500000) {
      fclose(fptr);
      XML_ParserFree(parser);
      cout << "CPepXMLAnalysis: " << sizeof(CPepXMLAnalysis) << endl;
      cout << "CPepXMLPeptide: " << sizeof(CPepXMLPeptide) << endl;
      cout << "CPepXMLPSM: " << sizeof(CPepXMLPSM) << endl;
      cout << "CPepXMLSearch: " << sizeof(CPepXMLSearch) << endl;
      cout << "CPepXMLSpectrum: " << sizeof(CPepXMLSpectrum) << endl;
      cout << "string: " << sizeof(string) << endl;
      cout << "vector: " << sizeof(vector<CPepXMLPSM>) << endl;
      return true;
    }
    */
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
		return false;
	}
	
	fclose(fptr);
	return true;
}

size_t PepXMLParser3::size(){
	return vSpectra.size();
}

bool PepXMLParser3::compareProt(const PepXMLProtTable& a, const PepXMLProtTable& b){
  return (a.name.compare(b.name)<0);
}

