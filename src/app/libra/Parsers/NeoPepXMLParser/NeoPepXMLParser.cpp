#include "NeoPepXMLParser.h"

using namespace std;

// Static callback handlers
static void CMzIdentML_startElementCallback(void *data, const XML_Char *el, const XML_Char **attr) {
  ((NeoPepXMLParser*)data)->startElement(el, attr);
}

static void CMzIdentML_endElementCallback(void *data, const XML_Char *el){
  ((NeoPepXMLParser*)data)->endElement(el);
}

static void CMzIdentML_charactersCallback(void *data, const XML_Char *s, int len){
  ((NeoPepXMLParser*)data)->characters(s, len);
}

NeoPepXMLParser::NeoPepXMLParser(){
  init();
}

NeoPepXMLParser::~NeoPepXMLParser(){
  XML_ParserFree(parser);
}

CnpxUIPSM& NeoPepXMLParser::operator[](const size_t& index){
  size_t pos = index;
  for (size_t i = 0; i<msms_pipeline_analysis.size(); i++){
    for (size_t j = 0; j<msms_pipeline_analysis[i].msms_run_summary.size(); j++){
      if(pos>=msms_pipeline_analysis[i].msms_run_summary[j].spectrum_query.size()){
        pos -= msms_pipeline_analysis[i].msms_run_summary[j].spectrum_query.size();
        continue;
      }
      psm.setPSM(msms_pipeline_analysis[i].msms_run_summary[j].spectrum_query[pos]);
      return psm;
    }
  }
  cerr << "ERROR, NeoPepXMLParser::operator[]: index out of bounds." << endl;
  exit(-50);
}

void NeoPepXMLParser::addMSMSPipelineAnalysis(std::string date, std::string summary_xml){
  CnpxMSMSPipelineAnalysis m;
  m.date.parseDateTime(date);
  m.summary_xml=summary_xml;
  msms_pipeline_analysis.push_back(m);
}

void NeoPepXMLParser::calcSize(){
  sz = 0;
  for (size_t i = 0; i<msms_pipeline_analysis.size(); i++){
    for (size_t j = 0; j<msms_pipeline_analysis[i].msms_run_summary.size(); j++){
      sz += msms_pipeline_analysis[i].msms_run_summary[j].spectrum_query.size();
    }
  }
}

void NeoPepXMLParser::characters(const XML_Char *s, int len) {
  /*
  switch (activeEl.back()){
  case PeptideSequence:
    char strbuf[1024];
    if (len>1024) {
      cout << "character buffer overrun" << endl;
      return;
    }
    strncpy(strbuf, s, len);
    strbuf[len] = '\0';
    sequenceCollection.peptide->back().peptideSequence.text += strbuf;
    break;
  default:
    //cout << "unprocessed characters: " << s << endl;
    break;
  }
  */
}

void NeoPepXMLParser::endElement(const XML_Char *el) {

  string s;
  for(int i=0;i<PEPXML_NUM_ELEMENTS;i++){
    if(isElement(elements[i].c_str(),el)){
      if(activeEl.back()!=(pepXMLElement)i) cout << "Error: unexpected end element: " << elements[i] << " should be " << elements[activeEl.back()] << endl;
      else activeEl.pop_back();
      break;
    }
  }

  if (isElement("interprophet_result", el)) {
    iProb= msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().interprophet_result.probability;
  } else if (isElement("msms_run_summary", el)) {
    if(rsFilter.size()>0){
      if(msms_pipeline_analysis.back().msms_run_summary.back().base_name.find(rsFilter)==string::npos) msms_pipeline_analysis.back().msms_run_summary.pop_back();
    }
  } else if (isElement("peptideprophet_result", el)) {
    pProb= msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().peptide_prophet_result.probability;
  } else if (isElement("search_hit", el)) {
    if (rsFilter.size() > 0) {
      if (msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().peptide.compare(shFilter) != string::npos) msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.pop_back();
    }
  } else if (isElement("spectrum_query", el)) {
    bool bPop=false;
    if(probFilter>-0.1){
      if(iProb>-0.1){
        if(iProb<probFilter && !bPop) {
          msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.pop_back();
          bPop=true;
        }
      } else if(pProb>-0.1){
        if (pProb < probFilter && !bPop) {
          msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.pop_back();
          bPop=true;
        }
      }
    }
    if (shFilter.size() > 0 && !bPop) {
      if (msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result[0].search_hit[0].peptide.compare(shFilter) != 0){
        msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.pop_back();
        bPop=true;
      }
    }
    iProb=-1;
    pProb=-1;
  }

}

void NeoPepXMLParser::init() {
  parser = XML_ParserCreate(NULL);
  XML_SetUserData(parser, this);
  XML_SetElementHandler(parser, CMzIdentML_startElementCallback, CMzIdentML_endElementCallback);
  XML_SetCharacterDataHandler(parser, CMzIdentML_charactersCallback);

  version = 22;
  probFilter=-1;
  rsFilter.clear();
  shFilter.clear();
  pProb=-1;
  iProb=-1;

  elements[pxAffectedChannel] = "affected_channel";
  elements[pxAlternativeProtein] = "alternative_protein";
  elements[pxAminoAcidModification] = "aminoacid_modification";
  elements[pxAminoAcidSubstitution] = "aminoacid_substitution";
  elements[pxAnalysisResult] = "analysis_result";
  elements[pxAnalysisSummary] = "analysis_summary";
  elements[pxAnalysisTimestamp] = "analysis_timestamp";
  elements[pxBin] = "bin";
  elements[pxContributingChannel] = "contributing_channel";
  elements[pxDatabaseRefreshTimestamp] = "database_refresh_timestamp";
  elements[pxDatasetDerivation] = "dataset_derivation";
  elements[pxDecoyAnalysis] = "decoy_analysis";
  elements[pxDecoyAnalysisSummary] = "decoy_analysis_summary";
  elements[pxDistributionPoint] = "distribution_point";
  elements[pxEnzymaticSearchConstraint] = "enzymatic_search_constraint";
  elements[pxErrorPoint] = "error_point";
  elements[pxFragmentMasses] = "fragment_masses";
  elements[pxInputfile] = "inputfile";
  elements[pxIntensity] = "intensity";
  elements[pxInteractSummary] = "interact_summary";
  elements[pxInterprophetResult] = "interprophet_result";
  elements[pxInterprophetSummary] = "interprophet_summary";
  elements[pxIsotopicContributions] = "isotopic_contributions";
  elements[pxLability] = "lability";
  elements[pxLibraResult] = "libra_result";
  elements[pxLibraSummary] = "libra_summary";
  elements[pxLinkedPeptide] = "linked_peptide";
  elements[pxMixture_Model] = "mixture_model";
  elements[pxMixturemodel] = "mixturemodel";
  elements[pxMixturemodelDistribution] = "mixturemodel_distribution";
  elements[pxModAminoAcidMass] = "mod_aminoacid_mass";
  elements[pxModAminoAcidProbability] = "mod_aminoacid_probability";
  elements[pxModificationInfo] = "modification_info";
  elements[pxModTerminalProbability] = "mod_terminal_probability";
  elements[pxMSMSPipelineAnalysis] = "msms_pipeline_analysis";
  elements[pxMSMSRunSummary] = "msms_run_summary";
  elements[pxNegmodelDistribution] = "negmodel_distribution";
  elements[pxParameter] = "parameter";
  elements[pxPeptideProphetResult] = "peptideprophet_result";
  elements[pxPeptideprophetSummary] = "peptideprophet_summary";
  elements[pxPepXMLQuantResult] = "pepxmlquant_result";
  elements[pxPoint] = "point";
  elements[pxPosmodelDistribution] = "posmodel_distribution";
  elements[pxPTMProphetResult] = "ptmprophet_result";
  elements[pxPTMProphetSummary] = "ptmprophet_summary";
  elements[pxQuanticResult] = "quantic_result";
  elements[pxQuanticSummary] = "quantic_summary";
  elements[pxROCDataPoint] = "roc_data_point";
  elements[pxROCErrorData] = "roc_error_data";
  elements[pxSampleEnzyme] = "sample_enzyme";
  elements[pxSearchDatabase] = "search_database";
  elements[pxSearchHit] = "search_hit";
  elements[pxSearchResult] = "search_result";
  elements[pxSearchScore] = "search_score";
  elements[pxSearchScoreSummary] = "search_score_summary";
  elements[pxSearchSummary] = "search_summary";
  elements[pxSpecificity] = "specificity";
  elements[pxSpectrumQuery] = "spectrum_query";
  elements[pxTerminalModification] = "terminal_modification";
  elements[pxXLink] = "xlink";
  elements[pxXLinkScore] = "xlink_score";
  elements[pxXpressLabelFreeResult] = "xpresslabelfree_result";
  elements[pxXpressLabelFreeSummary] = "xpresslabelfree_summary";
}

void NeoPepXMLParser::startElement(const XML_Char *el, const XML_Char **attr){

  //cout << el << endl; //for diagnostics

  //string s;
  if (isElement("affected_channel", el)){
    activeEl.push_back(pxAffectedChannel);
    CnpxAffectedChannel c;
    c.channel = atoi(getAttrValue("channel", attr));
    c.correction = atof(getAttrValue("correction", attr));
    msms_pipeline_analysis.back().analysis_summary.back().libra_summary.back().isotopic_contributions.back().contributing_channel.back().affected_channel.push_back(c);

  } else if (isElement("alternative_protein", el)){
    activeEl.push_back(pxAlternativeProtein);
    CnpxAlternativeProtein c;
    c.num_tol_term=atoi(getAttrValue("num_tol_term",attr));
    c.peptide_next_aa=getAttrValue("peptide_next_aa",attr);
    c.peptide_prev_aa = getAttrValue("peptide_prev_aa", attr);
    c.peptide_start_pos = atoi(getAttrValue("peptide_start_pos", attr));
    c.protein = getAttrValue("protein", attr);
    c.protein_descr = getAttrValue("protein_descr", attr);
    c.protein_link_pos_a = atoi(getAttrValue("protein_link_pos_a", attr));
    c.protein_link_pos_b = atoi(getAttrValue("protein_link_pos_b", attr));
    c.protein_mw = atof(getAttrValue("protein_mw", attr));
    switch (activeEl[activeEl.size() - 2]) {
    case pxSearchHit:
      msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().alternative_protein.push_back(c);
      break;
    case pxLinkedPeptide:
      msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().xlink.back().linked_peptide.back().alternative_protein.push_back(c);
      break;
    default:
      cout << "Error: stray alternative_protein element" << endl;
      break;
    }

  } else if (isElement("aminoacid_modification",el)){
    activeEl.push_back(pxAminoAcidModification);
    CnpxAminoAcidModification c;
    c.aminoacid=getAttrValue("aminoacid",attr);
    c.binary=getAttrValue("binary",attr);
    c.description=getAttrValue("description",attr);
    c.mass=(float)atof(getAttrValue("mass",attr));
    c.massdiff=(float)atof(getAttrValue("massdiff",attr));
    c.peptide_terminus=getAttrValue("peptide_terminus",attr);
    c.protein_terminus=getAttrValue("protein_terminus",attr);
    c.symbol=getAttrValue("symbol",attr);
    c.variable=getAttrValue("variable",attr);
    msms_pipeline_analysis.back().msms_run_summary.back().search_summary.back().aminoacid_modification.push_back(c);

  } else if (isElement("aminoacid_substitution", el)) {
    activeEl.push_back(pxAminoAcidSubstitution);
    CnpxAminoAcidSubstitution c;
    c.position = atoi(getAttrValue("position", attr));
    c.orig_aa = getAttrValue("orig_aa", attr);
    c.num_tol_term = atoi(getAttrValue("num_tol_term", attr));
    c.peptide_prev_aa = getAttrValue("peptide_prev_aa", attr);
    c.peptide_next_aa = getAttrValue("peptide_next_aa", attr);
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().modification_info.back().aminoacid_substitution.push_back(c);


  } else if (isElement("analysis_result", el)) {
    activeEl.push_back(pxAnalysisResult);
    CnpxAnalysisResult c;
    c.analysis = getAttrValue("analysis", attr);
    c.id = atoi(getAttrValue("id", attr));
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.push_back(c);

  } else if (isElement("analysis_summary", el)){
    activeEl.push_back(pxAnalysisSummary);
    CnpxAnalysisSummary c;
    c.time.parseDateTime(getAttrValue("time", attr));
    c.analysis=getAttrValue("analysis",attr);
    c.version=getAttrValue("version",attr);
    msms_pipeline_analysis.back().analysis_summary.push_back(c);

  } else if (isElement("analysis_timestamp", el)) {
    activeEl.push_back(pxAnalysisTimestamp);
    CnpxAnalysisTimestamp c;
    c.analysis = getAttrValue("analysis", attr);
    c.id = atoi(getAttrValue("id", attr));
    c.time.parseDateTime(getAttrValue("time", attr));
    msms_pipeline_analysis.back().msms_run_summary.back().analysis_timestamp.push_back(c);

  } else if (isElement("bin", el)) {
    activeEl.push_back(pxBin);
    npxBin c;
    c.pos_prob = atof(getAttrValue("pos_prob", attr));
    c.pos_prob = atof(getAttrValue("neg_prob", attr));
    if(strcmp(getAttrValue("value",attr),"true")==0) c.value=true;
    else c.value=false;
    msms_pipeline_analysis.back().analysis_summary.back().interprophet_summary.back().mixturemodel.back().bin.push_back(c);

  } else if (isElement("contributing_channel", el)){
    activeEl.push_back(pxContributingChannel);
    CnpxContributingChannel c;
    c.channel = atoi(getAttrValue("channel", attr));
    msms_pipeline_analysis.back().analysis_summary.back().libra_summary.back().isotopic_contributions.back().contributing_channel.push_back(c);

  } else if (isElement("database_refresh_timestamp",el)){
    activeEl.push_back(pxDatabaseRefreshTimestamp);
    CnpxDatabaseRefreshTimestamp c(true);
    c.database = getAttrValue("database", attr);
    c.min_num_enz_term = atoi(getAttrValue("min_num_enz_term", attr));
    msms_pipeline_analysis.back().msms_run_summary.back().analysis_timestamp.back().database_refresh_timestamp = c;

  } else if(isElement("dataset_derivation",el)){
    activeEl.push_back(pxDatasetDerivation);
    //TODO: Handle this

  } else if (isElement("decoy_analysis", el)){
    activeEl.push_back(pxDecoyAnalysis);
    CnpxDecoyAnalysis c;
    msms_pipeline_analysis.back().analysis_summary.back().decoy_analysis.push_back(c);

  } else if (isElement("decoy_analysis_summary", el)){
    activeEl.push_back(pxDecoyAnalysisSummary);
    CnpxDecoyAnalysisSummary c;
    c.decoy_ratio = atof(getAttrValue("decoy_ratio", attr));
    c.decoy_string=getAttrValue("decoy_string",attr);
    c.exclude_string = getAttrValue("exclude_string", attr);
    c.uniq_iproph_peps = getAttrValue("uniq_iproph_peps", attr);
    c.uniq_pproph_peps = getAttrValue("uniq_pproph_peps", attr);
    c.uniq_psm = getAttrValue("uniq_psm", attr);
    c.window_prob = getAttrValue("window_prob", attr);
    msms_pipeline_analysis.back().analysis_summary.back().decoy_analysis_summary.push_back(c);

  } else if (isElement("distribution_point", el)) {
    activeEl.push_back(pxDistributionPoint);
    CnpxDistributionPoint c;
    c.fvalue=atof(getAttrValue("fvalue",attr));
    c.obs_1_distr=atoi(getAttrValue("obs_1_distr",attr));
    c.model_1_pos_distr=atof(getAttrValue("model_1_pos_distr",attr));
    c.model_1_neg_distr=atof(getAttrValue("model_1_neg_distr",attr));
    c.obs_2_distr = atoi(getAttrValue("obs_2_distr", attr));
    c.model_2_pos_distr = atof(getAttrValue("model_2_pos_distr", attr));
    c.model_2_neg_distr = atof(getAttrValue("model_2_neg_distr", attr));
    c.obs_3_distr = atoi(getAttrValue("obs_3_distr", attr));
    c.model_3_pos_distr = atof(getAttrValue("model_3_pos_distr", attr));
    c.model_3_neg_distr = atof(getAttrValue("model_3_neg_distr", attr));
    c.obs_4_distr = atoi(getAttrValue("obs_4_distr", attr));
    c.model_4_pos_distr = atof(getAttrValue("model_4_pos_distr", attr));
    c.model_4_neg_distr = atof(getAttrValue("model_4_neg_distr", attr));
    c.obs_5_distr = atoi(getAttrValue("obs_5_distr", attr));
    c.model_5_pos_distr = atof(getAttrValue("model_5_pos_distr", attr));
    c.model_5_neg_distr = atof(getAttrValue("model_5_neg_distr", attr));
    c.obs_6_distr = atoi(getAttrValue("obs_6_distr", attr));
    c.model_6_pos_distr = atof(getAttrValue("model_6_pos_distr", attr));
    c.model_6_neg_distr = atof(getAttrValue("model_6_neg_distr", attr));
    c.obs_7_distr = atoi(getAttrValue("obs_7_distr", attr));
    c.model_7_pos_distr = atof(getAttrValue("model_7_pos_distr", attr));
    c.model_7_neg_distr = atof(getAttrValue("model_7_neg_distr", attr));
    switch (activeEl[activeEl.size() - 2]) {
    case pxPeptideprophetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().peptideprophet_summary.back().distribution_point.push_back(c);
      break;
    default:
      cout << "Error: stray distribution_point element" << endl;
      break;
    }

  } else if(isElement("enzymatic_search_constraint", el)){
    activeEl.push_back(pxEnzymaticSearchConstraint);
    CnpxEnzymaticSearchConstraint c;
    c.enzyme=getAttrValue("enzyme",attr);
    c.max_num_internal_cleavages=atoi(getAttrValue("max_num_internal_cleavages",attr));
    c.min_number_termini=atoi(getAttrValue("min_number_termini",attr));
    msms_pipeline_analysis.back().msms_run_summary.back().search_summary.back().enzymatic_search_constraint.push_back(c);

  } else if (isElement("error_point", el)){
    activeEl.push_back(pxErrorPoint);
    CnpxErrorPoint c;
    c.error=atof(getAttrValue("error",attr));
    c.min_prob= atof(getAttrValue("min_prob", attr));
    c.num_corr = atoi(getAttrValue("num_corr", attr));
    c.num_incorr = atoi(getAttrValue("num_incorr", attr));
    switch (activeEl[activeEl.size() - 3]) {
    case pxInterprophetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().interprophet_summary.back().roc_error_data.back().error_point.push_back(c);
      break;
    case pxPeptideprophetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().peptideprophet_summary.back().roc_error_data.back().error_point.push_back(c);
      break;
    case pxPTMProphetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().ptmprophet_summary.back().roc_error_data.back().error_point.push_back(c);
      break;
    default:
      cout << "Error: stray error_point element" << endl;
      break;
    }

  } else if (isElement("fragment_masses", el)){
    activeEl.push_back(pxFragmentMasses);
    CnpxFragmentMasses c;
    c.channel = atoi(getAttrValue("channel", attr));
    c.mz = atof(getAttrValue("mz", attr));
    c.offset = atof(getAttrValue("offset", attr));
    msms_pipeline_analysis.back().analysis_summary.back().libra_summary.back().fragment_masses.push_back(c);

  } else if (isElement("inputfile", el)){
    activeEl.push_back(pxInputfile);
    CnpxInputFile c;
    c.name=getAttrValue("name",attr);
    c.directory=getAttrValue("directory",attr);
    switch(activeEl[activeEl.size()-2]){
    case pxInteractSummary:
      msms_pipeline_analysis.back().analysis_summary.back().interact_summary.back().inputfile.push_back(c);
      break;
    case pxInterprophetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().interprophet_summary.back().inputfile.push_back(c);
      break;
    case pxPeptideprophetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().peptideprophet_summary.back().inputfile.push_back(c);
      break;
    case pxPTMProphetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().ptmprophet_summary.back().inputfile.push_back(c);
      break;
    case pxQuanticSummary:
      msms_pipeline_analysis.back().analysis_summary.back().quantic_summary.back().inputfile.push_back(c);
      break;
    default:
      cout << "Error: stray inputfile element" << endl;
      break;
    }

  } else if (isElement("intensity", el)) {
    activeEl.push_back(pxIntensity);
    CnpxIntensity c;
    c.channel = atoi(getAttrValue("channel", attr));
    c.absolute = atof(getAttrValue("absolute", attr));
    c.target_mass = atof(getAttrValue("target_mass", attr));
    c.normalized = atof(getAttrValue("normalized", attr));
    string s = getAttrValue("reject", attr);
    if(s.size()>0 && s[0]!='0') c.reject=true;
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().libra_result.intensity.push_back(c);

  } else if (isElement("interact_summary", el)){
    activeEl.push_back(pxInteractSummary);
    CnpxInteractSummary c;
    c.directory = getAttrValue("directory", attr);
    c.filename = getAttrValue("filename", attr);
    msms_pipeline_analysis.back().analysis_summary.back().interact_summary.push_back(c);

  } else if (isElement("interprophet_result", el)) {
    activeEl.push_back(pxInterprophetResult);
    CnpxInterprophetResult c(true);
    c.all_ntt_prob = getAttrValue("all_ntt_prob", attr);
    c.probability = atof(getAttrValue("probability", attr));
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().interprophet_result = c;

  } else if (isElement("interprophet_summary", el)){
    activeEl.push_back(pxInterprophetSummary);
    CnpxInterprophetSummary c;
    c.est_tot_num_correct_pep=atof(getAttrValue("est_tot_num_correct_pep",attr));
    c.est_tot_num_correct_psm = atof(getAttrValue("est_tot_num_correct_psm", attr));
    c.options=getAttrValue("options",attr);
    c.version=getAttrValue("version",attr);
    msms_pipeline_analysis.back().analysis_summary.back().interprophet_summary.push_back(c);

  } else if (isElement("isotopic_contributions", el)) {
    activeEl.push_back(pxIsotopicContributions);
    CnpxIsotopicContributions c;
    msms_pipeline_analysis.back().analysis_summary.back().libra_summary.back().isotopic_contributions.push_back(c);

  } else if (isElement("lability", el)) {
    activeEl.push_back(pxLability);
    CnpxLability c;
    c.numlosses = atoi(getAttrValue("numlosses", attr));
    c.pval = atof(getAttrValue("pval", attr));
    c.probability = atof(getAttrValue("probability", attr));
    c.oscore = atof(getAttrValue("oscore", attr));
    c.mscore = atof(getAttrValue("mscore", attr));
    c.cterm_score = atof(getAttrValue("cterm_score", attr));
    c.nterm_score = atof(getAttrValue("nterm_score", attr));
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().ptmprophet_result.back().lability.push_back(c);


  } else if (isElement("libra_result", el)) {
    activeEl.push_back(pxLibraResult);
    CnpxLibraResult c(true);
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().libra_result = c;

  } else if (isElement("libra_summary", el)){
    activeEl.push_back(pxLibraSummary);
    CnpxLibraSummary c;
    c.mass_tolerance = atof(getAttrValue("mass_tolerance", attr));
    c.centroiding_preference = atoi(getAttrValue("centroiding_preference", attr));
    c.normalization = atoi(getAttrValue("normalization", attr));
    c.output_type = atoi(getAttrValue("output_type", attr));
    c.channel_code = getAttrValue("channel_code", attr);
    msms_pipeline_analysis.back().analysis_summary.back().libra_summary.push_back(c);

  } else if (isElement("linked_peptide", el)){
    activeEl.push_back(pxLinkedPeptide);
    CnpxLinkedPeptide c;
    c.calc_neutral_pep_mass = atof(getAttrValue("calc_neutral_pep_mass", attr));
    c.num_tot_proteins = atoi(getAttrValue("num_tot_proteins", attr));
    c.peptide = getAttrValue("peptide", attr);
    c.peptide_next_aa = getAttrValue("peptide_next_aa", attr);
    c.peptide_prev_aa = getAttrValue("peptide_prev_aa", attr);
    c.peptide_start_pos = atoi(getAttrValue("peptide_start_pos", attr));
    c.protein = getAttrValue("protein", attr);
    c.protein_link_pos_a = atoi(getAttrValue("protein_link_pos_a", attr));
    c.designation = getAttrValue("designation", attr);
    c.complement_mass = atof(getAttrValue("complement_mass", attr));
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().xlink.back().linked_peptide.push_back(c);

  } else if(isElement("mixture_model",el)){
    activeEl.push_back(pxMixture_Model);
    CnpxMixture_Model c;
    c.comments=getAttrValue("comments",attr);
    c.est_tot_correct=atof(getAttrValue("est_tot_correct",attr));
    c.num_iterations=atoi(getAttrValue("num_iterations",attr));
    c.precursor_ion_charge = atoi(getAttrValue("precursor_ion_charge", attr));
    c.prior_probability = atof(getAttrValue("prior_proability", attr));
    c.tot_num_spectra = atoi(getAttrValue("tot_num_spectra", attr));
    msms_pipeline_analysis.back().analysis_summary.back().peptideprophet_summary.back().mixture_model.push_back(c);

  } else if (isElement("mixturemodel", el)){
    activeEl.push_back(pxMixturemodel);
    CnpxMixtureModel c;
    c.name = getAttrValue("name",attr);
    c.neg_bandwidth = (float)atof(getAttrValue("neg_bandwidth", attr));
    c.pos_bandwidth = (float)atof(getAttrValue("pos_bandwidth",attr));
    switch (activeEl[activeEl.size() - 2]) {
    case pxInterprophetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().interprophet_summary.back().mixturemodel.push_back(c);
      break;
    case pxMixture_Model:
      msms_pipeline_analysis.back().analysis_summary.back().peptideprophet_summary.back().mixture_model.back().mixturemodel.push_back(c);
      break;
    case pxPTMProphetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().ptmprophet_summary.back().mixturemodel.push_back(c);
      break;
    default:
      cout << "Error: stray mixturemodel element" << endl;
      break;
    }

  } else if (isElement("mixturemodel_distribution", el)) {
    activeEl.push_back(pxMixturemodelDistribution);
    CnpxMixtureModelDistribution c;
    c.name=getAttrValue("name",attr);
    if (activeEl[activeEl.size() - 2]==pxInterprophetSummary) {
      msms_pipeline_analysis.back().analysis_summary.back().interprophet_summary.back().mixturemodel_distribution.push_back(c);
    } else if (activeEl[activeEl.size() - 3] == pxPeptideprophetSummary) {
      msms_pipeline_analysis.back().analysis_summary.back().peptideprophet_summary.back().mixture_model.back().mixturemodel_distribution.push_back(c);
    } else {
      cout << "Error: stray mixturemodel_distribution element" << endl;
    }

  } else if (isElement("mod_aminoacid_mass", el)) {
    activeEl.push_back(pxModAminoAcidMass);
    CnpxModAminoAcidMass c;
    c.id = getAttrValue("id", attr);
    c.mass = atof(getAttrValue("mass", attr));
    c.position = atoi(getAttrValue("position", attr));
    c.source = getAttrValue("source", attr);
    c.staticMass = atof(getAttrValue("static", attr));
    c.variable = atof(getAttrValue("variable", attr));
    switch (activeEl[activeEl.size() - 3]) {
    case pxSearchHit:
      msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().modification_info.back().mod_aminoacid_mass.push_back(c);
      break;
    case pxLinkedPeptide:
      msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().xlink.back().linked_peptide.back().modification_info.back().mod_aminoacid_mass.push_back(c);
      break;
    default:
      cout << "Error: stray mod_aminoacid_mass element" << endl;
      break;
    }
   
  } else if (isElement("mod_aminoacid_probability", el)) {
    activeEl.push_back(pxModAminoAcidProbability);
    CnpxModAminoAcidProbability c;
    c.position = atoi(getAttrValue("position", attr));
    c.probability = atof(getAttrValue("probability", attr));
    c.oscore = atof(getAttrValue("oscore", attr));
    c.mscore = atof(getAttrValue("mscore", attr));
    c.direct_oscore = atof(getAttrValue("direct_oscore", attr));
    c.direct_mscore = atof(getAttrValue("direct_mscore", attr));
    c.cterm_score = atof(getAttrValue("cterm_score", attr));
    c.nterm_score = atof(getAttrValue("nterm_score", attr));
    c.shift = getAttrValue("shift",attr)[0];
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().ptmprophet_result.back().mod_amino_acid_probability.push_back(c);

  } else if (isElement("modification_info", el)) {
    activeEl.push_back(pxModificationInfo);
    CnpxModificationInfo c;
    c.modified_peptide = getAttrValue("modified_peptide", attr);
    c.mod_cterm_mass = atof(getAttrValue("mod_cterm_mass", attr));
    c.mod_nterm_mass = atof(getAttrValue("mod_nterm_mass", attr));
    switch (activeEl[activeEl.size() - 2]) {
    case pxSearchHit:
      msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().modification_info.push_back(c);
      break;
    case pxLinkedPeptide:
      msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().xlink.back().linked_peptide.back().modification_info.push_back(c);
      break;
    default:
      cout << "Error: stray modification_info element" << endl;
      break;
    }

  } else if (isElement("mod_terminal_probability", el)) {
    activeEl.push_back(pxModTerminalProbability);
    CnpxModTerminalProbability c;
    c.terminus = getAttrValue("position", attr)[0];
    c.probability = atof(getAttrValue("probability", attr));
    c.oscore = atof(getAttrValue("oscore", attr));
    c.mscore = atof(getAttrValue("mscore", attr));
    c.direct_oscore = atof(getAttrValue("direct_oscore", attr));
    c.direct_mscore = atof(getAttrValue("direct_mscore", attr));
    c.cterm_score = atof(getAttrValue("cterm_score", attr));
    c.nterm_score = atof(getAttrValue("nterm_score", attr));
    c.shift = getAttrValue("shift", attr)[0];
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().ptmprophet_result.back().mod_terminal_probability.push_back(c);

  } else if (isElement("msms_pipeline_analysis", el)){
    activeEl.push_back(pxMSMSPipelineAnalysis);
    CnpxMSMSPipelineAnalysis c;
    c.date.parseDateTime(getAttrValue("date",attr));
    c.summary_xml=getAttrValue("summary_xml",attr);
    msms_pipeline_analysis.push_back(c);

  } else if (isElement("msms_run_summary", el)){
    activeEl.push_back(pxMSMSRunSummary);
    CnpxMSMSRunSummary c;
    c.base_name=getAttrValue("base_name",attr);
    c.msDetector=getAttrValue("msDetector",attr);
    c.msIonization=getAttrValue("msIonization",attr);
    c.msManufacturer=getAttrValue("msManufacturer",attr);
    c.msMassAnalyzer=getAttrValue("msMassAnalyzer",attr);
    c.msModel=getAttrValue("msModel",attr);
    c.raw_data=getAttrValue("raw_data",attr);
    c.raw_data_type=getAttrValue("raw_data_type",attr);
    msms_pipeline_analysis.back().msms_run_summary.push_back(c);

  } else if (isElement("negmodel_distribution", el)) {
    activeEl.push_back(pxNegmodelDistribution);
    CnpxNegModelDistribution c;
    c.type = getAttrValue("type", attr);
    if (activeEl[activeEl.size() - 3]==pxInterprophetSummary) {
      msms_pipeline_analysis.back().analysis_summary.back().interprophet_summary.back().mixturemodel_distribution.back().negmodel_distribution.push_back(c);
    } else if (activeEl[activeEl.size() - 4] == pxPeptideprophetSummary){
      msms_pipeline_analysis.back().analysis_summary.back().peptideprophet_summary.back().mixture_model.back().mixturemodel_distribution.back().negmodel_distribution.push_back(c);
    } else {
      cout << "Error: stray negmodel_distribution element" << endl;
    }

  } else if (isElement("parameter", el)) {
    activeEl.push_back(pxParameter);
    CnpxParameter c;
    c.name = getAttrValue("name", attr);
    c.type = getAttrValue("type", attr);
    c.value = getAttrValue("value", attr);
    switch (activeEl[activeEl.size() - 2]) {
    case pxAnalysisResult:
      msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().parameter.push_back(c);
      break;
    case pxSearchScoreSummary:
      if (activeEl[activeEl.size() - 3] == pxPeptideProphetResult) {
        msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().peptide_prophet_result.search_score_summary.parameter.push_back(c);
      } else if (activeEl[activeEl.size() - 3] == pxInterprophetResult) {
        msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().interprophet_result.search_score_summary.parameter.push_back(c);
      } else if (activeEl[activeEl.size() - 3] == pxPepXMLQuantResult) {
        msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().pepxmlquant_result.search_score_summary.parameter.push_back(c);
      } else {
        cout << "Unknown location for parameter: " << pxSearchScoreSummary << endl;
      }
      break;
    case pxSearchSummary:
      msms_pipeline_analysis.back().msms_run_summary.back().search_summary.back().parameter.push_back(c);
      break;
    case pxNegmodelDistribution:
      if (activeEl[activeEl.size() - 5] == pxPeptideprophetSummary) {
        msms_pipeline_analysis.back().analysis_summary.back().peptideprophet_summary.back().mixture_model.back().mixturemodel_distribution.back().negmodel_distribution.back().parameter.push_back(c);
      } else if (activeEl[activeEl.size() - 4] == pxInterprophetSummary) {
        msms_pipeline_analysis.back().analysis_summary.back().interprophet_summary.back().mixturemodel_distribution.back().negmodel_distribution.back().parameter.push_back(c);
      } else {
        cout << "Unknown location for parameter: " << pxNegmodelDistribution << endl;
      }
      break;
    case pxPosmodelDistribution:
      if (activeEl[activeEl.size() - 5] == pxPeptideprophetSummary) {
        msms_pipeline_analysis.back().analysis_summary.back().peptideprophet_summary.back().mixture_model.back().mixturemodel_distribution.back().posmodel_distribution.back().parameter.push_back(c);
      } else if (activeEl[activeEl.size() - 4] == pxInterprophetSummary) {
        msms_pipeline_analysis.back().analysis_summary.back().interprophet_summary.back().mixturemodel_distribution.back().posmodel_distribution.back().parameter.push_back(c);
      } else {
        cout << "Unknown location for parameter: " << pxPosmodelDistribution << endl;
      }
      break;
    case pxPTMProphetResult:
      msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().ptmprophet_result.back().parameter.push_back(c);
      break;
    default:
      cout << "Error: stray parameter element: " << elements[activeEl[activeEl.size() - 2]] << endl;
      break;
    }

  } else if (isElement("peptideprophet_result",el)){
    activeEl.push_back(pxPeptideProphetResult);
    CnpxPeptideProphetResult c(true);
    c.all_ntt_prob = getAttrValue("all_ntt_prob", attr);
    c.analysis = getAttrValue("analysis", attr);
    c.probability = atof(getAttrValue("probability", attr));
    c.pep1_probability = atof(getAttrValue("pep1_probability",attr));
    c.pep2_probability = atof(getAttrValue("pep2_probability", attr));
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().peptide_prophet_result = c;

  } else if (isElement("peptideprophet_summary", el)) {
    activeEl.push_back(pxPeptideprophetSummary);
    CnpxPeptideprophetSummary c;
    c.author = getAttrValue("author",attr);
    c.est_tot_num_correct=atof(getAttrValue("est_tot_num_correct",attr));
    c.min_prob= atof(getAttrValue("min_prob", attr));
    c.options = getAttrValue("options", attr);
    c.type=getAttrValue("type",attr);
    c.version = getAttrValue("version", attr);
    msms_pipeline_analysis.back().analysis_summary.back().peptideprophet_summary.push_back(c);

  } else if (isElement("pepxmlquant_result", el)){
    activeEl.push_back(pxPepXMLQuantResult);
    CnpxPepXMLQuantResult c(true);
    c.area = atof(getAttrValue("area", attr));
    c.retention_time_sec = atof(getAttrValue("retention_time_sec", attr));
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().pepxmlquant_result = c;

  } else if (isElement("point", el)){
    activeEl.push_back(pxPoint);
    npxPointD c;
    npxPointM m;  
    switch (activeEl[activeEl.size() - 2]){
    case pxDecoyAnalysis:
      c.fdr_ip = atof(getAttrValue("fdr_ip", attr));
      c.fdr_ip_decoy = atof(getAttrValue("fdr_ip_decoy", attr));
      c.fdr_pp = atof(getAttrValue("fdr_pp", attr));
      c.fdr_pp_decoy = atof(getAttrValue("fdr_pp_decoy", attr));
      c.ip_decoy_uncert = atof(getAttrValue("ip_decoy_uncert", attr));
      c.ip_uncert = atof(getAttrValue("ip_uncert", attr));
      c.num_corr_ip = atof(getAttrValue("num_corr_ip", attr));
      c.num_corr_ip_decoy = atof(getAttrValue("num_corr_ip_decoy", attr));
      c.num_corr_pp = atof(getAttrValue("num_corr_pp", attr));
      c.num_corr_pp_decoy = atof(getAttrValue("num_corr_pp_decoy", attr));
      c.pp_decoy_uncert = atof(getAttrValue("pp_decoy_uncert", attr));
      c.pp_uncert = atof(getAttrValue("pp_uncert", attr));
      c.prob_cutoff = atof(getAttrValue("prob_cutoff", attr));
      msms_pipeline_analysis.back().analysis_summary.back().decoy_analysis.back().point.push_back(c);
      break;
    case pxMixturemodel:
      m.neg_dens = (float)atof(getAttrValue("neg_dens", attr));
      m.neg_obs_dens = (float)atof(getAttrValue("neg_obs_dens", attr));
      m.pos_dens = (float)atof(getAttrValue("pos_dens", attr));
      m.pos_obs_dens = (float)atof(getAttrValue("pos_obs_dens", attr));
      m.value = (float)atof(getAttrValue("value", attr));
      switch (activeEl[activeEl.size()-3]){
      case pxInterprophetSummary:
        msms_pipeline_analysis.back().analysis_summary.back().interprophet_summary.back().mixturemodel.back().point.push_back(m);
        break;
      case pxMixture_Model:
        msms_pipeline_analysis.back().analysis_summary.back().peptideprophet_summary.back().mixture_model.back().mixturemodel.back().point.push_back(m);
        break;
      case pxPTMProphetSummary:
        msms_pipeline_analysis.back().analysis_summary.back().ptmprophet_summary.back().mixturemodel.back().point.push_back(m);
        break;
      default:
        cout << "Error: stray point element" << endl;
        break;
      }
      break;
    default:
      cout << "Error: stray point element" << endl;
      break;
    }

  } else if (isElement("posmodel_distribution", el)) {
    activeEl.push_back(pxPosmodelDistribution);
    CnpxPosModelDistribution c;
    c.type = getAttrValue("type",attr);
    if (activeEl[activeEl.size() - 3] == pxInterprophetSummary) {
      msms_pipeline_analysis.back().analysis_summary.back().interprophet_summary.back().mixturemodel_distribution.back().posmodel_distribution.push_back(c);
    } else if (activeEl[activeEl.size() - 4] == pxPeptideprophetSummary) {
      msms_pipeline_analysis.back().analysis_summary.back().peptideprophet_summary.back().mixture_model.back().mixturemodel_distribution.back().posmodel_distribution.push_back(c);
    } else {
      cout << "Error: stray posmodel_distribution element" << endl;
    }

  } else if (isElement("ptmprophet_result", el)) {
    activeEl.push_back(pxPTMProphetResult);
    CnpxPTMProphetResult c(true);
    c.ptm = getAttrValue("ptm", attr);
    c.prior = getAttrValue("prior", attr);
    c.ptm_peptide = getAttrValue("ptm_peptide", attr);
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().ptmprophet_result.push_back(c);

  } else if (isElement("ptmprophet_summary", el)) {
    activeEl.push_back(pxPTMProphetSummary);
    CnpxPTMProphetSummary c;
    c.frag_ppm_tol = getAttrValue("frag_ppm_tol", attr);
    c.min_o = getAttrValue("min_o", attr);
    c.min_o_factors = getAttrValue("min_o_factors", attr);
    c.mod_string = getAttrValue("mod_string", attr);
    c.options = getAttrValue("options", attr);
    c.version = getAttrValue("version", attr);
    msms_pipeline_analysis.back().analysis_summary.back().ptmprophet_summary.push_back(c);

  } else if (isElement("quantic_result", el)) {
    activeEl.push_back(pxQuanticResult);
    CnpxQuanticResult c(true);
    c.antic = atof(getAttrValue("antic", attr));
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().quantic_result = c;

  } else if (isElement("quantic_summary", el)) {
    activeEl.push_back(pxQuanticSummary);
    CnpxQuanticSummary c;
    c.options = getAttrValue("options", attr);
    c.version = getAttrValue("version", attr);
    msms_pipeline_analysis.back().analysis_summary.back().quantic_summary.push_back(c);


  } else if (isElement("roc_error_data", el)){
    activeEl.push_back(pxROCErrorData);
    CnpxROCErrorData c;
    c.charge=getAttrValue("charge",attr);
    c.charge_est_correct=atof(getAttrValue("charge_est_correct",attr));
    switch (activeEl[activeEl.size() - 2]) {
    case pxInterprophetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().interprophet_summary.back().roc_error_data.push_back(c);
      break;
    case pxPeptideprophetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().peptideprophet_summary.back().roc_error_data.push_back(c);
      break;
    case pxPTMProphetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().ptmprophet_summary.back().roc_error_data.push_back(c);
      break;
    default:
      cout << "Error: stray roc_error_data element" << endl;
      break;
    }

  } else if (isElement("roc_data_point", el)){
    activeEl.push_back(pxROCDataPoint);
    CnpxROCDataPoint c;
    c.error = atof(getAttrValue("error", attr));
    c.min_prob = atof(getAttrValue("min_prob", attr));
    c.num_corr = atoi(getAttrValue("num_corr", attr));
    c.num_incorr = atoi(getAttrValue("num_incorr", attr));
    c.sensitivity = atof(getAttrValue("sensitivity", attr));
    switch (activeEl[activeEl.size() - 3]) {
    case pxInterprophetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().interprophet_summary.back().roc_error_data.back().roc_data_point.push_back(c);
      break;
    case pxPeptideprophetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().peptideprophet_summary.back().roc_error_data.back().roc_data_point.push_back(c);
      break;
    case pxPTMProphetSummary:
      msms_pipeline_analysis.back().analysis_summary.back().ptmprophet_summary.back().roc_error_data.back().roc_data_point.push_back(c);
      break;
    default:
      cout << "Error: stray roc_data_point element" << endl;
      break;
    }

  } else if (isElement("sample_enzyme", el)){
    activeEl.push_back(pxSampleEnzyme);
    CnpxSampleEnzyme c;
    c.description=getAttrValue("description",attr);
    c.fidelity=getAttrValue("fidelity",attr);
    if (strlen(getAttrValue("independent", attr))>0){
      if(getAttrValue("independent",attr)[0]=='0') c.independent=false;
      else c.independent=true;
    }
    c.name=getAttrValue("name",attr);
    msms_pipeline_analysis.back().msms_run_summary.back().sample_enzyme.push_back(c);

  } else if (isElement("search_database",el)){
    activeEl.push_back(pxSearchDatabase);
    CnpxSearchDatabase c;
    c.database_name=getAttrValue("database_name",attr);
    c.database_release_date.parseDateTime(getAttrValue("database_release_date",attr));
    c.database_release_identifier=getAttrValue("database_release_identifier",attr);
    c.local_path=getAttrValue("local_path",attr);
    c.orig_database_url=getAttrValue("orig_database_url",attr);
    c.size_in_db_entries=atoi(getAttrValue("size_in_db_entries",attr));
    c.size_of_residues=atoi(getAttrValue("size_of_residues",attr));
    c.type=getAttrValue("type",attr);
    c.URL=getAttrValue("URL",attr);
    msms_pipeline_analysis.back().msms_run_summary.back().search_summary.back().search_database.push_back(c);

  } else if(isElement("search_hit",el)){
    activeEl.push_back(pxSearchHit);
    CnpxSearchHit c;
    c.calc_neutral_pep_mass=atof(getAttrValue("calc_neutral_pep_mass",attr));
    c.calc_pI=atof(getAttrValue("calc_pI",attr));
    c.hit_rank=atoi(getAttrValue("hit_rank",attr));
    c.is_rejected=atoi(getAttrValue("is_rejected",attr));
    c.massdiff=atof(getAttrValue("massdiff",attr));
    c.num_matched_ions=atoi(getAttrValue("num_matched_ions",attr));
    c.num_matched_peptides=atoi(getAttrValue("num_matched_peptides",attr));
    c.num_missed_cleavages=atoi(getAttrValue("num_missed_cleavages",attr));
    c.num_tol_term=atoi(getAttrValue("num_tol_term",attr));
    c.num_tot_proteins=atoi(getAttrValue("num_tot_proteins",attr));
    c.peptide=getAttrValue("peptide",attr);
    c.peptide_next_aa=getAttrValue("peptide_next_aa",attr);
    c.peptide_prev_aa=getAttrValue("peptide_prev_aa",attr);
    c.peptide_start_pos = atoi(getAttrValue("peptide_start_pos", attr));
    c.protein=getAttrValue("protein",attr);
    c.protein_descr=getAttrValue("protein_descr",attr);
    c.protein_link_pos_a = atoi(getAttrValue("protein_link_pos_a", attr));
    c.protein_link_pos_b = atoi(getAttrValue("protein_link_pos_b", attr));
    c.protein_mw=atof(getAttrValue("protein_mw",attr));
    c.tot_num_ions=atoi(getAttrValue("tot_num_ions",attr));
    c.xlink_type=getAttrValue("xlink_type",attr);
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.push_back(c);

  } else if(isElement("search_result",el)){
    activeEl.push_back(pxSearchResult);
    CnpxSearchResult c;
    c.search_id=atoi(getAttrValue("search_id",attr));
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.push_back(c);

  } else if (isElement("search_score", el)) {
    activeEl.push_back(pxSearchScore);
    CnpxSearchScore c;
    c.name = getAttrValue("name", attr);
    c.value = getAttrValue("value", attr);
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().search_score.push_back(c);

  } else if (isElement("search_score_summary",el)){
    activeEl.push_back(pxSearchScoreSummary);
    CnpxSearchScoreSummary c(true);
    switch (activeEl[activeEl.size() - 2]) {
    case pxInterprophetResult:
      msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().interprophet_result.search_score_summary = c;
      break;
    case pxPeptideProphetResult:
      msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().peptide_prophet_result.search_score_summary = c;
      break;
    case pxPepXMLQuantResult:
      msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().pepxmlquant_result.search_score_summary = c;
      break;
    default:
      cout << "Error: stray search_score_summary element" << endl;
      break;
    }

  } else if (isElement("search_summary",el)){
    activeEl.push_back(pxSearchSummary);
    CnpxSearchSummary c;
    c.base_name=getAttrValue("base_name",attr);
    c.fragment_mass_type=getAttrValue("fragment_mass_type",attr);
    c.precursor_mass_type=getAttrValue("precursor_mass_type",attr);
    c.search_engine=getAttrValue("search_engine",attr);
    c.search_engine_version=getAttrValue("search_engine_version",attr);
    c.search_id=atoi(getAttrValue("search_id",attr));
    msms_pipeline_analysis.back().msms_run_summary.back().search_summary.push_back(c);

  } else if (isElement("specificity", el)){
    activeEl.push_back(pxSpecificity);
    CnpxSpecificity c;
    c.cut=getAttrValue("cut",attr);
    c.min_spacing=atoi(getAttrValue("min_spacing",attr));
    c.no_cut=getAttrValue("no_cut",attr);
    c.sense=getAttrValue("sense",attr);
    msms_pipeline_analysis.back().msms_run_summary.back().sample_enzyme.back().specificity.push_back(c);

  } else if(isElement("spectrum_query",el)){
    activeEl.push_back(pxSpectrumQuery);
    CnpxSpectrumQuery c;
    c.activation_method=getAttrValue("activation_method",attr);
    c.assumed_charge=atoi(getAttrValue("assumed_charge",attr));
    c.collision_energy=atof(getAttrValue("collison_energy",attr));
    c.compensation_voltage=atof(getAttrValue("compensation_voltage",attr));
    c.end_scan=atoi(getAttrValue("end_scan",attr));
    c.index=atoi(getAttrValue("index",attr));
    c.precursor_intensity=atof(getAttrValue("precursor_intensity",attr));
    c.precursor_neutral_mass=atof(getAttrValue("precursor_neutral_mass",attr));
    c.retention_time_sec=atof(getAttrValue("retention_time_sec",attr));
    c.search_specification=getAttrValue("search_specification",attr);
    c.spectrum=getAttrValue("spectrum",attr);
    c.spectrumNativeID=getAttrValue("spectrumNativeID",attr);
    c.start_scan=atoi(getAttrValue("start_scan",attr));
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.push_back(c);

  } else if (isElement("terminal_modification", el)){
    activeEl.push_back(pxTerminalModification);
    CnpxTerminalModification c;
    c.description = getAttrValue("description", attr);
    c.mass = (float)atof(getAttrValue("mass", attr));
    c.massdiff = (float)atof(getAttrValue("massdiff", attr));
    c.terminus = getAttrValue("terminus", attr);
    c.protein_terminus = getAttrValue("protein_terminus", attr);
    c.symbol = getAttrValue("symbol", attr);
    c.variable = getAttrValue("variable", attr);
    msms_pipeline_analysis.back().msms_run_summary.back().search_summary.back().terminal_modification.push_back(c);

  } else if (isElement("xlink", el)){
    activeEl.push_back(pxXLink);
    CnpxXLink c;
    c.identifier = getAttrValue("identifier",attr);
    c.mass = atof(getAttrValue("mass", attr));
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().xlink.push_back(c);

  } else if (isElement("xlink_score", el)) {
    activeEl.push_back(pxXLinkScore);
    CnpxXLinkScore c;
    c.name = getAttrValue("name", attr);
    c.value = getAttrValue("value", attr);
    c.type = getAttrValue("type",attr);
    switch (activeEl[activeEl.size() - 2]) {
    case pxLinkedPeptide:
      msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().xlink.back().linked_peptide.back().xlink_score.push_back(c);
      break;
    case pxXLink:
      msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().xlink.back().xlink_score.push_back(c);
      break;
    default:
      cout << "Error: stray xlink_score element: " << (int)activeEl[activeEl.size() - 2] << endl;
      break;
    }

  } else if (isElement("xpresslabelfree_result", el)) {
    activeEl.push_back(pxXpressLabelFreeResult);
    CnpxXpressLabelFreeResult c(true);
    c.charge = atoi(getAttrValue("charge", attr));
    c.first_scan = atoi(getAttrValue("first_scan", attr));
    c.last_scan = atoi(getAttrValue("last_scan", attr));
    c.first_scan_RT_seconds = getAttrValue("first_scan_RT_seconds", attr);
    c.last_scan_RT_seconds = getAttrValue("last_scan_RT_seconds", attr);
    c.precursor_mz = getAttrValue("precursor_mz", attr);
    c.peak_area = getAttrValue("peak_area", attr);
    c.peak_intensity = getAttrValue("peak_intensity", attr);
    c.peak_intensity_RT_seconds = getAttrValue("peak_intensity_RT_seconds", attr);
    c.peak_intensity_scan = atoi(getAttrValue("peak_intensity_scan", attr));
    msms_pipeline_analysis.back().msms_run_summary.back().spectrum_query.back().search_result.back().search_hit.back().analysis_result.back().expresslabelfree_result = c;

  } else if (isElement("xpresslabelfree_summary", el)) {
    activeEl.push_back(pxXpressLabelFreeSummary);
    CnpxXpressLabelFreeSummary c;
    c.masstol = getAttrValue("masstol", attr);
    c.ppmtol = getAttrValue("ppmtol", attr);
    c.min_num_chromatogram_points = getAttrValue("min_num_chromatogram_points", attr);
    c.min_num_isotope_peaks = getAttrValue("min_num_isotope_peaks", attr);
    c.author = getAttrValue("author", attr);
    c.version = getAttrValue("version", attr);
    msms_pipeline_analysis.back().analysis_summary.back().xpresslabelfree_summary.push_back(c);

  } else if (isElement("xpresslabelfree_timestamp", el)) {
    cout << "Ignoring: xpresslabelfree_timestamp" << endl;
    

  } else {
    cout << "WARNING: Element undefined: " << el << endl;
    exit(1);
  }
}

bool NeoPepXMLParser::read(const char* fn){
  XML_ParserFree(parser);
  parser = XML_ParserCreate(NULL);
  XML_SetUserData(parser, this);
  XML_SetElementHandler(parser, CMzIdentML_startElementCallback, CMzIdentML_endElementCallback);
  XML_SetCharacterDataHandler(parser, CMzIdentML_charactersCallback);

  // clear data
  msms_pipeline_analysis.clear();
  FILE* fptr = fopen(fn, "rb");
  if (fptr == NULL){
    cerr << "Error parse(): No open file." << endl;
    return false;
  }
  npxfseek(fptr, 0, SEEK_END);
  f_off iEOF=npxftell(fptr);
  fclose(fptr);

  fptr = fopen(fn, "rt");

  int iTmp;
  f_off prog=0;
  int iPercent = 0;
  printf("%2d%%", iPercent);
  fflush(stdout);

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
    prog += sizeof(buffer);
    iTmp = (int)((double)prog / iEOF * 100);
    if (iTmp>iPercent){
      iPercent = iTmp;
      printf("\b\b\b%2d%%", iPercent);
      fflush(stdout);
    }
  }
  success = success && (XML_Parse(parser, buffer, 0, true) != 0);
  cout << endl;

  if (!success) {
    XML_Error error = XML_GetErrorCode(parser);

    cerr << fn << "(" << XML_GetCurrentLineNumber(parser) << ") : error " << (int)error << ": ";
    switch (error) {
    case XML_ERROR_SYNTAX:
      cerr << "Syntax error parsing XML.";
      break;
    case XML_ERROR_INVALID_TOKEN:
      cerr << "XML invalid token.";
      break;
    case XML_ERROR_UNCLOSED_TOKEN:
      cerr << "XML unclosed token.";
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
    case XML_ERROR_BAD_CHAR_REF:
      cerr << "XML bad character reference.";
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

  if(msms_pipeline_analysis.size()==0){
    cerr << "PepXML file contains no MS/MS results." << endl;
    return false;
  }
  uiPipelines.set(&msms_pipeline_analysis);
  setRunSummaries(0);
  calcSize();

  /*
  fileFull = fn;
  filePath = fileFull;
  if (filePath.find_last_of("\\") != string::npos) filePath = filePath.substr(0, filePath.find_last_of("\\"));
  else if (filePath.find_last_of("/") != string::npos) filePath = filePath.substr(0, filePath.find_last_of("/"));
  else filePath.clear();
  fileBase = fileFull;
  if (fileBase.find_last_of("\\") != string::npos) fileBase = fileBase.substr(fileBase.find_last_of("\\") + 1, fileBase.size());
  else if (fileBase.find_last_of("/") != string::npos) fileBase = fileBase.substr(fileBase.find_last_of("/") + 1, fileBase.size());
  if (fileBase.find_last_of(".") != string::npos) fileBase = fileBase.substr(0, fileBase.find_last_of("."));
  */
  return true;
}

void NeoPepXMLParser::setFilterProbability(double probability){
  probFilter=probability;
}

void NeoPepXMLParser::setFilterRunSummary(string str) {
  rsFilter = str;
}

void NeoPepXMLParser::setFilterSearchHit(string str) {
  shFilter = str;
}

bool NeoPepXMLParser::setRunSummaries(const size_t pipeIndex){
  uiRunSummaries.set(NULL,0);
  if(pipeIndex>=msms_pipeline_analysis.size()) return false;
  uiRunSummaries.set(&msms_pipeline_analysis[pipeIndex].msms_run_summary,pipeIndex);
  if(uiRunSummaries.size()>0) setSpectra(pipeIndex,0);
  else uiSpectra.set(NULL,0,0);
  return true;
}

bool NeoPepXMLParser::setSpectra(const size_t pipeIndex, const size_t runIndex){
  uiSpectra.set(NULL,0,0);
  if (pipeIndex >= msms_pipeline_analysis.size()) return false;
  if(runIndex>= msms_pipeline_analysis[pipeIndex].msms_run_summary.size()) return false;
  uiSpectra.set(&msms_pipeline_analysis[pipeIndex].msms_run_summary[runIndex].spectrum_query,pipeIndex,runIndex);
  return true;
}

//returns the number of top ranked PSMs
size_t NeoPepXMLParser::size(){
  return sz;
}

string NeoPepXMLParser::versionNeo(){
  string s;
  s=NPX_VERSION;
  s+="\t";
  s+=NPX_DATE;
  return s;
}

bool NeoPepXMLParser::write(const char* fn, bool tabs){
  FILE* f = fopen(fn, "wt");
  if (f == NULL) return false;

  size_t i;

  fprintf(f, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  for(i=0;i<msms_pipeline_analysis.size();i++){
    if (tabs) msms_pipeline_analysis[i].write(f,0);
    else msms_pipeline_analysis[i].write(f);
  }
  fclose(f);
  return true;

}
