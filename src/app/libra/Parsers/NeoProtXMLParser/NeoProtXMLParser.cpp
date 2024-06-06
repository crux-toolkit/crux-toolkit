#include "NeoProtXMLParser.h"

using namespace std;

// Static callback handlers
static void CMzIdentML_startElementCallback(void *data, const XML_Char *el, const XML_Char **attr) {
  ((NeoProtXMLParser*)data)->startElement(el, attr);
}

static void CMzIdentML_endElementCallback(void *data, const XML_Char *el){
  ((NeoProtXMLParser*)data)->endElement(el);
}

static void CMzIdentML_charactersCallback(void *data, const XML_Char *s, int len){
  ((NeoProtXMLParser*)data)->characters(s, len);
}

NeoProtXMLParser::NeoProtXMLParser(){
  init();
}

NeoProtXMLParser::~NeoProtXMLParser(){
  XML_ParserFree(parser);
}

void NeoProtXMLParser::characters(const XML_Char *s, int len) {

}

void NeoProtXMLParser::endElement(const XML_Char *el) {

  string s;
  for (int i = 0; i<PROTXML_NUM_ELEMENTS; i++){
    if (isElement(elements[i].c_str(), el)){
      if (activeEl.back() != (protXMLElement)i) cout << "Error: unexpected end element: " << elements[i] << " should be " << elements[activeEl.back()] << endl;
      else activeEl.pop_back();
      break;
    }
  }

}

void NeoProtXMLParser::init() {
  parser = XML_ParserCreate(NULL);
  XML_SetUserData(parser, this);
  XML_SetElementHandler(parser, CMzIdentML_startElementCallback, CMzIdentML_endElementCallback);
  XML_SetCharacterDataHandler(parser, CMzIdentML_charactersCallback);

  version = 8;

  elements[prAffectedChannel] = "affected_channel";
  elements[prAnalysisResult] = "analysis_result";
  elements[prAnalysisSummary] = "analysis_summary";
  elements[prAnnotation] = "annotation";
  elements[prContributingChannel] = "contributing_channel";
  elements[prDatasetDerivation] = "dataset_derivation";
  elements[prDecoyAnalysis] = "decoy_analysis";
  elements[prDecoyAnalysisSummary] = "decoy_analysis_summary";
  elements[prErrorPoint] = "error_point";
  elements[prFragmentMasses] = "fragment_masses";
  elements[prIndistinguishablePeptide] = "indistinguishable_peptide";
  elements[prIndistinguishableProtein] = "indistinguishable_protein";
  elements[prIntensity] = "intensity";
  elements[prIsotopicContributions] = "isotopic_contributions";
  elements[prLibraResult] = "libra_result";
  elements[prLibraSummary] = "libra_summary";
  elements[prModAminoacidMass] = "mod_aminoacid_mass";
  elements[prModificationInfo] = "modification_info";
  elements[prNSPDistribution] = "nsp_distribution";
  elements[prNSPInformation] = "nsp_information";
  elements[prParameter] = "parameter";
  elements[prPeptide] = "peptide";
  elements[prPeptideParentProtein] = "peptide_parent_protein";
  elements[prPoint] = "point";
  elements[prProgramDetails] = "program_details";
  elements[prProtein] = "protein";
  elements[prProteinGroup] = "protein_group";
  elements[prProteinProphetDetails] = "proteinprophet_details";
  elements[prProteinSummary] = "protein_summary";
  elements[prProteinSummaryDataFilter] = "protein_summary_data_filter";
  elements[prProteinSummaryHeader] = "protein_summary_header";
  elements[prStPeterAnalysisSummary] = "StPeter_analysis_summary";
  elements[prStPeterQuant] = "StPeterQuant";
  elements[prStPeterQuantPeptide] = "StPeterQuant_peptide";

}

void NeoProtXMLParser::startElement(const XML_Char *el, const XML_Char **attr){

  //cout << el << endl; //for diagnostics

  //string s;
  if (isElement("affected_channel", el)){
    activeEl.push_back(prAffectedChannel);
    CnprAffectedChannel c;
    c.channel = atoi(getAttrValue("channel", attr));
    c.correction = atof(getAttrValue("correction", attr));
    protein_summary.analysis_summary.back().libra_summary.back().isotopic_contributions.back().contributing_channel.back().affected_channel.push_back(c);
  
  } else if (isElement("analysis_result", el)){
    activeEl.push_back(prAnalysisResult);
    CnprAnalysisResult c;
    c.analysis= getAttrValue("analysis", attr);
    c.id = atoi(getAttrValue("id", attr));
    protein_summary.protein_group.back().protein.back().analysis_result.push_back(c);

  } else if (isElement("analysis_summary", el)){
    activeEl.push_back(prAnalysisSummary);
    CnprAnalysisSummary c;
    c.analysis= getAttrValue("analysis", attr);
    c.id = atoi(getAttrValue("id", attr));
    c.time.parseDateTime(getAttrValue("time", attr));
    protein_summary.analysis_summary.push_back(c);

  } else if(isElement("annotation", el)){
    activeEl.push_back(prAnnotation);
    CnprAnnotation c;
    c.ensembl_name = getAttrValue("ensembl_name", attr);
    c.flybase = getAttrValue("flybase", attr);
    c.ipi_name = getAttrValue("ipi_name", attr);
    c.locus_link_name = getAttrValue("locus_link_name", attr);
    c.protein_description = getAttrValue("protein_description", attr);
    c.refseq_name = getAttrValue("refseq_name", attr);
    c.swissprot_name = getAttrValue("swissprot_name", attr);
    c.trembl_name = getAttrValue("trembl_name", attr);
    switch (activeEl[activeEl.size() - 2]){
    case prIndistinguishableProtein:
      protein_summary.protein_group.back().protein.back().indistinguishable_protein.back().annotation.push_back(c);
      break;
    case prProtein:
      protein_summary.protein_group.back().protein.back().annotation.push_back(c);
      break;
    default:
      cerr << "Unhandled annotation element inside " << elements[activeEl[activeEl.size() - 2]] << endl;
      exit(99);
    }

  } else if (isElement("contributing_channel", el)){
    activeEl.push_back(prContributingChannel);
    CnprContributingChannel c;
    c.channel = atoi(getAttrValue("channel", attr));
    protein_summary.analysis_summary.back().libra_summary.back().isotopic_contributions.back().contributing_channel.push_back(c);

  } else if (isElement("dataset_derivation", el)){
    activeEl.push_back(prDatasetDerivation);
    CnprDatasetDerivation c;
    c.generation_no=getAttrValue("generation_no",attr);
    protein_summary.dataset_derivation=c;

  } else if (isElement("decoy_analysis", el)){
    activeEl.push_back(prDecoyAnalysis);
    CnprDecoyAnalysis c;
    protein_summary.analysis_summary.back().decoy_analysis.push_back(c);

  } else if (isElement("decoy_analysis_summary", el)){
    activeEl.push_back(prDecoyAnalysisSummary);
    CnprDecoyAnalysisSummary c;
    c.decoy_ratio=atof(getAttrValue("decoy_ratio",attr));
    c.decoy_string=getAttrValue("decoy_string",attr);
    c.exclude_string = getAttrValue("exclude_string", attr);
    c.use_confidence = getAttrValue("use_confidence", attr);
    protein_summary.analysis_summary.back().decoy_analysis_summary.push_back(c);

  } else if (isElement("error_point", el)){
    activeEl.push_back(prErrorPoint);
    CnprErrorPoint c;
    c.error = atof(getAttrValue("error", attr));
    c.min_prob = atof(getAttrValue("min_prob", attr));
    c.num_corr = atoi(getAttrValue("num_corr", attr));
    c.num_incorr = atoi(getAttrValue("num_incorr", attr));
    protein_summary.protein_summary_header.program_details.proteinprophet_details.back().error_point.push_back(c);

  } else if (isElement("fragment_masses", el)){
    activeEl.push_back(prFragmentMasses);
    CnprFragmentMasses c;
    c.channel = atoi(getAttrValue("channel", attr));
    c.mz = atof(getAttrValue("mz", attr));
    protein_summary.analysis_summary.back().libra_summary.back().fragment_masses.push_back(c);

  } else if (isElement("indistinguishable_peptide", el)){
    activeEl.push_back(prIndistinguishablePeptide);
    CnprIndistinguishablePeptide c;
    c.peptide_sequence = getAttrValue("peptide_sequence", attr);
    c.charge=atoi(getAttrValue("charge",attr));
    c.calc_neutral_pep_mass=atof(getAttrValue("calc_neutral_pep_mass",attr));
    protein_summary.protein_group.back().protein.back().peptide.back().indistinguishable_peptide.push_back(c);

  } else if (isElement("indistinguishable_protein", el)){
    activeEl.push_back(prIndistinguishableProtein);
    CnprIndistinguishableProtein c;
    c.protein_name = getAttrValue("protein_name", attr);
    protein_summary.protein_group.back().protein.back().indistinguishable_protein.push_back(c);

  } else if (isElement("intensity", el)){
    activeEl.push_back(prIntensity);
    CnprIntensity c;
    c.channel = atoi(getAttrValue("channel", attr));
    c.mz = atof(getAttrValue("mz", attr));
    c.ratio = atof(getAttrValue("ratio", attr));
    c.error = atof(getAttrValue("error", attr));
    protein_summary.protein_group.back().protein.back().analysis_result.back().libra_result.back().intensity.push_back(c);

  } else if (isElement("isotopic_contributions", el)){
    activeEl.push_back(prIsotopicContributions);
    CnprIsotopicContributions c;
    protein_summary.analysis_summary.back().libra_summary.back().isotopic_contributions.push_back(c);

  } else if (isElement("libra_result", el)){
    activeEl.push_back(prLibraResult);
    CnprLibraResult c;
    c.number = atoi(getAttrValue("number", attr));
    protein_summary.protein_group.back().protein.back().analysis_result.back().libra_result.push_back(c);

  } else if (isElement("libra_summary", el)){
    activeEl.push_back(prLibraSummary);
    CnprLibraSummary c;
    c.centroiding_preference = atoi(getAttrValue("centroiding_preference", attr));
    c.channel_code = getAttrValue("channel_code", attr);
    c.mass_tolerance = atof(getAttrValue("mass_tolerance", attr));
    c.min_pep_prob = atof(getAttrValue("min_pep_prob", attr));
    c.min_pep_wt = atof(getAttrValue("min_pep_wt", attr));
    c.min_prot_prob = atof(getAttrValue("min_prot_prob", attr));
    c.normalization = atoi(getAttrValue("normalization", attr));
    c.output_type = atoi(getAttrValue("output_type", attr));
    c.version = getAttrValue("version", attr);
    protein_summary.analysis_summary.back().libra_summary.push_back(c);

  } else if (isElement("mod_aminoacid_mass", el)){
    activeEl.push_back(prModAminoacidMass);
    CnprModAminoacidMass c;
    c.position = atoi(getAttrValue("position", attr));
    c.mass = atof(getAttrValue("mass", attr));
    switch (activeEl[activeEl.size() - 3]){
    case prIndistinguishablePeptide:
      protein_summary.protein_group.back().protein.back().peptide.back().indistinguishable_peptide.back().modification_info.back().mod_aminoacid_mass.push_back(c);
      break;
    case prPeptide:
      protein_summary.protein_group.back().protein.back().peptide.back().modification_info.back().mod_aminoacid_mass.push_back(c);
      break;
    default:
      cerr << "Unhandled annotation element inside " << elements[activeEl[activeEl.size() - 2]] << endl;
      exit(99);
    }

  } else if (isElement("modification_info", el)){
    activeEl.push_back(prModificationInfo);
    CnprModificationInfo c;
    c.modified_peptide = getAttrValue("modified_peptide", attr);
    c.mod_cterm_mass = atof(getAttrValue("mod_cterm_mass", attr));
    c.mod_nterm_mass = atof(getAttrValue("mod_nterm_mass", attr));
    switch (activeEl[activeEl.size() - 2]){
    case prIndistinguishablePeptide:
      protein_summary.protein_group.back().protein.back().peptide.back().indistinguishable_peptide.back().modification_info.push_back(c);
      break;
    case prPeptide:
      protein_summary.protein_group.back().protein.back().peptide.back().modification_info.push_back(c);
      break;
    default:
      cerr << "Unhandled annotation element inside " << elements[activeEl[activeEl.size() - 2]] << endl;
      exit(99);
    }

  } else if (isElement("nsp_distribution", el)){
    activeEl.push_back(prNSPDistribution);
    CnprNSPDistribution c;
    c.alt_pos_to_neg_ratio = atof(getAttrValue("alt_pos_to_neg_ratio", attr));
    c.bin_no = atoi(getAttrValue("bin_no", attr));
    c.neg_freq = atof(getAttrValue("neg_freq", attr));
    c.nsp_lower_bound_excl = atof(getAttrValue("nsp_lower_bound_excl", attr));
    c.nsp_lower_bound_incl = atof(getAttrValue("nsp_lower_bound_incl", attr));
    c.nsp_upper_bound_excl = getAttrValue("nsp_upper_bound_excl", attr);
    c.nsp_upper_bound_incl = getAttrValue("nsp_upper_bound_incl", attr);
    c.pos_freq = atof(getAttrValue("pos_freq", attr));
    c.pos_to_neg_ratio = atof(getAttrValue("pos_to_neg_ratio", attr));
    protein_summary.protein_summary_header.program_details.proteinprophet_details.back().nsp_information.nsp_distribution.push_back(c);

  } else if (isElement("nsp_information", el)){
    activeEl.push_back(prNSPInformation);
    CnprNSPInformation c;
    c.neighboring_bin_smoothing = getAttrValue("neighboring_bin_smoothing", attr);
    protein_summary.protein_summary_header.program_details.proteinprophet_details.back().nsp_information=c;

  } else if(isElement("parameter", el)){
    activeEl.push_back(prParameter);
    CnprParameter c;
    c.name = getAttrValue("name", attr);
    c.value = getAttrValue("value", attr);
    c.type = getAttrValue("type", attr);
    switch(activeEl[activeEl.size()-2]){
    case prIndistinguishableProtein:
      protein_summary.protein_group.back().protein.back().indistinguishable_protein.back().parameter.push_back(c);
      break;
    case prPeptide:
      protein_summary.protein_group.back().protein.back().peptide.back().parameter.push_back(c);
      break;
    case prProtein:
      protein_summary.protein_group.back().protein.back().parameter.push_back(c);
      break;
    default:
      cerr << "Unhandled parameter element inside " << elements[activeEl[activeEl.size() - 2]] << endl;
      exit(99);
    }

  } else if (isElement("peptide", el)){
    activeEl.push_back(prPeptide);
    CnprPeptide c;
    c.calc_neutral_pep_mass = atof(getAttrValue("calc_neutral_pep_mass", attr));
    c.charge = atoi(getAttrValue("charge", attr));
    c.exp_sibling_ion_bin = atof(getAttrValue("exp_sibling_ion_bin", attr));
    c.exp_sibling_ion_instances = atof(getAttrValue("exp_sibling_ion_instances", attr));
    c.exp_tot_instances = atof(getAttrValue("exp_tot_instances", attr));
    c.fpkm_adjusted_probability = atof(getAttrValue("fpkm_adjusted_probability", attr));
    c.fpkm_bin = atoi(getAttrValue("fpkm_bin", attr));
    c.initial_probability = atof(getAttrValue("initial_probability", attr));
    c.is_contributing_evidence = getAttrValue("is_contributing_evidence", attr);
    c.is_nondegenerate_evidence = getAttrValue("is_nondegenerate_evidence", attr);
    c.max_fpkm = atof(getAttrValue("max_fpkm", attr));
    c.ni_adjusted_probability = atof(getAttrValue("ni_adjusted_probability", attr));
    c.nsp_adjusted_probability = atof(getAttrValue("nsp_adjusted_probability", attr));
    c.n_enzymatic_termini = atoi(getAttrValue("n_enzymatic_termini", attr));
    c.n_instances = atoi(getAttrValue("n_instances", attr));
    c.n_sibling_peptides = atof(getAttrValue("n_sibling_peptides", attr));
    c.n_sibling_peptides_bin = atoi(getAttrValue("n_sibling_peptides_bin", attr));
    c.peptide_group_designator = getAttrValue("peptide_group_designator", attr);
    c.peptide_sequence = getAttrValue("peptide_sequence", attr);
    c.weight = atof(getAttrValue("weight", attr));
    protein_summary.protein_group.back().protein.back().peptide.push_back(c);

  } else if (isElement("peptide_parent_protein", el)){
    activeEl.push_back(prPeptideParentProtein);
    CnprPeptideParentProtein c;
    c.protein_name = getAttrValue("protein_name", attr);
    protein_summary.protein_group.back().protein.back().peptide.back().peptide_parent_protein.push_back(c);

  } else if (isElement("point", el)){
    activeEl.push_back(prPoint);
    CnprPoint c;
    c.fdr_pp = atof(getAttrValue("fdr_pp", attr));
    c.fdr_pp_decoy = atof(getAttrValue("fdr_pp_decoy", attr));
    c.num_corr_pp = atof(getAttrValue("num_corr_pp", attr));
    c.num_corr_pp_decoy = atof(getAttrValue("num_corr_pp_decoy", attr));
    c.pp_decoy_uncert = atof(getAttrValue("pp_decoy_uncert", attr));
    c.pp_uncert = atof(getAttrValue("pp_uncert", attr));
    c.prob_cutoff = atof(getAttrValue("prob_cutoff", attr));
    protein_summary.analysis_summary.back().decoy_analysis.back().point.push_back(c);

  } else if (isElement("program_details", el)){
    activeEl.push_back(prProgramDetails);
    CnprProgramDetails c;
    c.analysis = getAttrValue("analysis", attr);
    c.time.parseDateTime(getAttrValue("time", attr));
    c.version = getAttrValue("version", attr);
    protein_summary.protein_summary_header.program_details=c;

  } else if (isElement("protein", el)){
    activeEl.push_back(prProtein);
    CnprProtein c;
    c.confidence = atof(getAttrValue("confidence", attr));
    c.group_sibling_id = getAttrValue("group_sibling_id",attr);
    c.n_indistinguishable_proteins = atoi(getAttrValue("n_indistinguishable_proteins", attr));
    c.pct_spectrum_ids = getAttrValue("pct_spectrum_ids", attr);
    c.percent_coverage = atof(getAttrValue("percent_coverage", attr));
    c.probability = atof(getAttrValue("probability", attr));
    c.protein_name = getAttrValue("protein_name", attr);
    c.subsuming_protein_entry = getAttrValue("subsuming_protein_entry", attr);
    c.total_number_distinct_peptides = atoi(getAttrValue("total_number_distinct_peptides", attr));
    c.total_number_peptides = atoi(getAttrValue("total_number_peptides", attr));
    c.unique_stripped_peptides = getAttrValue("unique_stripped_peptides", attr);
    protein_summary.protein_group.back().protein.push_back(c);

  } else  if (isElement("protein_group", el)){
    activeEl.push_back(prProteinGroup);
    CnprProteinGroup c;
    c.group_number = getAttrValue("group_number", attr);
    c.probability = atof(getAttrValue("probability", attr));
    c.pseudo_name = getAttrValue("pseudo_name", attr);
    protein_summary.protein_group.push_back(c);

  } else  if (isElement("proteinprophet_details", el)){
    activeEl.push_back(prProteinProphetDetails);
    CnprProteinProphetDetails c;
    c.degen_flag = getAttrValue("degen_flag", attr);
    c.final_peptide_wt_iters = getAttrValue("final_peptide_wt_iters", attr);
    c.fpkm_flag = getAttrValue("fpkm_flag", attr);
    c.groups_flag = getAttrValue("groups_flag", attr);
    c.initial_peptide_wt_iters = getAttrValue("initial_peptide_wt_iters", attr);
    c.nsp_distribution_iters = getAttrValue("nsp_distribution_iters", attr);
    c.nsp_flag = getAttrValue("nsp_flag", attr);
    c.occam_flag = getAttrValue("occam_flag", attr);
    c.run_options = getAttrValue("run_options", attr);
    protein_summary.protein_summary_header.program_details.proteinprophet_details.push_back(c);

  } else if (isElement("protein_summary", el)){
    activeEl.push_back(prProteinSummary);
    CnprProteinSummary c;
    c.summary_xml=getAttrValue("summary_xml",attr);
    c.xmlns = getAttrValue("xmlns", attr);
    c.xmlns_xsi = getAttrValue("xmlns:xsi", attr);
    c.xsi_schemaLocation = getAttrValue("xsi:schemaLocation", attr);
    protein_summary=c;

  } else if (isElement("protein_summary_data_filter", el)){
    activeEl.push_back(prProteinSummaryDataFilter);
    CnprProteinSummaryDataFilter c;
    c.false_positive_error_rate = atof(getAttrValue("false_positive_error_rate", attr));
    c.min_probability = atof(getAttrValue("min_probability", attr));
    c.predicted_num_correct = atof(getAttrValue("predicted_num_correct", attr));
    c.predicted_num_incorrect = atof(getAttrValue("predicted_num_incorrect", attr));
    c.sensitivity = atof(getAttrValue("sensitivity", attr));
    protein_summary.protein_summary_header.program_details.proteinprophet_details.back().protein_summary_data_filter.push_back(c);

  } else if (isElement("protein_summary_header", el)){
    activeEl.push_back(prProteinSummaryHeader);
    CnprProteinSummaryHeader c;
    c.initial_min_peptide_pro = atof(getAttrValue("initial_min_peptide_pro", attr));
    c.min_peptide_probability = atof(getAttrValue("min_peptide_probability", attr));
    c.min_peptide_weight = atof(getAttrValue("min_peptide_weight", attr));
    c.num_input_1_spectra = atoi(getAttrValue("num_input_1_spectra", attr));
    c.num_input_2_spectra = atoi(getAttrValue("num_input_2_spectra", attr));
    c.num_input_3_spectra = atoi(getAttrValue("num_input_3_spectra", attr));
    c.num_input_4_spectra = atoi(getAttrValue("num_input_4_spectra", attr));
    c.num_input_5_spectra = atoi(getAttrValue("num_input_5_spectra", attr));
    c.num_predicted_correct_prots = atof(getAttrValue("num_predicted_correct_prots", attr));
    c.organism = getAttrValue("organism", attr);
    c.reference_database = getAttrValue("reference_database", attr);
    c.residue_substitution_list = getAttrValue("residue_substitution_list", attr);
    c.sample_enzyme = getAttrValue("sample_enzyme", attr);
    c.source_files = getAttrValue("source_files", attr);
    c.source_files_alt = getAttrValue("source_files_alt", attr);
    c.source_file_xtn = getAttrValue("source_file_xtn", attr);
    c.total_no_spectrum_ids = atof(getAttrValue("total_no_spectrum_ids", attr));
    c.win_cyg_reference_database = getAttrValue("win_cyg_reference_database", attr);
    protein_summary.protein_summary_header=c;

  } else if (isElement("StPeter_analysis_summary", el)){
    activeEl.push_back(prStPeterAnalysisSummary);
    CnprStPeterAnalysisSummary c;
    c.probability = atof(getAttrValue("probability", attr));
    c.FDR = atof(getAttrValue("FDR", attr));
    c.tolerance = atof(getAttrValue("tolerance", attr));
    c.sampleLoad = atof(getAttrValue("sampleLoad", attr));
    c.version = getAttrValue("version", attr);
    c.degenerate_peptides = getAttrValue("degenerate_peptides", attr);
    protein_summary.analysis_summary.back().StPeter_analysis_summary.push_back(c);

  } else if (isElement("StPeterQuant", el)){
    activeEl.push_back(prStPeterQuant);
    CnprStPeterQuant c;
    c.counts = atof(getAttrValue("counts", attr));
    c.dCounts = atof(getAttrValue("dCounts", attr));
    c.dNSAF = atof(getAttrValue("dNSAF", attr));
    c.dSI = atof(getAttrValue("dSI", attr));
    c.dSIn = atof(getAttrValue("dSIn", attr));
    c.ng = atof(getAttrValue("ng", attr));
    c.ngC = atof(getAttrValue("ngC", attr));
    c.NSAF = atof(getAttrValue("NSAF", attr));
    c.SI = atof(getAttrValue("SI", attr));
    c.SIn = atof(getAttrValue("SIn", attr));
    protein_summary.protein_group.back().protein.back().analysis_result.back().StPeterQuant.push_back(c);

  } else if (isElement("StPeterQuant_peptide", el)){
    activeEl.push_back(prStPeterQuantPeptide);
    CnprStPeterQuantPeptide c;
    c.sequence = getAttrValue("sequence", attr);
    c.charge = atoi(getAttrValue("charge", attr));
    c.dSC = atof(getAttrValue("dSC", attr));
    c.dSI = atof(getAttrValue("dSI", attr));
    c.SC = atof(getAttrValue("SC", attr));
    c.SI = atof(getAttrValue("SI", attr));
    protein_summary.protein_group.back().protein.back().analysis_result.back().StPeterQuant.back().StPeterQuant_peptide.push_back(c);

  } else {
    cout << "WARNING: Element undefined: " << el << endl;
    //exit(1);
  }
}

bool NeoProtXMLParser::read(const char* fn){
  XML_ParserFree(parser);
  parser = XML_ParserCreate(NULL);
  XML_SetUserData(parser, this);
  XML_SetElementHandler(parser, CMzIdentML_startElementCallback, CMzIdentML_endElementCallback);
  XML_SetCharacterDataHandler(parser, CMzIdentML_charactersCallback);

  // clear data

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

string NeoProtXMLParser::versionNeo() {
  string s;
  s = NPR_VERSION;
  s += "\t";
  s += NPR_DATE;
  return s;
}

bool NeoProtXMLParser::write(const char* fn, bool tabs){
  FILE* f = fopen(fn, "wt");
  if (f == NULL) return false;

  fprintf(f, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  if(tabs) protein_summary.write(f,0);
  else protein_summary.write(f);

  fclose(f);
  return true;

}
