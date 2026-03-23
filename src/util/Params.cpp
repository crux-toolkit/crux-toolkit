#include "AminoAcidUtil.h"
#include "crux_version.h"
#include "model/Modification.h"
#include "parameter.h"
#include "Params.h"
#include "StringUtils.h"

#include <algorithm>

using namespace std;

/**
 * \file Params.cpp
 *
 * This module handles user-level parameters in Crux.  Required
 * arguments are given on the command line; optional parameters can be
 * specified either on the command line or in a parameter file.
 * Textual descriptions of each parameter are stored in the source
 * code, and these are used to automatically generate usage
 * statements, comments in output parameter files, and HTML
 * documentation.
 *
 * Following are the steps required to add a new argument or parameter
 * to an existing command:
 *
 * - Inside the appropriate Params*Options.cpp file (e.g., ParamsCometOptions.cpp
 *   for Comet-specific parameters, ParamsGeneralOptions.cpp for general options,
 *   etc.), add a call to Init{Bool,Int,String,Double,Arg}Param.  This will specify
 *   things like the parameter name, type, default value, what
 *   commands it works with, and whether it is visible to the end
 *   user.
 *
 * - If your parameter is optional, then in Params::Categorize(), add
 *   an "items.insert()" with the name of your new parameter. This
 *   controls what category the parameter gets printed in when the
 *   HTML documentation is created.  Note that if the parameter
 *   doesn't get added to a category in Params::Categorize(), it just
 *   appears under a generic category called "<command name> options"
 *
 * - In the .cpp file that contains the main() for the command that
 *   uses this parameter, find a call to either getArgs() (if your new
 *   parameter is required) or getOptions() (if it is optional), and
 *   add the name of your new parameter to the list of options there.

 * - In the same file, add a call to
 *   Params::Get{Bool,Int,String,Double}() to retrieve the value
 *   associated with the parameter.  In general, these methods can be
 *   used anywhere in the source code in order to retrieve parameters.
 *   However, it's good form, when feasible, to access parameters in
 *   the main() and then to pass them as arguments, rather than
 *   accessing them as globals within subroutines.
 *
 * - If you need to edit the textual description of the command
 *   itself, search in the same file for a call to getDescription().
 *   Portions of that description enclosed in "[[nohtml: XXX ]]" will
 *   be used for the usage statement, and portions in "[[html: XXX ]]"
 *   will be in the HTML docs.
 *
 */

static Params paramContainer_;

Params::Params() : finalized_(false) {
  /* generate_peptide arguments */
  InitArgParam("protein fasta file",
    "The name of the file in FASTA format from which to retrieve proteins.");
  InitArgParam("index name",
    "The desired name of the binary index.");
  InitArgParam("ms2 file",
    "The name of one or more files from which to parse the fragmentation "
    "spectra, in any of the file formats supported by ProteoWizard.");
  /* psm-convert arguments */
  InitArgParam("input PSM file",
    "The name of a PSM file in tab-delimited text, SQT, pepXML or mzIdentML format");
  InitArgParam("output format",
    "The desired format of the output file. Legal values are tsv, html, sqt, pin, "
    "pepxml, mzidentml.");
  /* get-ms2-spectrum */
  InitArgParam("scan number",
    "Scan number identifying the spectrum.");
  InitArgParam("output file",
    "File where spectrum will be written.");
  /* predict-peptide-ions */
  InitArgParam("peptide sequence",
    "The peptide sequence.");
  InitArgParam("charge state",
    "The charge state of the peptide.");
  /* hardklor arguments */
  InitArgParam("spectra",
    "The name of a file from which to parse high-resolution spectra. The file "
    "may be in MS1 (.ms1), binary MS1 (.bms1), compressed MS1 (.cms1), or mzXML "
    "(.mzXML) format.");
  /*Percolator arguments*/
  InitArgParam("peptide-spectrum matches",
    "One or more collections of target and decoy peptide-spectrum matches (PSMs). Input may "
    "be in one of four formats: PIN, SQT, pepXML, or [[html:<a href=\"../file-formats/txt-format.html\">]]"
    "Crux tab-delimited text[[html:</a>]]. "
    "Note that if the input is provided as SQT, pepXML, or Crux "
    "tab-delimited text, then a PIN file will be generated in the output directory "
    "prior to execution. "
    "Crux determines the format of the input file by examining its "
    "filename extension."
    "[[html:<br>For PIN files, target and decoy PSMs are assumed to appear in the "
    "same file.  For other file types, decoy PSMs can be provided to Percolator in two "
    "ways: either as a separate file or embedded within the same file as the target "
    "PSMs. Percolator will first search for target PSMs in a separate file. The "
    "decoy file name is constructed from the target name by replacing \"target\" with "
    "\"decoy\". For example, if search.target.txt is provided as input, then "
    "Percolator will search for a corresponding file named search.decoy.txt. If no "
    "decoy file is found, then Percolator will assume that the given input file "
    "contains a mix of target and decoy PSMs. Within this file, decoys are identified "
    "using a prefix (specified via --decoy-prefix) on the protein name.]]");
  /*make-pin arguments*/
  InitIntParam("max-charge-feature", 0, 0, BILLION,
    "Specifies the maximum charge state feature.  When set to zero, use the "
    "maximum observed charge state.",
    "Available for make-pin and percolator.", true);
  InitArgParam("psm results",
    "A collection of target and decoy peptide-spectrum matches (PSMs). Input may be in "
    "one of four formats: SQT, PepXML (obtained from SEQUEST), [[html:<a href=\""
    "../file-formats/txt-format.html\">]]Crux tab-delimited text[[html:</a>]], or list of files (when "
    "list-of-files=T)."
    "[[html:<br>Decoy PSMs can be provided to make-pin in two ways: either as a separate "
    "file or embedded within the same file as the target PSMs. make-pin will first search "
    "for the target PSMs in a separate file. The decoy file name is constructed from the "
    "target name by replacing \"target\" with \"decoy\". For example, if search.target.txt "
    "is provided as input, then make-pin will search for a corresponding file named "
    "search.decoy.txt. If no decoy file is found, then make-pin will assume that the "
    "given input file contains a mix of target and decoy PSMs. Within this file, decoys "
    "are identified using a prefix (specified via --decoy-prefix) on the protein name.]]");
  InitStringParam("decoy input", "",
    "make-pin can convert any file format in sqt, tab-delimited and pep.xml file "
    "to pin file ",
    "Argument, not option for make-pin", false);
  InitStringParam("output-file", "",
    "Path where pin file will be written instead of make-pin.pin.",
    "It is optional for make-pin", true);
  InitBoolParam("filestem-prefixes", false,
    "Prefix PSM IDs with filestems instead of target or decoy and file index.",
    "Available for make-pin", false);
  InitBoolParam("mod-symbols", false,
    "Print modification symbols instead of masses in peptide sequences.",
    "Available for make-pin", false);
  // pipeline arguments
  InitArgParam("mass spectra",
    "The name of the file(s) from which to parse the fragmentation spectra, in any of the "
    "[[html:<a href=\"http://proteowizard.sourceforge.net/tools.shtml\">]]file formats "
    "supported by ProteoWizard[[html:</a>]]. Alteratively, with Tide-search, these files "
    "may be binary spectrum files produced by a previous run of [[html:<code>]]crux "
    "tide-search[[html:</code>]] using the [[html:<code>]]store-spectra[[html:</code>]] "
    "parameter. Multiple files can be included on the command line (space delimited), "
    "prior to the name of the database.");
  InitArgParam("peptide source",
    "Either the name of a file in fasta format from which to retrieve proteins and "
    "peptides or an index created by a previous run of [[html:<code>]]crux tide-index"
    "[[html:</code>]] (for Tide searching).");
  InitGeneralOptions();
  InitCometOptions();
  InitAnalysisOptions();
  InitBullseyeOptions();
  InitDIAmeterOptions();
  Categorize();
}

Params::~Params() {
  for (map<string, Param*>::iterator i = params_.begin(); i != params_.end(); i++) {
     delete i->second;
  }
}

void Params::Categorize() {
  set<string> items;

  items.clear();
  items.insert("bullseye-max-mass");
  items.insert("bullseye-min-mass");
  items.insert("gap-tolerance");
  items.insert("max-persist");
  items.insert("persist-tolerance");
  items.insert("scan-tolerance");
  AddCategory("Identifying PPIDs in MS1 spectra", items);

  items.clear();
  items.insert("exact-match");
  items.insert("exact-tolerance");
  items.insert("retention-tolerance");
  AddCategory("Matching PPIDs to MS2 spectra", items);

  items.clear();
  items.insert("clip-nterm-methionine");
  items.insert("isotopic-mass");
  items.insert("max-length");
  items.insert("max-mass");
  items.insert("min-length");
  items.insert("min-mass");
  AddCategory("Peptide properties", items);

  items.clear();
  items.insert("cterm-peptide-mods-spec");
  items.insert("cterm-protein-mods-spec");
  items.insert("max-mods");
  items.insert("min-mods");
  items.insert("mod-precision");
  items.insert("mods-spec");
  items.insert("nterm-peptide-mods-spec");
  items.insert("nterm-protein-mods-spec");
  for (char c = 'A'; c <= 'Z'; c++) {
    items.insert(string(1, c));
  }
  items.insert("auto-modifications");
  items.insert("fixed_modification");
  items.insert("fixed_modification_protC");
  items.insert("fixed_modification_protN");
  items.insert("max_mods_per_peptide");
  items.insert("modification");
  items.insert("modification_protC");
  items.insert("modification_protN");
  AddCategory("Amino acid modifications", items);

  items.clear();
  items.insert("allow-dups");
  items.insert("decoy-format");
  items.insert("num-decoys-per-target");
  items.insert("keep-terminal-aminos");
  items.insert("seed");
  AddCategory("Decoy database generation", items);

  items.clear();
  items.insert("custom-enzyme");
  items.insert("digestion");
  items.insert("enzyme");
  items.insert("missed-cleavages");
  items.insert("kojak_enzyme");
  items.insert("max_miscleavages");
  AddCategory("Enzymatic digestion", items);

  items.clear();
  items.insert("auto-precursor-window");
  items.insert("max-precursor-charge");
  items.insert("min-precursor-charge");
  items.insert("precursor-window");
  items.insert("precursor-window-type");
  AddCategory("Precursor selection", items);

  items.clear();
  items.insert("auto-mz-bin-width");
  items.insert("compute-sp");
  items.insert("deisotope");
  items.insert("exact-p-value");
  items.insert("fragment-mass");
  items.insert("isotope-error");
  items.insert("max-ion-charge");
  items.insert("min-peaks");
  items.insert("mod-mass-format");
  items.insert("mz-bin-offset");
  items.insert("mz-bin-width");
  items.insert("peptide-centric-search");
  items.insert("remove-precursor-peak");
  items.insert("remove-precursor-tolerance");
  items.insert("scan-number");
  items.insert("skip-preprocessing");
  items.insert("spectrum-charge");
  items.insert("spectrum-max-mz");
  items.insert("spectrum-min-mz");
  items.insert("use-flanking-peaks");
  items.insert("use-neutral-loss-peaks");
  items.insert("score-function");
  items.insert("fragment-tolerance");
  items.insert("evidence-granularity");
  items.insert("top_count");
  items.insert("e_value_depth");
  items.insert("override-charges");
  AddCategory("Search parameters", items);

  items.clear();
  items.insert("picked-protein");
  items.insert("protein-enzyme");
  items.insert("protein-report-duplicates");
  items.insert("protein-report-fragments");
  AddCategory("Protein inference options", items);

  items.clear();
  items.insert("protein");
  items.insert("fido-alpha");
  items.insert("fido-beta");
  items.insert("fido-empirical-protein-q");
  items.insert("fido-fast-gridsearch");
  items.insert("fido-gamma");
  items.insert("fido-gridsearch-depth");
  items.insert("fido-gridsearch-mse-threshold");
  items.insert("fido-no-split-large-components");
  items.insert("fido-protein-truncation-threshold");
  AddCategory("Fido options", items);

  items.clear();
  // Kojak
  items.insert("cross_link");
  items.insert("mono_link");
  items.insert("mono_links_on_xl");
  items.insert("diff_mods_on_xl");
  AddCategory("Cross-linking parameters", items);

  items.clear();
  items.insert("decoy_search");
  items.insert("peff_format");
  items.insert("peff_obo");
  items.insert("decoy_filter");
  items.insert("max_peptide_mass");
  items.insert("min_peptide_mass");
  items.insert("truncate_prot_names");
  AddCategory("Database", items);

  items.clear();
  items.insert("num-threads");
  items.insert("num_threads");
  items.insert("threads");
  AddCategory("CPU threads", items);

  items.clear();
  items.insert("auto_peptide_mass_tolerance");
  items.insert("isotope_error");
  items.insert("mass_type_fragment");
  items.insert("mass_type_parent");
  items.insert("peptide_mass_tolerance");
  items.insert("peptide_mass_tolerance_lower");
  items.insert("peptide_mass_tolerance_upper");
  items.insert("peptide_mass_units");
  items.insert("precursor_tolerance_type");
  items.insert("ppm_tolerance_pre");
  items.insert("auto_ppm_tolerance_pre");
  items.insert("kojak_isotope_error");
  AddCategory("Masses", items);

  items.clear();
  items.insert("allowed_missed_cleavage");
  items.insert("num_enzyme_termini");
  items.insert("search_enzyme_number");
  items.insert("search_enzyme2_number");
  AddCategory("Search enzyme", items);

  items.clear();
  items.insert("auto_fragment_bin_tol");
  items.insert("fragment_bin_offset");
  items.insert("fragment_bin_tol");
  items.insert("theoretical_fragment_ions");
  items.insert("use_A_ions");
  items.insert("use_B_ions");
  items.insert("use_C_ions");
  items.insert("use_X_ions");
  items.insert("use_Y_ions");
  items.insert("use_Z_ions");
  items.insert("use_Z1_ions");
  items.insert("use_NL_ions");
  items.insert("auto_fragment_bin_size");
  items.insert("fragment_bin_size");
  items.insert("ion_series_A");
  items.insert("ion_series_B");
  items.insert("ion_series_C");
  items.insert("ion_series_X");
  items.insert("ion_series_Y");
  items.insert("ion_series_Z");
  AddCategory("Fragment ions", items);

  items.clear();
  items.insert("activation_method");
  items.insert("ms_level");
  items.insert("override_charge");
  items.insert("precursor_charge");
  items.insert("scan_range");
  AddCategory("mzXML/mzML parameters", items);

  items.clear();
  items.insert("clip_nterm_methionine");
  items.insert("decoy_prefix");
  items.insert("digest_mass_range");
  items.insert("equal_I_and_L");
  items.insert("mass_offsets");
  items.insert("max_duplicate_proteins");
  items.insert("max_fragment_charge");
  items.insert("max_index_runtime");
  items.insert("max_precursor_charge");
  items.insert("nucleotide_reading_frame");
  items.insert("num_results");
  items.insert("output_suffix");
  items.insert("peff_verbose_output");
  items.insert("peptide_length_range");
  items.insert("precursor_NL_ions");
  items.insert("skip_researching");
  items.insert("spectrum_batch_size");
  items.insert("text_file_extension");
  items.insert("explicit_deltacn");
  items.insert("old_mods_encoding");
  items.insert("resolve_fullpaths");
  items.insert("pinfile_protein_delimiter");
  AddCategory("Miscellaneous parameters", items);

  items.clear();
  items.insert("clear_mz_range");
  items.insert("minimum_intensity");
  items.insert("minimum_peaks");
  items.insert("remove_precursor_peak");
  items.insert("remove_precursor_tolerance");
  items.insert("kojak_instrument");
  items.insert("enrichment");
  items.insert("MS1_centroid");
  items.insert("MS2_centroid");
  items.insert("MS1_resolution");
  items.insert("MS2_resolution");
  items.insert("min_spectrum_peaks");
  items.insert("max_spectrum_peaks");
  items.insert("precursor_refinement");
  items.insert("prefer_precursor_pred");
  items.insert("spectrum_processing");
  AddCategory("Spectral processing", items);

  items.clear();
  for (int i = 1; i <= 9; i++) {
    items.insert("variable_mod0" + StringUtils::ToString(i));
  }
  for (int i = 10; i <= 15; i++) {
    items.insert("variable_mod" + StringUtils::ToString(i));
  }
  items.insert("auto_modifications");
  items.insert("max_variable_mods_in_peptide");
  items.insert("require_variable_mod");
  items.insert("protein_modlist_file");
  AddCategory("Variable modifications", items);

  items.clear();
  items.insert("add_Cterm_peptide");
  items.insert("add_Nterm_peptide");
  items.insert("add_Cterm_protein");
  items.insert("add_Nterm_protein");
  for (char c = 'A'; c <= 'Z'; c++) {
    string aaString = string(1, c);
    string aaName = AminoAcidUtil::GetName(c);
    aaName = aaName.empty() ? "user_amino_acid" : StringUtils::Replace(aaName, " ", "_");
    items.insert("add_" + aaString + "_" + aaName);
  }
  AddCategory("Static modifications", items);

  items.clear();
  items.insert("create_peptide_index");
  items.insert("create_fragment_index");
  items.insert("fragindex_max_fragmentmass");
  items.insert("fragindex_min_fragmentmass");
  items.insert("fragindex_min_ions_report");
  items.insert("fragindex_min_ions_score");
  items.insert("fragindex_num_spectrumpeaks");
  items.insert("fragindex_skipreadprecursors");
  AddCategory("Indexing", items);

  items.clear();
  items.insert("only-psms");
  items.insert("tdc");
  items.insert("search-input");
  AddCategory("General options", items);

  items.clear();
  items.insert("c-neg");
  items.insert("c-pos");
  items.insert("maxiter");
  items.insert("percolator-seed");
  items.insert("quick-validation");
  items.insert("subset-max-train");
  items.insert("test-each-iteration");
  items.insert("test-fdr");
  items.insert("train-fdr");
  items.insert("static");
  AddCategory("SVM training options", items);

  items.clear();
  items.insert("default-direction");
  items.insert("init-weights");
  items.insert("klammer");
  items.insert("output-weights");
  items.insert("override");
  items.insert("unitnorm");
  AddCategory("SVM feature input options", items);

  items.clear();
  items.insert("pm-charges");
  items.insert("pm-max-frag-mz");
  items.insert("pm-max-precursor-delta-ppm");
  items.insert("pm-max-precursor-mz");
  items.insert("pm-max-scan-separation");
  items.insert("pm-min-common-frag-peaks");
  items.insert("pm-min-frag-mz");
  items.insert("pm-min-peak-pairs");
  items.insert("pm-min-precursor-mz");
  items.insert("pm-min-scan-frag-peaks");
  items.insert("pm-pair-top-n-frag-peaks");
  items.insert("pm-top-n-frag-peaks");
  AddCategory("param-medic options", items);

  items.clear();
  items.insert("concat");
  items.insert("decoy-prefix");
  items.insert("decoy-xml-output");
  items.insert("feature-file-out");
  items.insert("file-column");
  items.insert("fileroot");
  items.insert("brief-output");
  items.insert("list-of-files");
  items.insert("mass-precision");
  items.insert("mzid-output");
  items.insert("num_output_lines");
  items.insert("output-dir");
  items.insert("output-file");
  items.insert("export_additional_pepxml_scores");
  items.insert("output_mzidentmlfile");
  items.insert("output_pepxmlfile");
  items.insert("output_percolatorfile");
  items.insert("output_sqtstream");
  items.insert("output_sqtfile");
  items.insert("output_txtfile");
  items.insert("overwrite");
  items.insert("parameter-file");
  items.insert("peptide-list");
  items.insert("pepxml-output");
  items.insert("pin-output");
  items.insert("mztab-output");
  items.insert("pout-output");
  items.insert("precision");
  items.insert("print-search-progress");
  items.insert("print_expect_score");
  items.insert("sample_enzyme_number");
  items.insert("show_fragment_ions");
  items.insert("spectrum-format");
  items.insert("spectrum-parser");
  items.insert("sqt-output");
  items.insert("store-index");
  items.insert("store-spectra");
  items.insert("temp-dir");
  items.insert("top-match");
  items.insert("txt-output");
  items.insert("use-z-line");
  items.insert("verbosity");
  items.insert("export_percolator");
  items.insert("export_pepXML");
  items.insert("export_mzID");
  AddCategory("Input and output", items);

}

bool Params::GetBool(const string& name) {
  return Require(name)->GetBool();
}

int Params::GetInt(const string& name) {
  return Require(name)->GetInt();
}

double Params::GetDouble(const string& name) {
  return Require(name)->GetDouble();
}

string Params::GetString(const string& name) {
  return Require(name)->GetString();
}

bool Params::GetBoolDefault(const string& name) {
  return Require(name)->GetBoolDefault();
}

int Params::GetIntDefault(const string& name) {
  return Require(name)->GetIntDefault();
}

double Params::GetDoubleDefault(const string& name) {
  return Require(name)->GetDoubleDefault();
}

string Params::GetStringDefault(const string& name) {
  return Require(name)->GetStringDefault();
}

const vector<string>& Params::GetStrings(const string& name) {
  Param* param = Require(name);
  if (!param->IsArgument()) {
    throw runtime_error("Parameter '" + name + "' is not an argument");
  }
  return ((ArgParam*)param)->GetStrings();
}

string Params::GetUsage(const string& name) {
  return Require(name)->GetUsage();
}

string Params::GetFileNotes(const string& name) {
  return Require(name)->GetFileNotes();
}

bool Params::IsVisible(const string& name) {
  return Require(name)->IsVisible();
}

bool Params::IsArgument(const string& name) {
  return Require(name)->IsArgument();
}

string Params::GetAcceptedValues(const string& name) {
  return Require(name)->GetAcceptedValues();
}

bool Params::IsDefault(const string& name) {
  return Require(name)->IsDefault();
}

bool Params::Exists(const string& name) {
  return paramContainer_.Get(name) != NULL;
}

void Params::Set(const string& name, bool value) {
  paramContainer_.CanModifyCheck();
  Param* param = Require(name);
  param->Set(value);
  param->ThrowIfInvalid();
}

void Params::Set(const string& name, int value) {
  paramContainer_.CanModifyCheck();
  Param* param = Require(name);
  param->Set(value);
  param->ThrowIfInvalid();
}

void Params::Set(const string& name, double value) {
  paramContainer_.CanModifyCheck();
  Param* param = Require(name);
  param->Set(value);
  param->ThrowIfInvalid();
}

void Params::Set(const string& name, const char* value) {
  Set(name, string(value));
}

void Params::Set(const string& name, const string& value) {
  paramContainer_.CanModifyCheck();
  Param* param = Require(name);
  param->Set(value);
  param->ThrowIfInvalid();
}

void Params::AddArgValue(const string& name, const string& value) {
  paramContainer_.CanModifyCheck();
  Param* param = Require(name);
  if (!param->IsArgument()) {
    throw runtime_error("Cannot add value to '" + name + "', it is not an argument");
  }
  ((ArgParam*)param)->AddValue(value);
}

void Params::Finalize() {
  paramContainer_.FinalizeParams();
}

void Params::Write(ostream* out, bool defaults) {
  if (out == NULL || !out->good()) {
    throw runtime_error("Bad file stream for writing parameter file");
  }

  *out << "# Crux parameter file (generated by Crux version " << CRUX_VERSION << ")" << endl
       << "# Full documentation available at http://cruxtoolkit.sourceforge.net/" << endl
       << "# comet_version 2016.01 rev. 1" << endl
       << "# Everything following the \'#\' symbol is treated as a comment." << endl
       << endl;

  for (vector<const Param*>::const_iterator i = Begin(); i != End(); i++) {
    string name = (*i)->GetName();
    // Print mods and Comet parameters later
    if (!(*i)->IsVisible() || name.find('_') != string::npos) {
      continue;
    }
    *out << (*i)->GetParamFileString(defaults) << endl;
  }


  // Print Comet parameters
  *out << "####################" << endl
       << "# Comet Parameters #" << endl
       << "####################" << endl;
  for (vector<const Param*>::const_iterator i = Begin(); i != End(); i++) {
    string name = (*i)->GetName();
    // Print mods and Comet parameters later
    if (!(*i)->IsVisible() || name.find('_') == string::npos) {
      continue;
    }
    *out << (*i)->GetParamFileString(defaults) << endl;
  }

  *out << "#" << endl
       << "# COMET_ENZYME_INFO _must_ be at the end of this parameters file" << endl
       << "#" << endl
       << "[COMET_ENZYME_INFO]" << endl;

  const vector<string>& cometEnzymes = get_comet_enzyme_info_lines();
  if (cometEnzymes.empty() || defaults) {
    *out << "0.  No_enzyme                      0  -          -" << endl
         << "1.  Trypsin                        1  KR         P" << endl
         << "2.  Trypsin/P                      1  KR         -" << endl
         << "3.  Lys_C                          1  K          P" << endl
         << "4.  Lys_N                          0  K          -" << endl
         << "5.  Arg_C                          1  R          P" << endl
         << "6.  Asp_N                          0  D          -" << endl
         << "7.  CNBr                           1  M          -" << endl
         << "8.  Glu_C                          1  DE         P" << endl
         << "9.  PepsinA                        1  FL         P" << endl
         << "10. Chymotrypsin                   1  FWYL       P" << endl;
    /*TODO: Put these back in after we figure out what to do with enzyme info
    *out << "11. Elastase                       1  ALIV       P" << endl
         << "12. Clostripain                    1  R          -" << endl
         << "13. Iodosobenzoate                 1  W          -" << endl
         << "14. Proline_Endopeptidase          1  P          -" << endl
         << "15. Staph_Protease                 1  E          -" << endl
         << "16. Modified_Chymotrypsin          1  FWYL       P" << endl
         << "17. Elastase_Trypsin_Chymotrypsin  1  ALIVKRWFY  P" << endl;
    */
  } else {
    for (vector<string>::const_iterator i = cometEnzymes.begin();
         i != cometEnzymes.end();
         i++) {
      *out << *i << endl;
    }
  }
}

map<string, Param*>::const_iterator Params::BeginAll() {
  return paramContainer_.params_.begin();
}

map<string, Param*>::const_iterator Params::EndAll() {
  return paramContainer_.params_.end();
}

vector<const Param*>::const_iterator Params::Begin() {
  return paramContainer_.paramsOrdered_.begin();
}

vector<const Param*>::const_iterator Params::End() {
  return paramContainer_.paramsOrdered_.end();
}

string Params::ProcessHtmlDocTags(string s, bool html) {
  // If html is true, instances of [[html:{text}]] become {text} and
  // instances of [[nohtml:{text}]] are removed.
  // If html is false, instances of [[nohtml:{text}]] become {text} and
  // instances of [[html:{text}]] are removed.

  const string OPEN_TAG = "[[";
  const string CLOSE_TAG = "]]";
  const string HTML_PREFIX = "html:";
  const string NO_HTML_PREFIX = "nohtml:";

  size_t pos, endPos = 0;
  while ((pos = s.find(OPEN_TAG, endPos)) != string::npos) {
    size_t prefixStart = pos + OPEN_TAG.length();
    if ((endPos = s.find(CLOSE_TAG, prefixStart)) == string::npos) {
      return s;
    }

    string fullOpen = OPEN_TAG;
    bool fullRemove;
    if (s.length() >= prefixStart + HTML_PREFIX.length() &&
        s.compare(prefixStart, HTML_PREFIX.length(), HTML_PREFIX) == 0) {
      fullOpen += HTML_PREFIX;
      fullRemove = !html;
    } else if (s.length() >= prefixStart + NO_HTML_PREFIX.length() &&
               s.compare(prefixStart, NO_HTML_PREFIX.length(), NO_HTML_PREFIX) == 0) {
      fullOpen += NO_HTML_PREFIX;
      fullRemove = html;
    } else {
      endPos = prefixStart;
      continue;
    }

    if (!fullRemove) {
      s.erase(pos, fullOpen.length());
      endPos -= fullOpen.length();
      s.erase(endPos, CLOSE_TAG.length());
    } else {
      s.erase(pos, endPos + CLOSE_TAG.length() - pos);
      endPos = pos;
    }
  }
  return s;
}

vector< pair< string, vector<string> > > Params::GroupByCategory(const vector<string>& options) {
  vector< pair< string, vector<string> > > groups;

  pair< string, vector<string> > uncategorizedPair = make_pair("", vector<string>(options));
  vector<string>& uncategorized = uncategorizedPair.second;

  // Iterate over all categories
  for (vector<ParamCategory>::const_iterator i = paramContainer_.categories_.begin();
       i != paramContainer_.categories_.end();
       i++) {
    bool any = false;
    // Iterate over each given option and check if it is in the category
    for (vector<string>::const_iterator j = options.begin(); j != options.end(); j++) {
      Param* param = Require(*j);
      // This option was in the category
      if (i->Items.find(param) != i->Items.end()) {
        if (!any) {
          any = true;
          groups.push_back(make_pair(i->Name, vector<string>()));
        }
        groups.back().second.push_back(*j);
        vector<string>::iterator iter;
        while ((iter = find(uncategorized.begin(), uncategorized.end(), *j)) !=
               uncategorized.end()) {
          uncategorized.erase(iter);
        }
      }
    }
  }

  if (!uncategorized.empty()) {
    groups.insert(groups.begin(), uncategorizedPair);
  }
  return groups;
}

void Params::InitBoolParam(
  const string& name,
  bool value,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  Add(new BoolParam(name, usage, fileNotes, visible, value));
}

void Params::InitIntParam(
  const string& name,
  int value,
  int min,
  int max,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  Add(new IntParam(name, usage, fileNotes, visible, value, min, max));
}

void Params::InitIntParam(
  const string& name,
  int value,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  Add(new IntParam(name, usage, fileNotes, visible, value));
}

void Params::InitDoubleParam(
  const string& name,
  double value,
  double min,
  double max,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  Add(new DoubleParam(name, usage, fileNotes, visible, value, min, max));
}

void Params::InitDoubleParam(
  const string& name,
  double value,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  Add(new DoubleParam(name, usage, fileNotes, visible, value));
}

void Params::InitStringParam(
  const string& name,
  const string& value,
  const string& validValues,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  Add(new StringParam(name, usage, fileNotes, visible, value,
                      StringUtils::Split(validValues, '|')));
}

void Params::InitStringParam(
  const string& name,
  const string& value,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  Add(new StringParam(name, usage, fileNotes, visible, value));
}

void Params::InitArgParam(
  const string& name,
  const string& usage
) {
  Add(new ArgParam(name, usage));
}

Param* Params::Require(const string& name) {
  Param* param = paramContainer_.Get(name);
  if (param == NULL) {
    throw runtime_error("Parameter '" + name + "' does not exist");
  }
  return param;
}

Param* Params::Get(const string& name) {
  map<string, Param*>::iterator i = params_.find(name);
  return (i == params_.end()) ? NULL : i->second;
}

void Params::Add(Param* param) {
  CanModifyCheck();
  param->ThrowIfInvalid();

  string paramName = param->GetName();
  if (!params_.insert(make_pair(paramName, param)).second) {
    throw runtime_error("Parameter '" + paramName + "' already exists");
  }
  if (!param->IsArgument()) {
    paramsOrdered_.push_back(param);
  }
}

void Params::AddCategory(const string& name, const set<string>& params) {
  // Validate passed in set
  for (set<string>::const_iterator i = params.begin(); i != params.end(); i++) {
    if (Get(*i) == NULL) {
      throw runtime_error("Parameter '" + *i + "' does not exist");
    }
  }

  ParamCategory* category = NULL;
  // Check if this category already exists
  for (vector<ParamCategory>::iterator i = categories_.begin(); i != categories_.end(); i++) {
    if (i->Name == name) {
      category = &*i;
      break;
    }
  }

  // Create new category
  if (!category) {
    categories_.push_back(ParamCategory(name));
    category = &(categories_.back());
  }

  // Loop over parameters and add them to the category if they are in the passed in set
  for (vector<const Param*>::const_iterator i = Begin(); i != End(); i++) {
    string paramName = (*i)->GetName();
    if (params.find(paramName) != params.end()) {
      // Check if this parameter has already been categorized
      for (vector<ParamCategory>::const_iterator j = categories_.begin();
           j != categories_.end();
           j++) {
        if (j->Items.find(*i) != j->Items.end()) {
          throw runtime_error("Parameter '" + paramName + "' has already been categorized");
        }
      }
      // Add parameter to category
      category->Items.insert(*i);
    }
  }
}

void Params::FinalizeParams() {
  if (finalized_) {
    return;
  }

  for (char c = 'A'; c <= 'Z'; c++) {
    string aa = string(1, c);
    double deltaMass = GetDouble(aa);
    if (deltaMass != 0) {
      ModificationDefinition::NewStaticMod(aa, deltaMass, ANY);
    }
  }

  if (GetString("enzyme") == "no-enzyme") {
    Set("digestion", "non-specific-digest");
    Set("missed-cleavages", 500);
  }

  string customEnzyme = GetString("custom-enzyme");
  if (!customEnzyme.empty()) {
    parse_custom_enzyme(customEnzyme);
    Set("enzyme", "custom-enzyme");
  }

  if (GetString("enzyme") == "no-enzyme") {
    Set("digestion", "non-specific-digest");
  } else if (GetString("digestion") == "non-specific-digest") {
    Set("enzyme", "no-enzyme");
  }

  double new_value = GetDouble("mz-bin-width");
// ***************************
#ifdef _MSC_VER
  // Peculiarities of Windows floating point handling
  // results in us getting 0.0 here rather than Nan
  // FIXME: is there a more portable way of checking
  // that a floating point value has not been set?
  if (new_value == 0.0) {
#else
  if (std::isnan(new_value)) {
#endif
    // If no width specified, choose based on mass type.
    if (get_mass_type_parameter("fragment-mass") == MONO) {
      new_value = BIN_WIDTH_MONO;
    } else {
      new_value = BIN_WIDTH_AVERAGE;
    }

    Set("mz-bin-width", new_value);
  }
// ***************************

  finalized_ = true;
}

void Params::CanModifyCheck() const {
  if (finalized_) {
    throw runtime_error("Parameters have been finalized and cannot be modified");
  }
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
