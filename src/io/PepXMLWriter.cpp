/**
 * \file PepXMLWriter.cpp
 * \brief Writes search results in the TPP-compliant .pep.xml format.
 */
#include <iostream>
#include "PepXMLWriter.h"
#include "util/AminoAcidUtil.h"
#include "util/crux-utils.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "model/MatchCollection.h"

using namespace std;

PepXMLWriter::PepXMLWriter():
  file_(NULL), current_index_(1), mass_precision_(Params::GetInt("mass-precision")),
  enzyme_(get_enzyme_type_parameter("enzyme")),
  precision_(Params::GetInt("precision")) {
}

PepXMLWriter::~PepXMLWriter() {
  closeFile();
}

/**
 * Open a file of the given name.  Replace an existing file if
 * overwrite is true, else exit if an existing file is found.
 */
void PepXMLWriter::openFile(const char* filename, bool overwrite) {
  file_ = create_file_in_path(filename, "", overwrite);
}

/**
 * Close the file, if open.
 */
void PepXMLWriter::closeFile() {
  if (file_) {
    fclose(file_);
    file_ = NULL;
  }
}

/**
 * Write all of the elements that preceed the spectrum_query
 * elements.
 * Requires OpenFile has been called without CloseFile.
 */
void PepXMLWriter::writeHeader() {
  MatchCollection::printXmlHeader(file_, "");
}

/**
 * Close the msms_run_summaryand msms_pipeline_analysis tags.
 * Requires OpenFile has been called without CloseFile.
 */
void PepXMLWriter::writeFooter() {
  closeSpectrumElement();
  MatchCollection::printXmlFooter(file_);
}

/**
 * Write the details for a PSM to be contained in a spectrum_query
 * element.  
 * Begins with the <spectrum_query> element if spectrum_scan_number is
 * different than in the last writePSM() call.  Writes one
 * <search_hit> in the list of <search_elements>.  To add more
 * <search_hits>, call writePSM() again with the same spectrum_scan_number
 *
 * Requires that the arrays pre_aas, post_aas,
 * protein_names, protein_descriptions are all num_proteins long.
 * Assumes scores is NUMBER_SCORER_TYPES long. 
 * Requires OpenFile has been called without CloseFile.
 */
void PepXMLWriter::writePSM(
  int spectrum_scan_number, ///< identifier for the spectrum
  const char* filename, ///< file the spectrum came from
  double spectrum_neutral_mass, ///< computed mass of the spectrum at this charge state
  int charge, ///< assumed charge state for the match
  int* PSM_ranks, ///< rank of this peptide for the spectrum
  const char* unmodified_peptide_sequence, ///< sequence with no mods
  const char* modified_peptide_sequence, ///< either with symbols or masses
  double peptide_mass, ///< mass of the peptide sequence
  int num_proteins, ///< proteins matched to this peptide
  const char* flanking_aas, ///< "XY, AB, " X and Y are the preceeding and
                            /// following aas in the first protein 
  vector<string>& protein_names, ///<
  vector<string>& protein_descriptions, ///<
  bool* scores_computed,
  double* scores, ///< indexed by score type
  unsigned cur_num_matches
) {
  string spectrum_title = getSpectrumTitle(spectrum_scan_number, filename, charge);
  //cerr<<"by_ion_fraction_matched: "<<by_ion_fraction_matched<<endl;
  // close the last spec element if this is a new spectrum and not the first
  if (!last_spectrum_printed_.empty() && last_spectrum_printed_ != spectrum_title) {
    closeSpectrumElement();
  }

  // print the spec info if this is a new spectrum
  if (last_spectrum_printed_ != spectrum_title) {
    printSpectrumElement(spectrum_scan_number, spectrum_title.c_str(), 
                         spectrum_neutral_mass, charge);
    last_spectrum_printed_ = spectrum_title;
  }
  // else, just add to the search_result list
  //cerr<<"3.LnExperimentSize: "<<current_ln_experiment_size<<endl;
  printPeptideElement(PSM_ranks, 
    unmodified_peptide_sequence,
    modified_peptide_sequence, 
    peptide_mass,
    spectrum_neutral_mass,
    num_proteins,
    flanking_aas,
    protein_names,
    protein_descriptions,
    scores_computed,
    scores,
    cur_num_matches);
}

/**
 * Write the <spectrum_query> element and the <search_result> tag.
 */
void PepXMLWriter::printSpectrumElement(int spectrum_scan_number, 
                                        const char* spectrum_title,
                                        double spectrum_neutral_mass, 
                                        int charge) {
  fprintf(file_, "    <spectrum_query spectrum=\"%s\" start_scan=\"%i\""
          " end_scan=\"%i\" precursor_neutral_mass=\"%.*f\""
          " assumed_charge=\"%i\" index=\"%i\">\n"
          "    <search_result>\n",
          spectrum_title,
          spectrum_scan_number,
          spectrum_scan_number,
          Params::GetInt("mass-precision"),
          spectrum_neutral_mass,
          charge,
          current_index_++);
}

string PepXMLWriter::getSpectrumTitle(int spectrum_scan_number, 
                                      const char* filename,
                                      int charge) {
  string formatted(filename);
  if (StringUtils::IEndsWith(formatted, ".pid")) {
    // Remove bullseye extension if it exists
    formatted.erase(formatted.length() - 4);
  }
  std::ostringstream spectrum_id;
  spectrum_id << formatted << "." << std::setw(5)  << std::setfill('0')
              << spectrum_scan_number << "." << std::setw(5) 
              << std::setfill('0') << spectrum_scan_number << "." << charge;
  return spectrum_id.str();
}

void PepXMLWriter::closeSpectrumElement() {
  fprintf(file_, "    </search_result>\n    </spectrum_query>\n");
}

/**
 * Print everything between the <search_hit> tags for this match
 */
void PepXMLWriter::printPeptideElement(int *ranks,
  const char* peptide_sequence,
  const char* modified_peptide_sequence,
  double peptide_mass,
  double spectrum_mass,
  int num_proteins,
  const char* flanking_aas,
  vector<string>& protein_names,
  vector<string>& protein_descriptions,
  bool* scores_computed,
  double* scores,
  unsigned current_num_matches  
) {
  // get values
  char flanking_aas_prev = flanking_aas[0];
  char flanking_aas_next = flanking_aas[1];
  string protein_id = (!protein_names.empty()) ? protein_names.front() : "";
  string protein_annotation = (!protein_descriptions.empty()) ?
    protein_descriptions.front() : "";
  int num_tol_term = get_num_terminal_cleavage(peptide_sequence,
                                               flanking_aas_prev,
                                               flanking_aas_next,
                                               enzyme_);
  int num_missed_cleavages = get_num_internal_cleavage(peptide_sequence,
                                                       enzyme_);

  // print <search_hit> tag
  fprintf(file_, "    <search_hit hit_rank=\"%i\" peptide=\"%s\" "
          "peptide_prev_aa=\"%c\" peptide_next_aa=\"%c\" protein=\"%s\" "
          "num_tot_proteins=\"%i\" ",
          ranks[XCORR], // -1 if unavailable, uses xcorr rank otherwise
          peptide_sequence,
          flanking_aas_prev,
          flanking_aas_next,
          protein_id.c_str(),
          num_proteins);
  if (scores_computed[BY_IONS_MATCHED]) {
    fprintf(file_, "num_matched_ions=\"%i\" ", (unsigned)scores[BY_IONS_MATCHED]);
  }
  if (scores_computed[BY_IONS_TOTAL]) {
    fprintf(file_, "tot_num_ions=\"%i\" ", (unsigned)scores[BY_IONS_TOTAL]);
  }

  if (flanking_aas_prev != 'X' && flanking_aas_next != 'X') {
    fprintf(file_, "calc_neutral_pep_mass=\"%.*f\" "
          "massdiff=\"%+.*f\" "
          "num_tol_term=\"%i\" num_missed_cleavages=\"%i\" "
          "num_matched_peptides=\"%i\""
          " is_rejected=\"%i\" ",
          mass_precision_,
          peptide_mass,
          mass_precision_,
          spectrum_mass - peptide_mass,
          num_tol_term, 
          num_missed_cleavages, 
          current_num_matches,
          0);
  } else {
    fprintf(file_, "calc_neutral_pep_mass=\"%.*f\" "
          "massdiff=\"%+.*f\" "
          "num_tol_term=\"\" num_missed_cleavages=\"%i\" "
          "num_matched_peptides=\"%i\""
          " is_rejected=\"%i\" ",
          mass_precision_,
          peptide_mass,
          mass_precision_,
          spectrum_mass - peptide_mass, 
          num_missed_cleavages, 
          current_num_matches,
          0);
  }

  fprintf(file_, "protein_descr=\"%s\">\n", protein_annotation.c_str());

  // print additonal proteins
  for(int prot_idx = 1; prot_idx < num_proteins; prot_idx++) {
    int flank_idx = strlen("XX,") * prot_idx;
    flanking_aas_prev = flanking_aas[flank_idx];
    flanking_aas_next = flanking_aas[flank_idx + 1];
    num_tol_term = get_num_terminal_cleavage(peptide_sequence,
                                             flanking_aas_prev,
                                             flanking_aas_next,
                                             enzyme_);
    protein_id = protein_names[prot_idx];
    protein_annotation = protein_descriptions[prot_idx];
    fprintf(file_, 
            "        <alternative_protein protein=\"%s\" "
            "protein_descr=\"%s\" "
            "num_tol_term=\"%i\"  peptide_prev_aa=\"%c\" "
            "peptide_next_aa=\"%c\"/> \n",
            protein_id.c_str(),
            protein_annotation.c_str(),
            num_tol_term, 
            flanking_aas_prev,
            flanking_aas_next);
  }

  // print modifications
  print_modifications_xml(modified_peptide_sequence, peptide_sequence, file_);

  // print scores
  printScores(scores, scores_computed, ranks);

  // print post-search (analysis) fields
  printAnalysis(scores, scores_computed);

  // close the search_hit tag
  fprintf(file_, "    </search_hit>\n");
}

void PepXMLWriter::printScores(
  double* scores, 
  bool* scores_computed,
  int* ranks
) {
  string ranks_to_string[2] = {"sprank", "xcorr_rank"};
  for (int i = 0; i < NUMBER_SCORER_TYPES; i++) {
    if (i == BY_IONS_MATCHED || i == BY_IONS_TOTAL) {
      continue;
    }
    if (scores_computed[i]) {
      fprintf(file_, "        <search_score name=\"%s\" value=\"%.*f\" />\n",
        scorer_type_to_string((SCORER_TYPE_T)i), precision_, scores[i]);
      if (i < 2) {
        fprintf(file_, "        <search_score name=\"%s\" value=\"%i\" />\n",
          ranks_to_string[i].c_str(), ranks[i]);
      }
    }
  }
}


void PepXMLWriter::printAnalysis(double* scores,
                                 bool* scores_computed) {
  // only write the analysis_result section for post-search psms
  bool post_search = false;
  if( scores_computed[PERCOLATOR_SCORE] 
      || scores_computed[QRANKER_SCORE]
      || scores_computed[DECOY_XCORR_QVALUE] 
      || scores_computed[LOGP_QVALUE_WEIBULL_XCORR] ) {
    post_search = true;
  }

  if( !post_search ) {
    return;
  }

  SCORER_TYPE_T score = INVALID_SCORER_TYPE;
  SCORER_TYPE_T qval = INVALID_SCORER_TYPE;
  SCORER_TYPE_T pep = INVALID_SCORER_TYPE;
  string score_name;
  
  if (scores_computed[PERCOLATOR_SCORE]) {
    score = PERCOLATOR_SCORE;
    qval = PERCOLATOR_QVALUE;
    pep = PERCOLATOR_PEP;
    score_name = "percolator_score";
  } else if (scores_computed[QRANKER_SCORE]) {
    score =  QRANKER_SCORE; 
    qval = QRANKER_QVALUE;
    pep = QRANKER_PEP;
    score_name = "q-ranker_score";
  } else if (scores_computed[LOGP_QVALUE_WEIBULL_XCORR]) {
    qval = LOGP_QVALUE_WEIBULL_XCORR;
    pep = LOGP_WEIBULL_PEP;
  }

  fprintf(file_, "<analysis_result analysis=\"peptideprophet\">\n");
  fprintf(file_, "<peptideprophet_result probability=\"%.*f\">\n", 
          precision_, (1 - scores[pep]));
  fprintf(file_, "<search_score_summary>\n");
  if (!score_name.empty()) {
    fprintf(file_, "<parameter name=\"%s\" value=\"%.*f\"/>\n",
            score_name.c_str(), precision_, scores[score]);
  }
  fprintf(file_, "<parameter name=\"%s\" value=\"%.*f\"/>\n", 
          "qvalue", precision_, scores[qval]);
  fprintf(file_, "<parameter name=\"%s\" value=\"%.*f\"/>\n", 
          "pep", precision_, scores[pep]);
  fprintf(file_, "</search_score_summary>\n");
  fprintf(file_, "</peptideprophet_result>\n");
  fprintf(file_, "</analysis_result>\n");
}

/**
 * \brief prints both variable and static modifications for 
 *  peptide sequence in xml format to the specificed output file
 */
void PepXMLWriter::print_modifications_xml(
  const char* mod_seq,
  const char* pep_seq,
  FILE* output_file
) {
  carp(CARP_DEBUG, "print_modifications_xml:%s %s", mod_seq, pep_seq);
  // variable modifications
  int mod_precision = Params::GetInt("mod-precision");
  map<int, double> var_mods = find_variable_modifications(mod_seq);
  if (!var_mods.empty()) {
    fprintf(output_file, "<modification_info modified_peptide=\"%s\">\n", mod_seq);
    for (map<int, double>::iterator it = var_mods.begin(); it != var_mods.end(); ++it) {
      fprintf(output_file, "<mod_aminoacid_mass position=\"%i\" mass=\"%.*f\"/>\n",
              it->first,   //index
              mod_precision, it->second); //mass
    }
    fprintf(output_file, "</modification_info>\n");
  }

  // static modifications
  map<int, double> static_mods = find_static_modifications(pep_seq);
  if (!static_mods.empty()) {
    carp(CARP_DEBUG, "<modification_info modified_peptide=\"%s\">", pep_seq);
    fprintf(output_file, "<modification_info modified_peptide=\"%s\">\n", pep_seq);
    for (map<int, double>::iterator it = static_mods.begin(); it != static_mods.end(); ++it) {
      fprintf(output_file, "<mod_aminoacid_mass position=\"%i\" mass=\"%.*f\"/>\n",
              it->first,   //index
              mod_precision, it->second); //mass
    }
    fprintf(output_file, "</modification_info>\n");
  }
}

/**
 * \brief takes an empty mapping of index to mass
 * and extract information from mod sequence fill
 * up map
 */
map<int, double> PepXMLWriter::find_variable_modifications(const char* mod_seq) {
  bool monoisotopic = get_mass_type_parameter("isotopic-mass") != AVERAGE;
  map<int, double> mods;
  int seq_index = 1;
  const char* amino = mod_seq;
  const char* end = NULL;
  const char* start = NULL;
  // Parse returned string to find modifications within brackets
  while (*(amino) != '\0' && *(amino+1) != '\0') {
    if (*(amino+1) =='[') {
      start = amino+2;
      end = amino+2;
      while (*end != ']') {
        end++;
      }
      char* mass  = (char*)mymalloc(sizeof(char)*(end-start+1));
      strncpy(mass, start, end-start);
      mass[end-start] = '\0';
      mods[seq_index] = atof(mass) + AminoAcidUtil::GetMass(*amino, monoisotopic);
      amino = end;
      free(mass);
    } else if (*(amino+1) < 'A' || *(amino+1) > 'Z') { // a mod symbol
      double mass = 0; // sum up all adjacent symbols
      end = amino + 1;
      while( *end < 'A' || *end > 'Z' ) {
        mass += get_mod_mass_from_symbol(*end);
        end++;
      }
      mods[seq_index] = mass + AminoAcidUtil::GetMass(*amino, monoisotopic);
    }
    seq_index++;
    amino++;
  }
  return mods;
}

/**
 * \brief takes an empty mapping of index to mass
 * of static mods and a full mapping of var mods
 * to fill up the mapping of static mods
 */
map<int, double> PepXMLWriter::find_static_modifications(
  const char* peptide_sequence
) {
  map<int, double> static_mods;
  const char* seq_iter = peptide_sequence;
  MASS_TYPE_T isotopic_type = get_mass_type_parameter("isotopic-mass");
  for (size_t i = 0; peptide_sequence[i] != '\0'; i++) {
    char aa = peptide_sequence[i];
    double sum = AminoAcidUtil::GetMass(aa, isotopic_type == MONO);
    int modCount = 0;
    vector<const ModificationDefinition*> staticMods = ModificationDefinition::StaticMods(aa);
    for (vector<const ModificationDefinition*>::const_iterator j = staticMods.begin();
         j != staticMods.end();
         j++) {
      if ((*j)->Position() == ANY ||
          ((*j)->Position() == PEPTIDE_N && i == 0) ||
          ((*j)->Position() == PEPTIDE_C && peptide_sequence[i + 1] == '\0')) {
        ++modCount;
        sum += (*j)->DeltaMass();
      }
    }
    if (modCount > 0) {
      static_mods[i + 1] = sum;
    }
  }
  return static_mods;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
