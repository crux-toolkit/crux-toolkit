#include "SQTWriter.h"

using namespace Crux;

SQTWriter::SQTWriter() {

}

SQTWriter::~SQTWriter() {
  closeFile();
}

void SQTWriter::openFile(string filename) {
  file_ = create_file(filename.c_str(), get_boolean_parameter("overwrite"));
  if (!file_) {
    carp(CARP_FATAL, "Error creating file '%s'.", filename.c_str());
  }
}

void SQTWriter::closeFile() {
  if (file_ && file_->is_open()) {
    file_->close();
    delete file_;
    file_ = NULL;
  }
}

void SQTWriter::writeHeader(
  string database,
  int num_proteins,
  bool is_decoy
) {
  if (!file_->is_open()) {
    return;
  }

  streamsize old_precision = file_->precision();

  time_t hold_time = time(0);

  if (file_exists(database.c_str()) && is_directory(database.c_str())) {
    database = Index::getBinaryFastaName(database.c_str());
  }

  MASS_TYPE_T mass_type = get_mass_type_parameter("isotopic-mass");
  char precursor_masses[64];
  mass_type_to_string(mass_type, precursor_masses);

  mass_type = get_mass_type_parameter("fragment-mass");
  char fragment_masses[64];
  mass_type_to_string(mass_type, fragment_masses);

  double tol = get_double_parameter("precursor-window");
  double frag_mass_tol = get_double_parameter("mz-bin-width") / 2.0;

  SCORER_TYPE_T score = get_scorer_type_parameter("prelim-score-type");
  string prelim_score_type(scorer_type_to_string(score));
  score = get_scorer_type_parameter("score-type");
  string score_type(scorer_type_to_string(score));

  *file_ << "H\tSQTGenerator Crux" << endl
         << "H\tSQTGeneratorVersion 1.0" << endl
         << "H\tComment Crux was written by..." << endl
         << "H\tComment ref..." << endl
         << "H\tStartTime\t" << string(ctime(&hold_time))
         << "H\tEndTime                               " << endl
         << "H\tDatabase\t" << database << endl;

  if (is_decoy) {
    *file_ << "H\tComment\tDatabase shuffled; these are decoy matches" << endl;
  }

  *file_ << "H\tDBSeqLength\t?" << endl
         << "H\tDBLocusCount\t" << num_proteins << endl
         << "H\tPrecursorMasses\t" << precursor_masses << endl
         << "H\tFragmentMasses\t" << fragment_masses << endl //?????????
         << fixed
         << "H\tAlg-PreMasTol\t" << setprecision(1) << tol << endl
         << "H\tAlg-FragMassTol\t" << setprecision(2) << frag_mass_tol << endl
         << "H\tAlg-XCorrMode\t0" << endl;
  file_->unsetf(ios_base::fixed);
  file_->precision(old_precision);

  *file_ << "H\tComment\tpreliminary algorithm " << prelim_score_type << endl
         << "H\tComment\tfinal algorithm " << score_type << endl;

  int aa = 0;
  char aa_str[2];
  aa_str[1] = '\0';
  int alphabet_size = (int)'A' + ((int)'Z'-(int)'A');
  MASS_TYPE_T isotopic_type = get_mass_type_parameter("isotopic-mass");

  *file_ << fixed << setprecision(3);
  for(aa = (int)'A'; aa < alphabet_size -1; aa++){
    aa_str[0] = (char)aa;
    double mod = get_double_parameter(aa_str);
    if( mod != 0 ){
      //      double mass = mod + get_mass_amino_acid(aa, isotopic_type);
      double mass = get_mass_amino_acid(aa, isotopic_type);
      *file_ << "H\tStaticMod\t" << aa_str << "=" << mass << endl;
    }
  }
  file_->unsetf(ios_base::fixed);
  file_->precision(old_precision);

  // print dynamic mods, if any
  // format DiffMod <AAs><symbol>=<mass change>
  AA_MOD_T** aa_mod_list = NULL;
  int num_mods = get_all_aa_mod_list(&aa_mod_list);
  int mod_idx = 0;
  *file_ << fixed << setprecision(2);
  for(mod_idx = 0; mod_idx < num_mods; mod_idx++){
    
    AA_MOD_T* aamod = aa_mod_list[mod_idx];
    char* aa_list_str = aa_mod_get_aa_list_string(aamod);
    char aa_symbol = aa_mod_get_symbol(aamod);
    double mass_dif = aa_mod_get_mass_change(aamod);

    *file_ << "H\tDiffMod\t" << aa_list_str << aa_symbol << "="
           << (mass_dif >= 0 ? "+" : "-") << mass_dif << endl;
    free(aa_list_str);
  }
  file_->unsetf(ios_base::fixed);
  file_->precision(old_precision);

  num_mods = get_c_mod_list(&aa_mod_list);
  for(mod_idx = 0; mod_idx < num_mods; mod_idx++){
    AA_MOD_T* aamod = aa_mod_list[mod_idx];
    char aa_symbol = aa_mod_get_symbol(aamod);

    *file_ << "H\tComment\tMod " << aa_symbol
           << " is a C-terminal modification" << endl;
  }

  num_mods = get_n_mod_list(&aa_mod_list);
  for(mod_idx = 0; mod_idx < num_mods; mod_idx++){
    AA_MOD_T* aamod = aa_mod_list[mod_idx];
    char aa_symbol = aa_mod_get_symbol(aamod);

    *file_ << "H\tComment\tMod " << aa_symbol
           << " is a N-terminal modification" << endl;
  }

  //for letters in alphabet
  //  double mod = get_double_parameter(letter);
  //  if mod != 0
  //     double mass = mod + getmass(letter);
  //     fprintf(output, "H\tStaticMod\t%s=%.3f\n", letter, mass);
  //  fprintf(output, "H\tStaticMod\tC=160.139\n");
  *file_ << "H\tAlg-DisplayTop\t" << get_int_parameter("top-match") << endl;
          //          get_int_parameter("max-sqt-result")); 
  // this is not correct for an sqt from analzyed matches

  ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
  DIGEST_T digestion = get_digest_type_parameter("digestion");
  char* enz_str = enzyme_type_to_string(enzyme);
  char* dig_str = digest_type_to_string(digestion);
  string custom_str;
  if( enzyme == CUSTOM_ENZYME){
    string rule = get_string_parameter_pointer("custom-enzyme");
    custom_str = ", custom pattern: " + rule;
  }
  *file_ << "H\tEnzymeSpec\t" << enz_str << "-" << dig_str << custom_str << endl;
  free(enz_str);
  free(dig_str);

  // write a comment that says what the scores are
  *file_ << "H\tLine fields: S, scan number, scan number,"
         << "charge, 0, server, precursor mass, 0, 0, number of matches" << endl
         << "H\tLine fields: M, rank by xcorr score, rank by sp score, "
         << "peptide mass, deltaCn, xcorr score, sp score, number ions matched, "
         << "total ions compared, sequence" << endl;
}

void SQTWriter::writeSpectrum(
  Spectrum* spectrum,
  SpectrumZState& z_state,
  int num_matches
) {
  if (!file_->is_open()) {
    return;
  }
  
  *file_ << "S"
         << "\t" << spectrum->getFirstScan()
         << "\t" << spectrum->getLastScan()
         << "\t" << z_state.getCharge()
         << "\t0.0" // process time
         << "\tserver"
         << "\t" << setprecision(get_int_parameter("mass-precision"))
         << z_state.getSinglyChargedMass()
         << "\t0.0" // tic
         << "\t" << setprecision(get_int_parameter("precision"))
         << "0.0" // lowest sp
         << "\t" << num_matches
         << endl;
}

void SQTWriter::writePSM(
  Peptide* peptide,
  FLOAT_T xcorr_score,
  int xcorr_rank,
  FLOAT_T sp_score,
  int sp_rank,
  FLOAT_T delta_cn,
  int b_y_matched,
  int b_y_total,
  bool is_decoy
) {
  if (!file_->is_open()) {
    return;
  }

  int length = peptide->getLength();
  peptide->getModifiedAASequence();
  MODIFIED_AA_T* mod_seq =
    copy_mod_aa_seq(peptide->getModifiedAASequence(), length);
  if (!mod_seq) {
    return;
  }

  Protein* protein = peptide->getPeptideSrc()->getParentProtein();
  string seq_str;
  if (protein->isPostProcess()) {
    PostProcessProtein* post_process_protein = (PostProcessProtein*)protein;
    seq_str =
      post_process_protein->getNTermFlankingAA() + "." +
      string(post_process_protein->getSequencePointer()) +
      "." + post_process_protein->getCTermFlankingAA();
  } else {
    char* seq = peptide->getSequenceSqt();
    seq_str = seq;
    free(seq);
  }

  *file_ << "M"
         << "\t" << xcorr_rank
         << "\t" << sp_rank
         << "\t" << setprecision(get_int_parameter("mass-precision"))
         << peptide->getPeptideMass() + MASS_PROTON
         << "\t" << delta_cn
         << "\t" << setprecision(get_int_parameter("precision"))
         << xcorr_score
         << "\t" << sp_score
         << "\t" << b_y_matched
         << "\t" << b_y_total
         << "\t" << seq_str
         << endl;

  for (PeptideSrcIterator iter = peptide->getPeptideSrcBegin();
       iter != peptide->getPeptideSrcEnd();
       ++iter) {
    PeptideSrc* peptide_src = *iter;
    Protein* protein = peptide_src->getParentProtein();
    char* protein_id = protein->getId();
    string protein_id_str(protein_id);
    free(protein_id);
    if (is_decoy && protein->getDatabase()->getDecoyType() == NO_DECOYS) {
      protein_id_str = get_string_parameter_pointer("decoy-prefix") + protein_id_str;
    }
    *file_ << "L"
           << "\t" << protein_id_str
           << endl;
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

