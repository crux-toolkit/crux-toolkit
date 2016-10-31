#include "SQTWriter.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "util/GlobalParams.h"

using namespace Crux;

SQTWriter::SQTWriter() {
}

SQTWriter::~SQTWriter() {
  closeFile();
}

void SQTWriter::openFile(string filename) {
  file_ = FileUtils::GetWriteStream(filename, Params::GetBool("overwrite"));
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

  MASS_TYPE_T mass_type = GlobalParams::getIsotopicMass();
  char precursor_masses[64];
  mass_type_to_string(mass_type, precursor_masses);

  mass_type = get_mass_type_parameter("fragment-mass");
  char fragment_masses[64];
  mass_type_to_string(mass_type, fragment_masses);

  double tol = Params::GetDouble("precursor-window");
  double frag_mass_tol = Params::GetDouble("mz-bin-width") / 2.0;

  string prelim_score_type = Params::GetString("prelim-score-type");
  string score_type = Params::GetString("score-type");

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

  for (char aa = 'A'; aa <= 'Z'; aa++) {
    vector<const ModificationDefinition*> staticMods = ModificationDefinition::StaticMods(aa);
    for (vector<const ModificationDefinition*>::const_iterator i = staticMods.begin();
         i != staticMods.end();
         i++) {
      *file_ << "H\tStaticMod\t" << aa << '='
             << StringUtils::ToString((*i)->DeltaMass(), Params::GetInt("mod-precision"))
             << endl;
    }
  }

  // print dynamic mods, if any
  // format DiffMod <AAs><symbol>=<mass change>
  vector<const ModificationDefinition*> varMods = ModificationDefinition::VarMods();
  for (vector<const ModificationDefinition*>::const_iterator i = varMods.begin();
       i != varMods.end();
       i++) {
    char symbol = (*i)->Symbol();
    *file_ << "H\tDiffMod\t" << StringUtils::Join((*i)->AminoAcids()) << symbol << "=";
    if ((*i)->DeltaMass() >= 0) {
      *file_ << '+';
    }
    *file_ << StringUtils::ToString((*i)->DeltaMass(), Params::GetInt("mod-precision"))
           << endl;
    switch ((*i)->Position()) {
    case PEPTIDE_N:
    case PROTEIN_N:
      *file_ << "H\tComment\tMod " << symbol << " is an N-terminal modification" << endl;
      break;
    case PEPTIDE_C:
    case PROTEIN_C:
      *file_ << "H\tComment\tMod " << symbol << " is an C-terminal modification" << endl;
      break;
    }
  }

  *file_ << "H\tAlg-DisplayTop\t" << Params::GetInt("top-match") << endl;
  // this is not correct for an sqt from analzyed matches

  ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
  DIGEST_T digestion = get_digest_type_parameter("digestion");
  const char* enz_str = enzyme_type_to_string(enzyme);
  const char* dig_str = digest_type_to_string(digestion);
  string custom_str;
  if (enzyme == CUSTOM_ENZYME) {
    string rule = Params::GetString("custom-enzyme");
    custom_str = ", custom pattern: " + rule;
  }
  *file_ << "H\tEnzymeSpec\t" << enz_str << "-" << dig_str << custom_str << endl;


  *file_ << "H\tLine fields: S, scan number, scan number, "
         << "charge, 0, server, experimental mass, total ion intensity, "
         << "lowest Sp, number of matches" << endl
         << "H\tLine fields: M, rank by xcorr score, rank by sp score, "
         << "peptide mass, deltaCn, xcorr score, sp score, number ions matched, "
         << "total ions compared, sequence, validation status" << endl;
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
         << "\t"
         << StringUtils::ToString(
            z_state.getSinglyChargedMass(), Params::GetInt("mass-precision"))
         << "\t";

  // print total ion intensity if exists...
  if (spectrum->hasTotalEnergy()) {
    *file_ << StringUtils::ToString(spectrum->getTotalEnergy(), 2);
  }
  *file_ << "\t";

  // print lowest sp if exists
  if (spectrum->hasLowestSp()) {
    *file_ << StringUtils::ToString(spectrum->getLowestSp(), Params::GetInt("precision"));
  }
  *file_ << "\t";

  if (num_matches != 0) {
    *file_ << num_matches;
  }
  *file_ << endl;

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
  MODIFIED_AA_T* mod_seq = copy_mod_aa_seq(peptide->getModifiedAASequence(), length);
  if (!mod_seq) {
    return;
  }

  Protein* protein = peptide->getPeptideSrc()->getParentProtein();
  string seq_str = peptide->getSequenceSqt();

  *file_ << "M"
         << "\t" << xcorr_rank
         << "\t" << sp_rank
         << "\t" << StringUtils::ToString(
                    peptide->calcModifiedMass() + MASS_PROTON, Params::GetInt("mass-precision"))
         << "\t" << StringUtils::ToString(delta_cn, 2)
         << "\t" << StringUtils::ToString(xcorr_score, Params::GetInt("precision"))
         << "\t" << StringUtils::ToString(sp_score, Params::GetInt("precision"))
         << "\t" << b_y_matched
         << "\t" << b_y_total
         << "\t" << seq_str
         << "\tU"
         << endl;

  for (PeptideSrcIterator iter = peptide->getPeptideSrcBegin();
       iter != peptide->getPeptideSrcEnd();
       ++iter) {
    PeptideSrc* peptide_src = *iter;
    Protein* protein = peptide_src->getParentProtein();
    string protein_id_str = protein->getId();
    if (is_decoy && protein->getDatabase()->getDecoyType() == NO_DECOYS) {
      protein_id_str = Params::GetString("decoy-prefix") + protein_id_str;
    }
    *file_ << "L" << "\t" << protein_id_str << endl;
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

