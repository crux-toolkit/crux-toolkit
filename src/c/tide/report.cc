// Tahmina Baker

// This file contains implementations for the Reporter, TextReporter and 
// PBReporter classes. Please see report.h for details.
  
#include <iomanip>
#include "report.h"
#include "peptides.pb.h"

using namespace std;
using google::protobuf::int64;

///////////////////////////////////////////////////////////////////////////////
// TextReporter class methods
///////////////////////////////////////////////////////////////////////////////
void TextReporter::ReportSpectrum(const Spectrum *spectrum, int charge, 
                                  int spectrum_index, pb::Stats* stats) {
  spectrum_ = spectrum;
  charge_ = charge;
}

void TextReporter::ReportMatch(int score, const Peptide& peptide) {
  cout << spectrum_->SpectrumNumber() << "\t"
       << fixed << setprecision(2) << spectrum_->PrecursorMZ() << "\t"
       << charge_ << "\t"
       << setprecision(6) << (score/100000000.0) << "\t"
#ifdef DEBUG
       << peptide.PB()->id() << "\t"
#endif
       << peptide.SeqWithMods() << endl;

#if 0
    TheoreticalPeakSetMakeAll workspace(2000);
    peptide.ComputeTheoreticalPeaks(&workspace);
    TheoreticalPeakArr charge1(2000), charge2(2000);
    workspace.GetPeaks(&charge1, NULL, &charge2, NULL, NULL);

    cout << "CHARGE 1:\n";
    for (int i = 0; i < charge1.size(); ++i)
      cout << "TPeak[" << charge1[i].Bin() << "] = "
           << charge1[i].Type() << endl;

    cout << "CHARGE 2:\n";
    for (int i = 0; i < charge2.size(); ++i)
      cout << "TPeak[" << charge2[i].Bin() << "] = "
           << charge2[i].Type() << endl;    
#endif
}

// TextReporter writes as it reports matches, so nothing to do here.
void TextReporter::WriteReport() {
} 


///////////////////////////////////////////////////////////////////////////////
// PBReporter class methods
///////////////////////////////////////////////////////////////////////////////
PBReporter::PBReporter(const string& filename, pb::Header &results_header)
                    : results_writer_(AbsPath(filename), results_header) {
  CHECK(results_writer_.OK());
}

void PBReporter::ReportSpectrum(const Spectrum *spectrum, int charge, 
                                int spectrum_index, pb::Stats* stats) {
  spectrum_ = spectrum;
  pb_results_.Clear();
  pb_results_.set_charge(charge);
  pb_results_.set_spectrum_index(spectrum_index);
  pb_results_.mutable_stats()->CopyFrom(*stats);
  pb::Spectrum* results_spectrum = pb_results_.mutable_spectrum();
  results_spectrum->set_spectrum_number(spectrum->SpectrumNumber());
  results_spectrum->set_precursor_m_z(spectrum->PrecursorMZ());
  results_spectrum->set_rtime(spectrum->RTime());
  for (int i = 0; i < spectrum->NumChargeStates(); ++i)
    results_spectrum->add_charge_state(spectrum->ChargeState(i));
}

void PBReporter::ReportMatch(int score, const Peptide& peptide) {
  pb::Match* pb_match = pb_results_.add_matches();
  pb_match->set_xcorr(score/100000000.0);
  
  // We need rebuild a pb::peptide from the Peptide object passed in
  pb::Peptide* matching_peptide = pb_match->mutable_peptide();
  matching_peptide->set_id(peptide.Id());
  matching_peptide->set_mass(peptide.Mass());
  matching_peptide->set_length(peptide.Len());
  if (peptide.HasAuxLocationsIndex())
    matching_peptide->set_aux_locations_index(peptide.AuxLocationsIndex());
  
  // Copy over all the modifications for this Peptide
  const ModCoder::Mod* mods;
  int num_mods = peptide.Mods(&mods);
  for (int i = 0; i < num_mods; i++) {
    matching_peptide->add_modifications(mods[i]);
  }

  // Copy over the Peptide's first location within the first protein
  pb::Location* first_location = matching_peptide->mutable_first_location();
  first_location->set_protein_id(peptide.FirstLocProteinId());
  first_location->set_pos(peptide.FirstLocPos());  
}

// During WriteReport, PBReporter writes the reported spectrum and matches
// to the Results protocol buffer.
void PBReporter::WriteReport() {
  CHECK(results_writer_.Write(&pb_results_));
}
