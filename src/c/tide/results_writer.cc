// Tahmina Baker
//
// This file contains implementations for the classes that inherit from
// ResultsWriter (e.g. TextResultsWriter). Please see results_writer.h 
// for details.

#include <sstream>
#include "results_writer.h"


///////////////////////////////////////////////////////////////////////////////
// ResultsWriterUtils methods.
///////////////////////////////////////////////////////////////////////////////
char ResultsWriterUtils::GetAABefore(int first_aa_pos, 
                                     const pb::Protein& protein) {
  if (first_aa_pos > 0)
    return *(protein.residues().data() + first_aa_pos - 1); 
  return '-';
}

char ResultsWriterUtils::GetAAAfter(int last_aa_pos, 
                                    const pb::Protein& protein) {
  int aa_after_index = last_aa_pos + 1;
  if (aa_after_index < protein.residues().length())
    return *(protein.residues().data() + aa_after_index);
  return '-';
}

string ResultsWriterUtils::GetSeq(const pb::Peptide& pb_peptide, 
                                  const ProteinVec& proteins) {
  Peptide peptide(pb_peptide, proteins);
  if (FLAGS_show_mods)
    return peptide.SeqWithMods();
  return peptide.Seq();
}

int ResultsWriterUtils::GetNumProteins(const pb::Peptide& pb_peptide, 
                                       const AuxLocVec& aux_locs) {
  if (FLAGS_show_all_proteins) {
    if (pb_peptide.has_aux_locations_index()) {
      int aux_loc_idx = pb_peptide.aux_locations_index();
      const pb::AuxLocation* aux_loc = aux_locs[aux_loc_idx];
      return aux_loc->location_size();
    }
  }

  return 1;
}


///////////////////////////////////////////////////////////////////////////////
// TextResultsWriter methods.
///////////////////////////////////////////////////////////////////////////////
TextResultsWriter::TextResultsWriter(const ProteinVec& proteins, 
                                     const AuxLocVec& aux_locs,
                                     ostream& out_stream)
  : proteins_(proteins), aux_locs_(aux_locs), out_stream_(out_stream) {
  // Get the spectrum fields
  ResultsFieldMap spectrum_field_map;
  spectrum_field_map.insert(make_pair("spectrum_num", &spectrum_num_writer_));
  spectrum_field_map.insert(make_pair("mz", &mz_writer_));
  spectrum_field_map.insert(make_pair("charge", &charge_writer_));
  spectrum_field_map.insert(make_pair("rtime", &rtime_writer_));
  ParseFieldNames(FLAGS_spectrum_fields, ',', spectrum_field_map, 
                  spectrum_fields_);

  // Get the match fields
  ResultsFieldMap match_field_map;
  match_field_map.insert(make_pair("xcorr", &xcorr_writer_));
  match_field_map.insert(make_pair("sequence", &sequence_writer_));
  ParseFieldNames(FLAGS_match_fields, ',', match_field_map, 
                  match_fields_);

  // Get the protein fields
  ResultsFieldMap protein_field_map;
  protein_field_map.insert(make_pair("protein_name", &protein_name_writer_));
  protein_field_map.insert(make_pair("pos", &position_writer_));
  protein_field_map.insert(make_pair("aa_before", &aa_before_writer_));
  protein_field_map.insert(make_pair("aa_after", &aa_after_writer_));
  ParseFieldNames(FLAGS_protein_fields, ',', protein_field_map,
                  protein_fields_);
}

void TextResultsWriter::ParseFieldNames(const string& fields_str, char delim, 
                                        ResultsFieldMap& results_field_map, 
                                        vector<TextFieldWriter*>& fields_vec) {
  int cut_at = 0;
  string str = fields_str;
  while(!str.empty() && (cut_at != str.npos)) {
    string str_curr_field;
    if ((cut_at = str.find_first_of(delim)) == str.npos) {
      str_curr_field = str; // We've reached the last field
    } else {
      str_curr_field = str.substr(0, cut_at);
      str = str.substr(cut_at+1);
    }

    if (!str_curr_field.empty()) {
      TextFieldWriter *writer = results_field_map[str_curr_field];
      CHECK(NULL != writer) << str_curr_field << " is not a valid field.\n";
      fields_vec.push_back(writer);
    }
  }
}

void TextResultsWriter::WriteResults(pb::Results& pb_result) {
  WriteSpectrumFields(pb_result.spectrum(), pb_result.charge());
  for (int match = 0; match < pb_result.matches_size(); match++) {
    WriteMatchFields(pb_result.matches(match));

    const pb::Peptide& pb_peptide = pb_result.matches(match).peptide();
    int peptide_length = pb_peptide.length();

    // Write the first location by default
    const pb::Location& first_location = pb_peptide.first_location();
    const pb::Protein& first_protein = 
      *(proteins_[first_location.protein_id()]);
    WriteProteinFields(first_protein, first_location, peptide_length);

    // Write the auxiliary locations only if the user has asked for it
    if (FLAGS_show_all_proteins) {
      if (pb_peptide.has_aux_locations_index()) {
        int aux_loc_idx = pb_peptide.aux_locations_index();
        const pb::AuxLocation* aux_loc = aux_locs_[aux_loc_idx];
        for (int loc_idx = 0; loc_idx < aux_loc->location_size(); ++loc_idx) {
          const pb::Location& aux_location = aux_loc->location(loc_idx);
          const pb::Protein& aux_protein = 
            *(proteins_[aux_location.protein_id()]);
          WriteProteinFields(aux_protein, aux_location, peptide_length);
        }
      }
    }
  }
}

void TextResultsWriter::WriteSpectrumFields(const pb::Spectrum& spectrum, 
                                            int charge) {
  for (int i = 0; i < spectrum_fields_.size(); i++) {
    if (i > 0)
      out_stream_ << "\t";
    spectrum_fields_[i]->WriteSpectrumField(spectrum, charge, out_stream_);    
  }
  out_stream_ << endl;
}

void TextResultsWriter::WriteMatchFields(const pb::Match& match) {
  out_stream_ << "\t"; 
  for (int i = 0; i < match_fields_.size(); i++) {
    if (i > 0)
      out_stream_ << "\t";
    match_fields_[i]->WriteMatchField(proteins_, match, out_stream_);
  }
  out_stream_ << endl;
}

void TextResultsWriter::WriteProteinFields(const pb::Protein& protein, 
                                           const pb::Location& location, 
                                           int peptide_len) {
  out_stream_ << "\t\t"; 
  for (int i = 0; i < protein_fields_.size(); i++) {
    if (i > 0)
      out_stream_ << "\t";
    protein_fields_[i]->WriteProteinField(protein, location, peptide_len, 
                                          out_stream_);    
  }
  out_stream_ << endl;
}


///////////////////////////////////////////////////////////////////////////////
// SqtResultsWriter methods.
///////////////////////////////////////////////////////////////////////////////

SqtResultsWriter::SqtResultsWriter(const ProteinVec& proteins, 
                                   const AuxLocVec& aux_locs,
                                   SpectrumCollection& spectra, 
                                   const pb::ModTable* mod_table,
                                   ostream& out_stream,
                                   pb::Header& results_header, 
                                   const string& command_line) 
  : proteins_(proteins), aux_locs_(aux_locs), spectra_(spectra),
  out_stream_(out_stream) {
  WriteHeader(mod_table, results_header, command_line);
}

void SqtResultsWriter::WriteHeader(const pb::ModTable* mod_table, 
                                   pb::Header& results_header,
                                   const string& command_line) {
  out_stream_ << "H\tSQTGenerator\tTide" << endl;

  if (mod_table) {
    for (int i = 0; i < mod_table->static_mod_size(); ++i) {
      char aa = mod_table->static_mod(i).amino_acids()[0];
      double delta = mod_table->static_mod(i).delta();
      out_stream_ << "H\tStaticMod\t" << aa << "=" << delta << endl;
    }

    for (int i = 0; i < mod_table->variable_mod_size(); ++i) {
      char aa = mod_table->variable_mod(i).amino_acids()[0];
      double delta = mod_table->variable_mod(i).delta();
      out_stream_ << "H\tDiffMod\t" << aa << "=" << delta << endl;
    }
  }

  // Print the index command line, which can be found in the header of the
  // proteins protocol buffer file (.protix file)
  for (int i = 0; i < results_header.source_size(); i++) {
    pb::Header_Source& source = *(results_header.mutable_source(i));
    pb::Header* header = source.mutable_header();
    if (header->file_type() == pb::Header::RAW_PROTEINS) {
      out_stream_ << "H\tCommandLineIndex\t" << header->command_line() << endl;
      break;
    }
  }

  out_stream_ << "H\tCommandLineSearch\t" 
              << results_header.command_line() 
              << endl;

  out_stream_ << "H\tCommandLineResults\t" << command_line << endl;

  pb::Header_ResultsHeader& results_specific_header = 
    *(results_header.mutable_results_header());
  out_stream_ << "H\tMassWindow\t" << results_specific_header.mass_window() << endl;
  out_stream_ << "H\tTopMatches\t" << results_specific_header.top_matches() << endl;

  out_stream_ << "H\tShowAllProteins\t" << FLAGS_show_all_proteins << endl;
  out_stream_ << "H\tShowMods\t" << FLAGS_show_mods << endl;

  pb::Header_PeptidesHeader& peptides_specific_header = 
    *(results_specific_header.mutable_peptides_header());
  out_stream_ << "H\tMinMass\t" << peptides_specific_header.min_mass() << endl;
  out_stream_ << "H\tMaxMass\t" << peptides_specific_header.max_mass() << endl;
  out_stream_ << "H\tMinLength\t" << peptides_specific_header.min_length() << endl;
  out_stream_ << "H\tMaxLength\t" << peptides_specific_header.max_length() << endl;
  out_stream_ << "H\tEnzyme\t" << peptides_specific_header.enzyme() << endl;
  out_stream_ << "H\tFullDigestion\t" << peptides_specific_header.full_digestion() << endl;
  out_stream_ << "H\tMaxMissedCleavages\t" << peptides_specific_header.max_missed_cleavages() << endl;
  out_stream_ << "H\tMonoisotopicPrecursor\t" << peptides_specific_header.monoisotopic_precursor() << endl;
}

void SqtResultsWriter::WriteResults(pb::Results& pb_result) {
  // Call this BEFORE WriteSpectrumLine and WriteMatchLine because
  // they both rely on ranked sp scores to be available
  GenerateSpScoreInfo(pb_result);

  WriteSpectrumLine(pb_result);

  for (int match = 0; match < pb_result.matches_size(); match++) {
    WriteMatchLine(pb_result, match);
  
    const pb::Peptide& pb_peptide = pb_result.matches(match).peptide();
    const pb::Location& first_location = pb_peptide.first_location();
    const pb::Protein& first_protein = 
      *(proteins_[first_location.protein_id()]);
    WriteLocusLine(first_protein);
    if (FLAGS_show_all_proteins) {
      if (pb_peptide.has_aux_locations_index()) {
        int aux_loc_idx = pb_peptide.aux_locations_index();
        const pb::AuxLocation* aux_loc = aux_locs_[aux_loc_idx];
        for (int loc_idx = 0; loc_idx < aux_loc->location_size(); ++loc_idx) {
          const pb::Location& aux_location = aux_loc->location(loc_idx);
          const pb::Protein& aux_protein = 
            *(proteins_[aux_location.protein_id()]);
          WriteLocusLine(aux_protein);
        }
      }
    }
  }
}

void SqtResultsWriter::GenerateSpScoreInfo(pb::Results& pb_result) {
  sp_scores_.clear();
  smallest_sp_score_ = 0.0;
  total_ion_intensity_ = 0.0;

  Spectrum *spec = (*spectra_.Spectra())[pb_result.spectrum_index()];
  SpScorer sp_scorer = SpScorer(proteins_, *spec, pb_result.charge(),
                                spectra_.FindHighestMZ());

  for (int match = 0; match < pb_result.matches_size(); match++) {
    const pb::Peptide& pb_peptide = pb_result.matches(match).peptide();
    SpScorer::SpScoreData sp_score_data;
    sp_scorer.Score(pb_peptide, sp_score_data);
    sp_scores_.push_back(sp_score_data);
  }

  // Assign rankings to the sp scores stored in sp_score_data_
  // Also get the smallest sp score
  sp_scorer.RankSpScores(sp_scores_, &smallest_sp_score_);

  // Get the total ion intensity for the sp processed spectrum
  total_ion_intensity_ = sp_scorer.TotalIonIntensity();
}

void SqtResultsWriter::WriteSpectrumLine(pb::Results& pb_result) {
  const pb::Spectrum& spectrum = pb_result.spectrum();
  
  double singly_charged_mass = (spectrum.precursor_m_z() * 
                                pb_result.charge()) - 
                               (pb_result.charge() - 1)
                               *MassConstants::proton;

  int candidates = pb_result.matches_size(); // default if no stats
  if (pb_result.has_stats() && pb_result.stats().has_count())
    candidates = pb_result.stats().count();

  out_stream_ << "S\t" << spectrum.spectrum_number() << "\t" 
              << spectrum.spectrum_number() << "\t" 
              << pb_result.charge() << "\t" 
              << "0\t0" << "\t"  // Elapsed time and hostname - don't have it
              << singly_charged_mass << "\t" 
              << total_ion_intensity_ << "\t" // Total ion intensity
              << smallest_sp_score_ << "\t"
              << candidates << endl;
}

void SqtResultsWriter::WriteMatchLine(pb::Results& pb_result, int match) {
  const pb::Match& pb_match = pb_result.matches(match);
  double delta_cn = ResultsWriterUtils::GetDeltaCn(pb_match.xcorr(),
                    pb_result.matches(0).xcorr());
  const pb::Peptide& pb_peptide = pb_match.peptide();
  const pb::Location& first_location = pb_peptide.first_location();
  const pb::Protein& first_protein = *(proteins_[first_location.protein_id()]);
  
  string sequence = ResultsWriterUtils::GetSeq(pb_peptide, proteins_);
  
  int first_aa_pos = first_location.pos();
  char aa_before = ResultsWriterUtils::GetAABefore(first_aa_pos, 
                                                   first_protein);

  int last_aa_pos = first_location.pos() + pb_peptide.length() - 1;
  char aa_after = ResultsWriterUtils::GetAAAfter(last_aa_pos,
                                                 first_protein);

  out_stream_ << "M\t" 
              << match + 1 << "\t" // rank by Xcorr
              << sp_scores_[match].sp_rank << "\t"
              << pb_peptide.mass() + MassConstants::proton << "\t"
              << delta_cn << "\t" // 1.0 - (xcorr_current/xcorr_top)
              << pb_match.xcorr() << "\t"
              << sp_scores_[match].sp_score << "\t"
              << sp_scores_[match].matched_ions << "\t"
              << sp_scores_[match].total_ions << "\t"
              << aa_before << "." << sequence << "." << aa_after << "\t"
              << "U" // Validation status -> U = unknown
              << endl;
}

void SqtResultsWriter::WriteLocusLine(const pb::Protein& protein) {
  out_stream_ << "L\t" 
              << protein.name()
              << endl;
}

///////////////////////////////////////////////////////////////////////////////
// PepXMLResultsWriter methods.
///////////////////////////////////////////////////////////////////////////////

PepXMLResultsWriter::PepXMLResultsWriter(const ProteinVec& proteins, 
                                   const AuxLocVec& aux_locs,
                                   ostream& out_stream) 
  : proteins_(proteins), aux_locs_(aux_locs), base_writer_(out_stream, 0),
  index_(0) {
  out_stream << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  out_stream << "<?xml-stylesheet type=\"text/xsl\"href=\""
             << "http://regis-web.systemsbiology.net/pepXML_std.xsl\"?>"
             << endl;
  base_writer_.StartTag("msms_pipeline_analysis");
  base_writer_.FinishStartTag();
}

void PepXMLResultsWriter::WriteSpectrum(pb::Results& pb_result, XMLWriter& xml) {
  xml.StartTag("spectrum_query");
  const pb::Spectrum& spectrum = pb_result.spectrum();
  stringstream str_stream;
  str_stream << "000." << spectrum.spectrum_number() << "." 
             << spectrum.spectrum_number() << "." << pb_result.charge();
  xml.KeyVal("spectrum", str_stream.str());
  xml.KeyVal("start_scan", spectrum.spectrum_number());
  xml.KeyVal("end_scan", spectrum.spectrum_number());
  neutral_mass_ = ((spectrum.precursor_m_z() - MassConstants::proton)
		   * pb_result.charge());
  xml.KeyVal("precursor_neutral_mass", neutral_mass_);
  xml.KeyVal("assumed_charge", pb_result.charge());
  xml.KeyVal("index", ++index_);
  xml.KeyVal("retention_time_sec", "");
  xml.FinishStartTag();
}

void PepXMLResultsWriter::WriteMatch(const pb::Match& pb_match, int match,
				     double delta_cn, XMLWriter& xml) {
  const pb::Peptide& pb_peptide = pb_match.peptide();
  const pb::Location& first_location = pb_peptide.first_location();
  const pb::Protein& first_protein = *(proteins_[first_location.protein_id()]);

  xml.StartTag("search_hit");
  xml.KeyVal("hit_rank", match+1);
  string peptide = ResultsWriterUtils::GetSeq(pb_peptide, proteins_);
  xml.KeyVal("peptide", peptide);

  int first_aa_pos = first_location.pos();
  char aa_before = ResultsWriterUtils::GetAABefore(first_aa_pos, 
                                                   first_protein);
  xml.KeyVal("peptide_prev_aa", aa_before);
  int last_aa_pos = first_location.pos() + pb_peptide.length() - 1;
  char aa_after = ResultsWriterUtils::GetAAAfter(last_aa_pos,
                                                 first_protein);
  xml.KeyVal("peptide_next_aa", aa_after);
  xml.KeyVal("protein", first_protein.name());
  xml.KeyVal("num_tot_proteins", ResultsWriterUtils::GetNumProteins(pb_peptide,
             aux_locs_));
  xml.KeyVal("num_matched_ions", ""); // DON'T HAVE
  xml.KeyVal("tot_num_ions", ""); // DON'T HAVE
  xml.KeyVal("calc_neutral_pep_mass", pb_peptide.mass());
  xml.KeyVal("massdiff", pb_peptide.mass() - neutral_mass_);
  xml.KeyVal("num_tol_term", ""); // DON'T HAVE
  xml.KeyVal("num_missed_cleavages", ""); // DON'T HAVE
  xml.KeyVal("is_rejected", 0);

  if (pb_peptide.modifications_size() > 0) {
    XMLWriter mod_info(xml);
    mod_info.StartTag("modification_info");
    //TODO
  }
  
  XMLWriter score(xml);
  score.StartTag("search_score");
  score.KeyVal("xcorr", pb_match.xcorr());
  score.CloseTag();
  score.StartTag("search_score");
  score.KeyVal("deltacn", delta_cn);
  score.CloseTag();

  // Set "deltacnstar=0" to indicate that we are not checking to see if the 
  // peptide corresponding to the current xcorr is homologous to the peptide 
  // corresponding to the best xcorr.
  score.StartTag("search_score");
  score.KeyVal("deltacnstar", 0);
  score.CloseTag();

  // DON'T HAVE
  score.StartTag("search_score");
  score.KeyVal("spscore", "");
  score.CloseTag();

  // DON'T HAVE
  score.StartTag("search_score");
  score.KeyVal("sprank", "");
  score.CloseTag();
}

void PepXMLResultsWriter::WriteResults(pb::Results& pb_result) {
  if (pb_result.matches_size() == 0)
    return;
  XMLWriter spectrum_query(base_writer_);
  WriteSpectrum(pb_result, spectrum_query);
  for (int match = 0; match < pb_result.matches_size(); match++) {
    const pb::Match& pb_match = pb_result.matches(match);
    XMLWriter search_result(spectrum_query);
    search_result.StartTag("search_result");
    
    // delta_cn = 1.0 - (xcorr_current/xcorr_top)
    double delta_cn = ResultsWriterUtils::GetDeltaCn(pb_match.xcorr(),
                      pb_result.matches(0).xcorr());

    XMLWriter sub(search_result);
    WriteMatch(pb_match, match, delta_cn, sub);
  }
}
