// Tahmina Baker
//
// This file contains classes for writing search results that are read from a 
// results protocol buffer file.
//
// The ResultsWriter class is a generic base class for displaying the results.
// Any format specific results writer (e.g. text, pepxml, sqt) should inherit
// from this class.
//
// The TextResultsWriter class is a text format specific results writer. It
// uses the TextFieldWriter and all classes that inherit from it to display
// the results fields in the order specified by the user.
//
// The SqtResultsWriter class writes results out in the SQT format. Details 
// on the format can be found at:
// http://noble.gs.washington.edu/proj/crux/sqt-format.html
// NOTE: We currently do not have all the fields for the SQT format. We write
// out a 0 for these fields.
//
// The PepXmlResultsWriter class writes results out in the pep.xml format. 
// Details on the schema can be found at:
// http://sashimi.svn.sourceforge.net/viewvc/sashimi/trunk/trans_proteomic_pipeline/schema/
// NOTE: We leave blank fields we currently don't have for the pep.xml format.


#ifndef RESULTS_WRITER_H
#define RESULTS_WRITER_H

#include <iomanip>
#include "results.pb.h"
#include "records.h"
#include "crux_sp_spectrum.h"
#include "sp_scorer.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK((x))

DECLARE_string(match_fields);
DECLARE_string(spectrum_fields);
DECLARE_string(protein_fields);
DECLARE_bool(show_mods);
DECLARE_bool(show_all_proteins);


///////////////////////////////////////////////////////////////////////////////
// Common utility methods needed by results writers.
///////////////////////////////////////////////////////////////////////////////
namespace ResultsWriterUtils {
  char GetAABefore(int first_aa_pos, const pb::Protein& protein);
  char GetAAAfter(int last_aa_pos, const pb::Protein& protein);
  string GetSeq(const pb::Peptide& pb_peptide, const ProteinVec& proteins);
  int GetNumProteins(const pb::Peptide& pb_peptide, const AuxLocVec& aux_locs);
  inline double GetDeltaCn(double xcorr_n, double xcorr_top) {
    if (xcorr_top > 0)
      return 1.0 - (xcorr_n/xcorr_top);
    return 0.0;
  }
};

///////////////////////////////////////////////////////////////////////////////
// Base class to display results.
///////////////////////////////////////////////////////////////////////////////
class ResultsWriter {
public:
  virtual void WriteResults(pb::Results& pb_result) = 0;
  virtual ~ResultsWriter() {};
};


///////////////////////////////////////////////////////////////////////////////
// Classes for text specific results writer. The main class for writing text
// results is TextResultsWriter below. The rest of the classes are used for
// writing individual text fields in the order specified by the user.
///////////////////////////////////////////////////////////////////////////////
class TextFieldWriter {
  public:
    // Subclasses are to override exactly one of these
    virtual void WriteSpectrumField(const pb::Spectrum& spectrum, int charge, 
                                    ostream& out_stream_) {
    }
    virtual void WriteMatchField(const ProteinVec& proteins,
                                 const pb::Match& match, 
                                 ostream& out_stream_) {
    }
    virtual void WriteProteinField(const pb::Protein& protein, 
                                   const pb::Location& location, 
                                   int peptide_len, ostream& out_stream_) {
    }
};

class TextSpectrumNumWriter : public TextFieldWriter {
  public:
    void WriteSpectrumField(const pb::Spectrum& spectrum, int charge,
                            ostream& out_stream) {
      out_stream << spectrum.spectrum_number();
    }
};

class TextMzWriter : public TextFieldWriter {
  public:
    void WriteSpectrumField(const pb::Spectrum& spectrum, int charge,
                            ostream& out_stream) {
      out_stream << fixed << setprecision(2) 
                 << spectrum.precursor_m_z();
    }
};

class TextChargeWriter : public TextFieldWriter {
  public:
    void WriteSpectrumField(const pb::Spectrum& spectrum, int charge,
                            ostream& out_stream) {
      out_stream << charge;
    }
};

class TextRTimeWriter : public TextFieldWriter {
  public:
    void WriteSpectrumField(const pb::Spectrum& spectrum, int charge,
                            ostream& out_stream) {
      out_stream << spectrum.rtime();
    }
};

class TextXcorrWriter : public TextFieldWriter {
  public:
    void WriteMatchField(const ProteinVec& proteins, const pb::Match& match, 
                         ostream& out_stream) {
      out_stream << setprecision(6) << match.xcorr();
    }
};

class TextSequenceWriter : public TextFieldWriter {
  public:
    void WriteMatchField(const ProteinVec& proteins, const pb::Match& match, 
                         ostream& out_stream) {
      out_stream << ResultsWriterUtils::GetSeq(match.peptide(), proteins);
    }
};

class TextProteinNameWriter : public TextFieldWriter {
  public:
    void WriteProteinField(const pb::Protein& protein, 
                           const pb::Location& location, 
                           int peptide_len, ostream& out_stream) {
      out_stream << protein.name();
    }
};

class TextPositionWriter : public TextFieldWriter {
  public:
    void WriteProteinField(const pb::Protein& protein, 
                           const pb::Location& location, 
                           int peptide_len, ostream& out_stream) {
      out_stream << location.pos();
    }
};

class TextAABeforeWriter : public TextFieldWriter {
  public:
    void WriteProteinField(const pb::Protein& protein, 
                           const pb::Location& location, 
                           int peptide_len, ostream& out_stream) {
      int first_aa_pos = location.pos();
      out_stream << ResultsWriterUtils::GetAABefore(first_aa_pos, protein);
    }
};

class TextAAAfterWriter : public TextFieldWriter {
  public:
    void WriteProteinField(const pb::Protein& protein, 
                           const pb::Location& location, 
                           int peptide_len, ostream& out_stream) {
      int last_aa_pos = location.pos() + peptide_len - 1;
      out_stream << ResultsWriterUtils::GetAAAfter(last_aa_pos, protein);
    }
};

// The main class for writing text output
class TextResultsWriter : public ResultsWriter {
 public:
  TextResultsWriter(const ProteinVec& proteins, const AuxLocVec& aux_locs, 
                    ostream& out_stream);

  void WriteResults(pb::Results& pb_result);

 private:
  typedef map<string, TextFieldWriter*> ResultsFieldMap;
  
  void ParseFieldNames(const string& fields_str, char delim,
                       ResultsFieldMap& results_field_map, 
                       vector<TextFieldWriter*>& fields_vec);
  void WriteSpectrumFields(const pb::Spectrum& spectrum, int charge);
  void WriteMatchFields(const pb::Match& match);
  void WriteProteinFields(const pb::Protein& protein, 
                          const pb::Location& location, 
                          int peptide_len);

  vector<TextFieldWriter*> spectrum_fields_;
  vector<TextFieldWriter*> match_fields_;
  vector<TextFieldWriter*> protein_fields_;
  
  TextSpectrumNumWriter spectrum_num_writer_;
  TextMzWriter mz_writer_;
  TextChargeWriter charge_writer_;
  TextRTimeWriter rtime_writer_;
  TextXcorrWriter xcorr_writer_;
  TextSequenceWriter sequence_writer_;
  TextProteinNameWriter protein_name_writer_;
  TextPositionWriter position_writer_;
  TextAABeforeWriter aa_before_writer_;
  TextAAAfterWriter aa_after_writer_;

  const ProteinVec& proteins_;
  const AuxLocVec& aux_locs_;

  ostream& out_stream_;
};

///////////////////////////////////////////////////////////////////////////////
// Classes for sqt specific results writer.
///////////////////////////////////////////////////////////////////////////////

class SqtResultsWriter : public ResultsWriter {
 public:
  SqtResultsWriter(const ProteinVec& proteins, const AuxLocVec& aux_locs,
                   SpectrumCollection& spectra, const pb::ModTable* mod_table,
                   ostream& out_stream, pb::Header& results_header, 
                   const string& command_line);
  void WriteResults(pb::Results& pb_result);

 private:
  void WriteHeader(const pb::ModTable* mod_table, 
                   pb::Header& results_header,
                   const string& command_line);
  void WriteSpectrumLine(pb::Results& pb_result);
  void WriteMatchLine(pb::Results& pb_result, int match);
  void WriteLocusLine(const pb::Protein& protein);
  void GenerateSpScoreInfo(pb::Results& pb_result);
  
  const ProteinVec& proteins_;
  const AuxLocVec& aux_locs_;
  SpectrumCollection& spectra_;

  vector<SpScorer::SpScoreData> sp_scores_;
  double smallest_sp_score_;
  double total_ion_intensity_;

  ostream& out_stream_;
};

///////////////////////////////////////////////////////////////////////////////
// Classes for .pep.xml specific results writer.
///////////////////////////////////////////////////////////////////////////////
class XMLWriter {
 public:
  XMLWriter(ostream& out_stream, int indent)
    : out_stream_(out_stream), indent_(indent), closed_(true) {
  }

  XMLWriter(XMLWriter& parent)
    : out_stream_(parent.out_stream_),
      indent_(parent.indent_ + 2),
      closed_(true) {
    parent.FinishStartTag();
  }

  ~XMLWriter() {
    if (!closed_)
      CloseTag();
  }

  void StartTag(string tag) {
    closed_ = false;
    finished_start_tag_ = false;
    tag_ = tag;
    out_stream_ << string(indent_, ' ') << '<' << tag;
  }

  template<class C> void KeyVal(string key, C val) {
    out_stream_ << ' ' << key << "=\"" << val << '"';
  }

  void FinishStartTag() {
    if (finished_start_tag_)
      return;
    out_stream_ << ">\n";
    finished_start_tag_ = true;
  }

  void CloseTag() {
    if (!finished_start_tag_) {
      out_stream_ << "/>\n";
    } else {
      out_stream_ << string(indent_, ' ') << "</" << tag_ << ">\n";
    }
    closed_ = true;
  }

 private:
  ostream& out_stream_;
  int indent_;
  bool closed_;
  bool finished_start_tag_;
  string tag_;
};

class PepXMLResultsWriter : public ResultsWriter {
 public:
  PepXMLResultsWriter(const ProteinVec& proteins, const AuxLocVec& aux_locs,
		      ostream& out_stream);
  void WriteResults(pb::Results& pb_result);

 private:
  void WriteSpectrum(pb::Results& pb_result, XMLWriter& xml);
  void WriteMatch(const pb::Match& pb_match, int match,
		  double delta_cn, XMLWriter& xml);
  // void WriteProteinLocation(const pb::Protein& protein);

  const ProteinVec& proteins_;
  const AuxLocVec& aux_locs_;

  XMLWriter base_writer_;
  double neutral_mass_;
  int index_;
};

#endif // RESULTS_WRITER_H
