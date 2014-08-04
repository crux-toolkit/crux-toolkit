// Benjamin Diament
//
// Add to the index of peptide records the pre-computed theoretical peaks.
// We store the TheoreticalPeakSetDiff (q.v.) for each peptide. 
//
// Example command-line:
// peptide_peaks --proteins=<raw_proteins.proto input file> \
//               --in_peptides=<peptides.proto input file> \
//               --out_peptides=<peptides.proto output file> \
// 
// Rather than store the code for each peak in the diff set, we store the
// delta between each pair. This gives us a bit of compression since the
// values are stored as varint. The peak locations then have to be restored
// at search time.
//
// TODO 248: We're only doing this to guarantee the exact same results as Crux
// used to return, but perhaps the diffs don't really add useful info, in which 
// case we could eliminate them.

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include "records.h"
#include "peptide.h"
#include "theoretical_peak_set.h"
#include "abspath.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)
/*
static void AddPeaksToPB(pb::Peptide* peptide, const TheoreticalPeakArr* peaks,
			 int charge, bool neg) {
  int last_code = 0;
  TheoreticalPeakArr::const_iterator i = peaks->begin();
  for (; i != peaks->end(); ++i) {
    int delta = i->Code() - last_code;
    last_code = i->Code();
    if (neg) {
      if (charge == 1) {
        peptide->add_neg_peak1(delta);
      } else {
        peptide->add_neg_peak2(delta);
      }
    } else {
      if (charge == 1) {
        peptide->add_peak1(delta);
      } else {
        peptide->add_peak2(delta);
      }
    }
  }
}
*/

void AddTheoreticalPeaks(const vector<const pb::Protein*>& proteins,
			 const string& input_filename,
			 const string& output_filename) {
  pb::Header orig_header, new_header;
  HeadedRecordReader reader(input_filename, &orig_header);
  CHECK(orig_header.file_type() == pb::Header::PEPTIDES);
  CHECK(orig_header.has_peptides_header());
//  MassConstants::Init(&orig_header.peptides_header().mods());
  new_header.set_file_type(pb::Header::PEPTIDES);
  pb::Header_PeptidesHeader* subheader = new_header.mutable_peptides_header();
  subheader->CopyFrom(orig_header.peptides_header());
  subheader->set_has_peaks(true);
  pb::Header_Source* source = new_header.add_source();
  source->mutable_header()->CopyFrom(orig_header);
  source->set_filename(AbsPath(input_filename));
  HeadedRecordWriter writer(output_filename, new_header);
  CHECK(reader.OK());
  CHECK(writer.OK());

  pb::Peptide pb_peptide;
//  const int workspace_size = 2000; // More than sufficient for theor. peaks.
//  TheoreticalPeakSetDiff workspace(workspace_size);
  while (!reader.Done()) {
    reader.Read(&pb_peptide);
/*    Peptide peptide(pb_peptide, proteins);
    workspace.Clear();
    peptide.ComputeTheoreticalPeaks(&workspace);
    TheoreticalPeakArr peaks_charge_1(2000);
    TheoreticalPeakArr peaks_charge_2(2000);
    TheoreticalPeakArr negs_charge_1(2000);
    TheoreticalPeakArr negs_charge_2(2000);
    workspace.GetPeaks(&peaks_charge_1, &negs_charge_1,
		       &peaks_charge_2, &negs_charge_2, NULL);
    AddPeaksToPB(&pb_peptide, &peaks_charge_1, 1, false);
    AddPeaksToPB(&pb_peptide, &peaks_charge_2, 2, false);
    AddPeaksToPB(&pb_peptide, &negs_charge_1, 1, true);
    AddPeaksToPB(&pb_peptide, &negs_charge_2, 2, true);
*/    CHECK(writer.Write(&pb_peptide));
  }
  CHECK(reader.OK());
}
