// Benjamin Diament
//
// WARNING: When we filter modified peptides, we end up generating a new aux
// location for every single modification of the peptide that makes it through
// the filter. If you don't want to do this, please re-write the test!!

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <bitset>
#include <gflags/gflags.h>
#include "records.h"
#include "records_to_vector-inl.h"
#include "raw_proteins.pb.h"
#include "peptides.pb.h"
#include "mass_constants.h"
#include "abspath.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)

DEFINE_string(proteins, "", "File of raw proteins, raw_proteins.proto");
DEFINE_string(in_peptides, "", "File of peptides to filter");
DEFINE_string(out_peptides, "", "File of peptides to create");
DEFINE_string(in_aux_locations, "", "File of auxiliary locations associated with in_peptides");
DEFINE_string(out_aux_locations, "", "File of auxiliary locations associated with out_peptides");
DEFINE_int32(max_missed_cleavages, 0, "Number of missed cleavages to allow");
DEFINE_string(digestion, "full", "May be 'full' or 'partial'");

bool FewMissedCleavages(const char* residues, int len) {
  int n = 0;
  for (int i = 0; i < len-1; ++i)
    if ((residues[i] == 'K' || residues[i] == 'R') && residues[i+1] != 'P')
      if (++n > FLAGS_max_missed_cleavages)
	return false;

  return true;
}

// FULLY TRYPTIC, ALLOW MISSED CLEAVAGES
bool FilterFull(const pb::Peptide& peptide, const pb::Protein& protein, int pos) {
  int len = peptide.length();
  const char* residues = protein.residues().data() + pos;
  if (pos + len < protein.residues().size()) { // not at protein end
    if (residues[len-1] != 'K' && residues[len-1] != 'R')
      return false;
    if (residues[len] == 'P')
      return false;
  }
  if (pos > 0) { // not at protein start
    if (residues[0] == 'P')
      return false;
    if (residues[-1] != 'K' && residues[-1] != 'R')
      return false;
  }

  return FewMissedCleavages(residues, len);
}

// PARTIALLY TRYPTIC, ALLOW MISSED CLEAVAGES
bool FilterPartial(const pb::Peptide& peptide, const pb::Protein& protein, int pos) {
  int len = peptide.length();
  const char* residues = protein.residues().data() + pos;
  if (pos == 0 || pos + len == protein.residues().size()) // at protein ends
    return FewMissedCleavages(residues, len);
  if ((residues[len-1] == 'K' || residues[len-1] == 'R') && residues[len] != 'P')
    return FewMissedCleavages(residues, len);
  if ((residues[-1] == 'K' || residues[-1] == 'R') && residues[0] != 'P')
    return FewMissedCleavages(residues, len);

  return false; // no tryptic edge
}

typedef bool (*FilterFunc)(const pb::Peptide& peptide, const pb::Protein& protein, int pos);

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  vector<const pb::Protein*> proteins;
  cerr << "Reading Proteins" << endl;
  CHECK(ReadRecordsToVector<pb::Protein>(&proteins, FLAGS_proteins));
  
  vector<const pb::AuxLocation*> aux_locs;
  pb::Header orig_aux_locs_header;
  cerr << "Reading auxiliary locations" << endl;
  CHECK(ReadRecordsToVector<pb::AuxLocation>(&aux_locs, 
    FLAGS_in_aux_locations, &orig_aux_locs_header));

  pb::Header orig_header, new_header;
  HeadedRecordReader reader(FLAGS_in_peptides, &orig_header);
  CHECK(orig_header.file_type() == pb::Header::PEPTIDES);
  CHECK(orig_header.has_peptides_header());
  CHECK(orig_header.peptides_header().enzyme() == "none"
   || orig_header.peptides_header().enzyme() == "trypsin");

  new_header.set_file_type(pb::Header::PEPTIDES);
  pb::Header_PeptidesHeader* subheader = new_header.mutable_peptides_header();
  subheader->CopyFrom(orig_header.peptides_header());
  // describe new filtering criteria in output file header
  subheader->set_enzyme("trypsin");
  subheader->set_full_digestion(false);
  subheader->set_max_missed_cleavages(FLAGS_max_missed_cleavages);
  pb::Header_Source* source = new_header.add_source();
  source->mutable_header()->CopyFrom(orig_header);
  source->set_filename(AbsPath(FLAGS_in_peptides));

  HeadedRecordWriter peptide_writer(FLAGS_out_peptides, new_header);
  CHECK(peptide_writer.OK());

  // Create the new auxiliary locations header
  pb::Header out_aux_locs_header;
  out_aux_locs_header.set_file_type(pb::Header::AUX_LOCATIONS);  
  
  // Add the new peptides file as a source
  pb::Header_Source* aux_locs_source = out_aux_locs_header.add_source();
  aux_locs_source->set_filename(AbsPath(FLAGS_out_peptides));
  aux_locs_source->mutable_header()->CopyFrom(new_header);

  // Add the old auxiliary locations file as a source
  aux_locs_source = out_aux_locs_header.add_source();
  source->set_filename(AbsPath(FLAGS_in_aux_locations));
  aux_locs_source->mutable_header()->CopyFrom(orig_aux_locs_header);

  // Create the new auxiliary locations writer
  HeadedRecordWriter aux_loc_writer(FLAGS_out_aux_locations, out_aux_locs_header);
  CHECK(aux_loc_writer.OK());


  CHECK(FLAGS_digestion == "full" || FLAGS_digestion == "partial");
  FilterFunc filter_func = (FLAGS_digestion == "full") ? &FilterFull : &FilterPartial;

  pb::Peptide peptide;
  int count = 0;
  int id = 0;
  int out_count = 0;
  int out_aux_loc_idx = 0;
  while (!reader.Done()) {
    reader.Read(&peptide);    

    // The first location can be found as part of the pb::peptide. We need to
    // check to see if this first location makes it through the filter.
    ++count;
    const pb::Protein* first_protein = proteins[peptide.first_location().protein_id()];
    int first_pos = peptide.first_location().pos();
    if (filter_func(peptide, *first_protein, first_pos)) {
      ++out_count;
    } else {
      peptide.clear_first_location();
    }

    // The rest of the locations are in the aux locations (if the current
    // peptide has any). We need to go through each aux location to see if it 
    // makes it through the filter. Furthermore, if the original first location
    // from above did NOT make it through the filter, we need to find a NEW 
    // first location; in this case, the FIRST AUX LOCATION becomes the NEW 
    // FIRST LOCATION.
    pb::AuxLocation out_aux_loc;
    if (peptide.has_aux_locations_index()) {
      // Grab the index into the auc_locs vector
      int aux_loc_idx = peptide.aux_locations_index();
      
      // Use the index to get the set of aux locations from the aux_locs vector
      const pb::AuxLocation* aux_loc = aux_locs[aux_loc_idx];
      
      // Go through each aux location in the set of aux locations
      for (int loc = 0; loc < aux_loc->location_size(); ++loc) {
        ++count;

        // Check to see if this aux location makes it through the filter
        const pb::Protein* protein = proteins[aux_loc->location(loc).protein_id()];
        int pos = aux_loc->location(loc).pos();
        if (filter_func(peptide, *protein, pos)) {
          ++out_count;

          // This aux location made it through the filter. Check to see if we
          // should add it as the NEW first location, or just as a regular 
          // aux location.
          if (!peptide.has_first_location()) {
            // If the peptide does not have a first location, set this aux
            // location as the new first location.
            peptide.mutable_first_location()->set_protein_id(aux_loc->location(loc).protein_id());
            peptide.mutable_first_location()->set_pos(aux_loc->location(loc).pos());
          } else {
            // The peptide has a first location already, so just add this
            // aux location as a regular aux location.
            out_aux_loc.add_location()->CopyFrom(aux_loc->location(loc));
          }
        }
      }
    }

    // For a peptide to have made it through the filter, it MUST have a first
    // location within a protein.
    if (peptide.has_first_location()) {
      if (out_aux_loc.location_size() > 0) {
        // The peptide has some aux locations that made it through the filter,
        // and it has a new aux location index we need to point to.
        peptide.set_aux_locations_index(out_aux_loc_idx++);
        aux_loc_writer.Write(&out_aux_loc);          
      } else {
        // The peptide does NOT have any aux locations, so clear the aux
        // locations index since it is no longer valid
        peptide.clear_aux_locations_index();
      }

      peptide.set_id(id++);
      peptide_writer.Write(&peptide);
    }
    if (count % 100000 == 0)
      cerr << "Read " << count << " peptides/locations; "
      << "Wrote " << out_count << " peptides/locations.\n";
  }
  cout << "Totals: Read " << count << " peptides/locations; "
       << "Wrote " << out_count << " peptides/locations.\n";

  CHECK(reader.OK());

  return 0;
}
