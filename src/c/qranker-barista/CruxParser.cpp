/*************************************************************************//**
 * \file CruxParser.cpp
 * AUTHOR: Manijeh Naseri and Sean McIlwain
 * CREATE DATE: 06/06/2012
 * \brief Parses tab-delimited search result files
 ****************************************************************************/
#include "CruxParser.h"
#include "Peptide.h"
#include "MatchFileReader.h"
#include "PeptideSrc.h"
#include "Match.h"
#include "MatchColumns.h"


using namespace std; 
/******************************/


/**
 * Default constructor
 */
CruxParser :: CruxParser() 
  : SQTParser(){
}

/**
 * Default destructor
 */
CruxParser :: ~CruxParser(){
}

 
/**
 * Parse tab delimited file.
 * Generates the same QRanker internal tables.
 * Set the matches in the tab delimited file. 
 */
void CruxParser ::readMatches(
  MatchFileReader& reader, ///< Reader for the delimited file.
  int final_hits,  ///< Total number of matches
  enzyme enz, ///< Enzyme used in search
  bool decoy ///< Are all the matches decoy?
  ) {
  
  //set the matches to zero 
  int hits_read = 0;

  // read file line by line and set the variables  
  while (reader.hasNext()) {  
    
    int current_scan = reader.getInteger(SCAN_COL);
    int current_charge = reader.getInteger(CHARGE_COL);
    sqt_match matches;
    matches.scan = current_scan;
    matches.charge = current_charge;
    matches.precursor_mass = reader.getDouble(SPECTRUM_NEUTRAL_MASS_COL);
    
    if (!reader.empty(DISTINCT_MATCHES_SPECTRUM_COL)) {
      matches.num_sequence_comparisons = reader.getInteger(DISTINCT_MATCHES_SPECTRUM_COL);
    } else if (!reader.empty(MATCHES_SPECTRUM_COL)) {
      matches.num_sequence_comparisons = reader.getInteger(MATCHES_SPECTRUM_COL);
    } else {
      carp_once(CARP_WARNING, "Empty Matches/Spectrum col");
      matches.num_sequence_comparisons = 0;
    }
    hits_read = 0;

    //iterate over all rows in delimited file, fill-in sqt match structure.
    while (reader.hasNext() && 
           reader.getInteger(SCAN_COL) == current_scan && 
           reader.getInteger(CHARGE_COL) == current_charge) {

      //fill in sqt_match structure with correct fields from delimited file.
      matches.xcorr_rank.push_back(reader.getInteger(XCORR_RANK_COL));
      matches.sp_rank.push_back(reader.getInteger(SP_RANK_COL));
      matches.delta_cn.push_back(reader.getDouble(DELTA_CN_COL));
      matches.sp_score.push_back(reader.getDouble(SP_SCORE_COL));
      matches.calc_mass.push_back(reader.getDouble(PEPTIDE_MASS_COL));
      matches.xcorr_score.push_back(reader.getDouble(XCORR_SCORE_COL));
      matches.num_ions_matched.push_back(reader.getInteger(BY_IONS_MATCHED_COL));
      matches.num_total_ions.push_back(reader.getFloat(BY_IONS_TOTAL_COL)); 
        
      //set sequence_id 
      string sequence = reader.getString(SEQUENCE_COL);
      //set flanking_aa 
      string flanking_aa = reader.getString(FLANKING_AA_COL); 
      
      if (flanking_aa.length() != 2)
      {
        carp(CARP_DEBUG,
            "Flanking AA value length is expected to be 2, "
            "but was %d", flanking_aa.length());
      }

      ostringstream oss;

      oss << flanking_aa[0] << "." << sequence << "." << flanking_aa[1];
      string sqt_sequence = oss.str();
      matches.peptides.push_back(sqt_sequence);
      
      //read protein Id from search result 
      vector<string> protein_ids;
      reader.getStringVectorFromCell(
        PROTEIN_ID_COL, 
        protein_ids, 
        ','
      ); 
      matches.num_proteins_in_match.push_back(protein_ids.size());
      for (size_t idx = 0;idx < protein_ids.size(); idx++) {
        string protein_id = protein_ids.at(idx);
        //change the format of protein_id to be same as SQTParser
        size_t protein_idx_pos = protein_id.rfind('(');
        if (protein_idx_pos != string::npos) {
          string pep_pos= protein_id.substr(protein_idx_pos);
          pep_pos=pep_pos.substr(1,pep_pos.size()-2);
          matches.peptide_pos.push_back(atoi(pep_pos.c_str()));
          protein_id.erase(protein_idx_pos);
        } else {
          matches.peptide_pos.push_back(-1);
        }
        matches.proteins.push_back(protein_id);
      }       
      hits_read++;
      reader.next();  
    }
  
    //call SQTParser method to generate internal tables
    add_matches_to_tables(
      matches, 
      decoy_prefix,
      hits_read, 
      final_hits,
      decoy
    );

    //call SQTParser method to generate feature vectors
    extract_features(
      matches,  
      hits_read, 
      final_hits,
      enz 
    );
    
  }
}

 /*
  *gets the path of delimited file 
  *\returns true if it can open the file 
  */
bool CruxParser :: read_search_results(
  string& cur_fname, /// < current file path to parse
  bool decoy ///< is this a file of decoys?
  ) {

  //read file 
  MatchFileReader reader(cur_fname);

  //TODO (SJM) - test to see if we opened the file, if not, WARNING, then return false
 
  readMatches(
    reader, 
    fhps,  
    e,
    decoy
  );
  return true; 
}

/**
 *\returns txt file extension  
 */
std::string CruxParser::get_parser_extension() {
  return ".txt";
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
