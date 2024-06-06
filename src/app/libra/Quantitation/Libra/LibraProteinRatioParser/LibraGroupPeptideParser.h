#ifndef LIBRA_GROUP_PEP_PARSER_H
#define LIBRA_GROUP_PEP_PARSER_H

#include "Parsers/Parser/Parser.h"

#include <vector>
#include <map>

#define FALSE 0
#define TRUE 1
#define dbug FALSE

#define MAX_LOG 3

class LibraGroupPeptideParser : public Parser {

 public:

  double minIntensity; // minimum intensity to use in peptide removal stage
  double tmpInf;   // dummy value used for infinite errors
  double tmpNoQuantitation; // dummy value used for no quantitation (all peptides were removed while processing)

  std::string current_peptide;
  std::string current_protein;

  LibraGroupPeptideParser(const char* xmlfile, Array<const char*>* peptides,
			  double minprob, int norm_channel);

  LibraGroupPeptideParser(Array<const char*>* xmlfiles, Array<const char*>* peptides, 
			  double minprob, int norm_channel);

  LibraGroupPeptideParser(const char* xmlfile, Array<const char*>* peptides,
			  double minprob, int norm_channel, double minimumThreshholdIntensity);

  LibraGroupPeptideParser(Array<const char*>* xmlfiles, Array<const char*>* peptides, 
			  double minprob, int norm_channel, double minimumThreshholdIntensity);
  
  LibraGroupPeptideParser(Array<const char*>* xmlfiles, const char * quantitationFileName,
			  double minprob, int norm_channel, double minimumThreshholdIntensity);


  ~LibraGroupPeptideParser();

  void setFilter(Tag* tag);

  double getRatioLogSum(int index);
  double getRatioSquareSum(int index);
  int getRatioNum();
  int getRatioNumFiltered();
  //  void getRatio();
  int getNumChannels();
  double getChannelMass(int index);

  Array<Tag*>* getProtXMLSummaryTags(double minpepprob, double minpepwt, double minprotprob);
  Array<Tag*>* getProtXMLTags(const char* prot_name, Array<const char*>* peptides);
  Array<Tag*>* getProtXMLTags();

  /**
   * a map to hold a protein's peptide's channel intensities
   * key = ratio_number (=peptide number within a protein), 
   * value = vector of peptide's channel intensities 
   */
  std::map<int, std::vector<double> > peptide_intensities_map;

  /**
   * a map to hold a protein's peptide's pre-normalized channel intensities
   * key = ratio_number (=peptide number within a protein), 
   * value = vector of peptide's pre-normalized channel intensities 
   */
  std::map<int, std::vector<double> > peptide_intensities_prenorm_map;

  /**
   * a map to hold a peptide's residue sequence
   * key = ratio_number (=peptide number within a protein), value = sequence
   */
  std::map<int, std::string > peptide_sequence_map;

  /**
   * A map to hold a protein's peptide's outlier flags 
   * key = ratio_number (=peptide number within a protein), value = vector of peptide's channel intensities.
   * A value of Y means kept, a value of N means removed from protein quantitation
   */
  std::map<int, std::string > peptide_keptFlag_map;


 protected:

  void parse(const char * xmlfile);
  Boolean peptideListMember(const char* pep);
  Boolean possiblePeptideListMember(const char* data);

  struct psm {
    std::string sequence;
    std::string protein;
    int num_channels;
    double* intensities;
  };

  Array<const char *>* peptides_;
  map <string, vector<psm> > all_peptides_map_; // used when parsing the entire pepXML file(s) only once


  //  int num_channels_;

  /**
   * user selected channel to use as reference for normalization
   * (0 means no normalization).
   */
  int norm_channel_;

  Array<double>* channel_masses_; // array of m/z's from the condition file

  double* ratio_log_sum_;

  double* ratio_sum_;

  double* protein_mass_array_; // array of found reagent m/z's for given peptide
  double* protein_ratio_array_;
  double* protein_stdev_array_;

  double* protein_se_array_; // holds standard error array (which is standard dev / sqrt(n))
  double* protein_variance_array_;
  double* protein_ratio_array_wrt_refchannel_;
  double* protein_error_array_wrt_refchannel_;

  int ratio_num_; // number of peptides for a given protein


  /** number of protein's peptides that pass all filters (such as outlier 
   * removal and minimum intensity filters)  
   */
  int ratio_num_filtered_;

  double min_probability_;

  Array<Tag*>* summary_tags_;

  //double inverse_ratio_square_sum_;
  //  int ratio_num_;
  //Boolean heavy2light_;

  Array<const char *>* xmlfiles_;
  const char * quantitation_file_name_;

  Boolean single_input_;
  Boolean parse_all_;

  void handleSingleXmlFileAttributes( const char * xmlFile );

  /**
   * initialize class attributes
   */
  void initializeClassAttributes( Array<const char*>* peptides, double minprob,
				  int norm_channel, double thresh);

  /**
   * initialize storage classes used in parse
   *@param number of reagent lines in a spectrum
   */
  void clearProteinRatioData();
  void initializeStorageWithinParse ( int nChan );
  void storeAndNormalizeValues( double* intensities, int nChan );
  void calculateBasicStats( int nChan );

  /**
   * Re-calculate Ratios by not including peptides that deviate from the mean
   * by more than 2 sigma.
   * If other attributes have been set, such as a minimum threshhold intensity,
   * that's handled here too.
   * Stores kept or removed flag for each peptide.  
   *@param number of reagent lines in a spectrum
   */  
  void recalculateRatios( int nChan );


  /**  
   * Appends peptide and protein quantitation to existing output file
   * called quantitation.tsv.  The outfile is created by LibraProteinRatioParser
   * and appended to here.  Started to create a class to handle the File
   * writing, but haven't finished.  It's the QuantitationFile in this package.
   * wrote ability to read a name other than quantitation.tsv in the LibraConditionHandler2
   * so that can be read used by these classes in the future.
   *@param number of reagent lines in a spectrum
   */  
  void writePeptideAndProteinQuantitationToOutfile( int nChan );

};


#endif
