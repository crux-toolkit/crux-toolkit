/**
 * \file XLinkMatch.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Object for Defining a Match in an xlink search
 *****************************************************************************/
#ifndef XLINKMATCH_H_
#define XLINKMATCH_H_

#include "crux-utils.h"
#include "XLink.h"
#include "XLinkMatchCollection.h"
#include <string>
#include "Match.h"


class XLinkMatch : public Crux::Match {

 protected:
  XLinkMatchCollection* parent_; ///< Owner of this match
  FLOAT_T pvalue_; ///< p-value of the match
  bool mass_calculated_[NUMBER_MASS_TYPES]; ///< is mass calculated?
  FLOAT_T mass_[NUMBER_MASS_TYPES]; ///<calculated mass

 public:
  
  /**
   *  Constructor for XLinkMatch()
   */
  XLinkMatch();

  /**
   * Default destructor for XLinkMatch
   */
  virtual ~XLinkMatch();

  virtual XLINKMATCH_TYPE_T getCandidateType() = 0;
  virtual int getNumMissedCleavages() = 0;
  virtual bool isModified() = 0;
  virtual std::string getSequenceString() = 0;
  virtual FLOAT_T calcMass(MASS_TYPE_T mass_type) = 0;
  virtual XLinkMatch* shuffle() = 0;
  virtual void predictIons(IonSeries* ion_series, int charge)=0;
  virtual std::string getIonSequence(Ion* ion)=0;
  virtual Crux::Peptide* getPeptide(int peptide_idx)=0;

  /**
   * \returns the mass of the match
   */
  FLOAT_T getMass(
    MASS_TYPE_T mass_type /// MONO or AVERAGE?
  );

  std::string getCandidateTypeString();

  
  void decrementPointerCount();
  
  /**
   * computes the pvalue for this match using the provided weibull paramters
   */
  void computeWeibullPvalue(
    FLOAT_T shift, ///< shift parameter for weibull
    FLOAT_T eta, ///< eta parameter for weibull
    FLOAT_T beta ///< beta parameter for weibull
    );

  /**
   * \returns the mass error in part-per-million (ppm)
   */
  FLOAT_T getPPMError();

  /**
   * sets the XLinkMatchCollection owner of the match
   */
  void setParent(XLinkMatchCollection* parent);

  /**
   * \returns the protein id string for this match
   */
  virtual std::string getProteinIdString();
  
  /**
   * \returns the protein id string for this match.  In the case
   * of a xlinked peptide, reportes the position in the protein
   */
  virtual std::string getProteinIdXString();

  /**
   *\returns the flanking amino acids for the match
   */
  virtual std::string getFlankingAAString();

  
  /**
   * Print one field in the tab-delimited output file, based on column index.
   * overridden from Match
   */
  virtual void printOneMatchField(
    int      column_idx,             ///< Index of the column to print. -in
    MatchCollection* collection,  ///< collection holding this match -in 
    MatchFileWriter*    output_file,            ///< output stream -out
    int      scan_num,               ///< starting scan number -in
    FLOAT_T  spectrum_precursor_mz,  ///< m/z of spectrum precursor -in
    int      num_target_matches,            ///< target matches in spectrum -in
    int      num_decoy_matches,      ///< decoy matches (if any) for this spectrum -in
    int      b_y_total,              ///< total b/y ions -in
    int      b_y_matched             ///< Number of b/y ions matched. -in
  );    

  /**
   *\returns the string value of the given candidate type
   */
  static std::string getCandidateTypeString(
    XLINKMATCH_TYPE_T candidate ///< candidate
  );

  /**
   * \returns the candidate type from a string value
   */
  static XLINKMATCH_TYPE_T getCandidateType(
    std::string& candidate ///< candidate in string format
  );


};

#endif

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
