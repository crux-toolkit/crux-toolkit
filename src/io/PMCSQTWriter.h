#ifndef PMCSQTWRITER_H
#define PMCSQTWRITER_H

#include "model/PeptideMatch.h"
#include "model/ProteinMatchCollection.h"
#include "PSMWriter.h"
#include "model/SpectrumMatch.h"
#include "SQTWriter.h"

using namespace std;

class PMCSQTWriter : public SQTWriter, public PSMWriter {

 public:

  void openFile(
    CruxApplication* application,
    string filename,
    MATCH_FILE_TYPE type
  );

  void write(
    MatchCollection* collection,
    string database
  );

  void closeFile();

  /**
   * Writes the data in a ProteinMatchCollection to the currently open file
   */
  void write(
    ProteinMatchCollection* collection, ///< collection to be written
    string database ///< the database name
//    int top_match ///< the top matches to output
  );

  string getSpectrumTitle(
    int spectrum_scan_number,
    int charge
  );

 protected:

  /**
   * Writes the PSMs in a ProteinMatchCollection to the currently open file
   */
  void writePSMs(
    ProteinMatchCollection* collection ///< collection to be written
//    int top_match //< the top matches to output
  );

};

#endif // PMCSQTWRITER_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

