#ifndef PMCSQTWRITER_H
#define PMCSQTWRITER_H

#include "PeptideMatch.h"
#include "ProteinMatchCollection.h"
#include "SpectrumMatch.h"
#include "SQTWriter.h"

using namespace std;

class PMCSQTWriter : public SQTWriter {

 public:
  /**
   * Writes the data in a ProteinMatchCollection to the currently open file
   */
  void write(
    ProteinMatchCollection* collection, ///< collection to be written
    string database, ///< the database name
    int top_match ///< the top matches to output
  );

 protected:

  /**
   * Writes the PSMs in a ProteinMatchCollection to the currently open file
   */
  void writePSMs(
    ProteinMatchCollection* collection, ///< collection to be written
    int top_match //< the top matches to output
  );

};

#endif // PMCSQTWRITER_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

