#ifndef PMCPEPXMLWRITER_H
#define PMCPEPXMLWRITER_H

#include "carp.h"
#include "model/Peptide.h"
#include "model/PeptideMatch.h"
#include "PepXMLWriter.h"
#include "model/Protein.h"
#include "model/ProteinMatch.h"
#include "model/ProteinMatchCollection.h"
#include "model/Spectrum.h"
#include "model/SpectrumMatch.h"

using namespace std;

class PMCPepXMLWriter : public PepXMLWriter {

 public:
  /**
   * Writes the data in a ProteinMatchCollection to the currently open file
   */
  void write(
    ProteinMatchCollection* collection ///< collection to be written
  );

 protected:

  /**
   * Writes the PSMs in a ProteinMatchCollection to the currently open file
   */
  void writePSMs(
    ProteinMatchCollection* collection ///< collection to be written
  );

};

#endif // PMCPEPXMLWRITER_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

