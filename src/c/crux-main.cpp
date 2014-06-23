/**
 * \file crux-main.cpp
 * AUTHOR: Barbara Frewen
 * CREATE DATE: November 24, 2008
 * \brief The starting point for the main crux program.
 *
 * Usage is "crux [command] [options] [arguments]" where command
 * is one of the primary crux commands.
 **/

#include "crux-main.h"
#include "crux-utils.h" // Need to get definition of NUM_FEATURES.

#include "CruxApplicationList.h"
#include "CreateIndex.h"
#include "MatchSearch.h"
#include "SequestSearch.h"
#include "ComputeQValues.h"
#include "ComputeQValuesLegacy.h"
#include "QRanker.h"
#include "Barista.h"
#include "PrintProcessedSpectra.h"
#include "GeneratePeptides.h"
#include "GetMs2Spectrum.h"
#include "PredictPeptideIons.h"
#include "SearchForXLinks.h"
#include "ExtractColumns.h"
#include "SpectralCounts.h"
#include "ExtractRows.h"
#include "PrintVersion.h"
#include "StatColumn.h"
#include "SortColumn.h"
#include "CruxHardklorApplication.h"
#include "CruxBullseyeApplication.h"
#include "GenerateDecoys.h"
#include "PercolatorApplication.h"
#include "MakePinApplication.h"
#include "TideIndexApplication.h"
#include "ReadSpectrumRecordsApplication.h"
#include "TideSearchApplication.h"
#include "CometApplication.h"
/**
 * The starting point for crux.  Prints a general usage statement when
 * given no arguments.  Runs one of the crux commands, including
 * printing the current version number.
 */
int main(int argc, char** argv){
  try {
#ifdef _MSC_VER
    // Turn off auto-tranlation of line-feed to 
    // carriage-return/line-feed
    _set_fmode(_O_BINARY);
#endif 

    CruxApplicationList applications("crux");

    // Primary commands
    applications.addMessage(applications.getListName() +
      " supports the following primary commands:");
    applications.add(new CruxBullseyeApplication());
    applications.add(new CreateIndex());
    applications.add(new TideIndexApplication());
    applications.add(new SequestSearch());
    applications.add(new MatchSearch());
    applications.add(new TideSearchApplication());
    applications.add(new ReadSpectrumRecordsApplication());
    applications.add(new CometApplication());
    applications.add(new PercolatorApplication());
    applications.add(new QRanker());
    applications.add(new Barista());
    applications.add(new SearchForXLinks());
    applications.add(new SpectralCounts());

    // Utilities
    applications.addMessage(applications.getListName() +
      " supports the following utility commands:");
    applications.add(new ComputeQValuesLegacy()); // deprecated name
    applications.add(new ComputeQValues());
    applications.add(new MakePinApplication());
    applications.add(new PredictPeptideIons());
    applications.add(new CruxHardklorApplication());
    applications.add(new PrintProcessedSpectra());
    applications.add(new GeneratePeptides());
    applications.add(new GenerateDecoys());
    applications.add(new GetMs2Spectrum());
    applications.add(new PrintVersion());

    // Utilities for processing tab-delimited text files
    applications.add(new ExtractColumns());
    applications.add(new ExtractRows());
    applications.add(new StatColumn());
    applications.add(new SortColumn());

    int ret = applications.main(argc, argv);
    return ret;
  } catch (const std::exception& e) {
    carp(CARP_FATAL, "An exception occurred: %s", e.what());
    return 1;
  } catch (...) {
    carp(CARP_FATAL, "An unknown exception occurred.");
    return 2;
  }

}// end main

