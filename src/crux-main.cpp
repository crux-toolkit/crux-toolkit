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
#include "util/crux-utils.h" // Need to get definition of NUM_FEATURES.

#include "app/xlink/xlink_assign_ions.h"
#include "app/xlink/xhhc_score_peptide_spectrum.h"
#include "app/xlink/xlink_search.h"
#include "app/CruxApplicationList.h"
#include "app/ComputeQValues.h"
#include "app/ComputeQValuesLegacy.h"
#include "app/CreateDocs.h"
#include "app/qranker-barista/QRanker.h"
#include "app/qranker-barista/Barista.h"
#include "app/PrintProcessedSpectra.h"
#include "app/GeneratePeptides.h"
#include "app/GetMs2Spectrum.h"
#include "app/LocalizeModification.h"
#include "app/ParamMedicApplication.h"
#include "app/Pipeline.h"
#include "app/PredictPeptideIons.h"
#include "app/xlink/SearchForXLinks.h"
#include "app/ExtractColumns.h"
#include "app/SpectralCounts.h"
#include "app/ExtractRows.h"
#include "app/PrintVersion.h"
#include "app/StatColumn.h"
#include "app/SortColumn.h"
#include "app/hardklor/CruxHardklorApplication.h"
#include "app/bullseye/CruxBullseyeApplication.h"
#include "app/PercolatorApplication.h"
#include "app/MakePinApplication.h"
#include "app/TideIndexApplication.h"
#include "app/ReadSpectrumRecordsApplication.h"
#include "app/ReadTideIndex.h"
#include "app/TideSearchApplication.h"
#include "app/CometApplication.h"
#include "app/PSMConvertApplication.h"
#include "app/CascadeSearchApplication.h"
#include "app/AssignConfidenceApplication.h"
#include "app/SubtractIndexApplication.h"
/**
 * The starting point for crux.  Prints a general usage statement when
 * given no arguments.  Runs one of the crux commands, including
 * printing the current version number.
 */
int main(int argc, char** argv) {
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
    applications.add(new TideIndexApplication());
    applications.add(new TideSearchApplication());
    applications.add(new ReadSpectrumRecordsApplication());
    applications.add(new ReadTideIndex());
    applications.add(new CometApplication());
    applications.add(new PercolatorApplication());
    applications.add(new QRanker());
    applications.add(new Barista());
    applications.add(new SearchForXLinks());
    applications.add(new SpectralCounts());
    applications.add(new PipelineApplication());
    applications.add(new CascadeSearchApplication());
    applications.add(new AssignConfidenceApplication());

    // Utilities
    applications.addMessage(applications.getListName() +
      " supports the following utility commands:");
    applications.add(new MakePinApplication());
    applications.add(new PredictPeptideIons());
    applications.add(new CruxHardklorApplication());
    applications.add(new ParamMedicApplication());
    applications.add(new PrintProcessedSpectra());
    applications.add(new GeneratePeptides());
    applications.add(new GetMs2Spectrum());
    applications.add(new CreateDocs());
    applications.add(new PrintVersion());
    applications.add(new PSMConvertApplication());
    applications.add(new SubtractIndexApplication());
    applications.add(new XLinkAssignIons());
    applications.add(new XLinkScoreSpectrum());
    applications.add(new LocalizeModificationApplication());

    // Utilities for processing tab-delimited text files
    applications.add(new ExtractColumns());
    applications.add(new ExtractRows());
    applications.add(new StatColumn());
    applications.add(new SortColumn());

    int ret = applications.main(argc, argv);
    google::protobuf::ShutdownProtobufLibrary();
    return ret;
  } catch (const std::exception& e) {
    carp(CARP_FATAL, "An exception occurred: %s", e.what());
  } catch (...) {
    carp(CARP_FATAL, "An unknown exception occurred.");
  }
}// end main

