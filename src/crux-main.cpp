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

#include "app/CruxApplicationList.h"
#include "app/ComputeQValues.h"
#include "app/ComputeQValuesLegacy.h"
#include "app/CreateDocs.h"
#include "app/PrintProcessedSpectra.h"
#include "app/GeneratePeptides.h"
#include "app/GetMs2Spectrum.h"
#include "app/LocalizeModification.h"
#include "app/ParamMedicApplication.h"
#include "app/Pipeline.h"
#include "app/PredictPeptideIons.h"
#include "app/SpectralCounts.h"
#include "app/PrintVersion.h"
#include "app/hardklor/CruxHardklorApplication.h"
#include "app/bullseye/CruxBullseyeApplication.h"
#include "app/PercolatorApplication.h"
#include "app/MakePinApplication.h"
#include "app/TideIndexApplication.h"
#include "app/ReadSpectrumRecordsApplication.h"
#include "app/ReadTideIndex.h"
#include "app/TideSearchApplication.h"
#include "app/CometApplication.h"
#include "app/CometIndexApplication.h"
#include "app/PSMConvertApplication.h"
#include "app/CascadeSearchApplication.h"
#include "app/AssignConfidenceApplication.h"
#include "app/SubtractIndexApplication.h"
#include "app/SpectrumConvertApplication.h"
#include "app/DIAmeterApplication.h"
#include "app/KojakApplication.h"

/* Code addded by Rufino*/
#include "app/CruxApplication.h"

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
    applications.add(new SpectrumConvertApplication());
    applications.add(new ReadSpectrumRecordsApplication());
    applications.add(new ReadTideIndex());
    applications.add(new CometApplication());
    applications.add(new CometIndexApplication());
    applications.add(new PercolatorApplication());
    applications.add(new SpectralCounts());
    applications.add(new PipelineApplication());
    applications.add(new CascadeSearchApplication());
    applications.add(new AssignConfidenceApplication());
    applications.add(new KojakApplication());

    applications.add(new DIAmeterApplication());

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
    applications.add(new LocalizeModificationApplication());


    /* Code added by Rufino */
    /* Pointers to the first and last element in the applications list, passed to the static function setApplicationsList
    that will be used in CruxApplication.cpp to locate a specific parameter assinged to an appplication. */
    CruxApplication::setApplicationsList(applications.begin(), applications.end());
    

    int ret = applications.main(argc, argv);
    google::protobuf::ShutdownProtobufLibrary();
    return ret;
  } catch (const std::exception& e) {
    carp(CARP_FATAL, "An exception occurred: %s", e.what());
  } catch (...) {
    carp(CARP_FATAL, "An unknown exception occurred.");
  }
}// end main

