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
#include "Percolator.h"
#include "QRanker.h"
#include "Barista.h"
#include "PrintProcessedSpectra.h"
#include "SearchForXLinks.h"
#include "ExtractColumns.h"
#include "SpectralCounts.h"
#include "ExtractRows.h"
#include "PrintVersion.h"
#include "StatColumn.h"
#include "SortColumn.h"

/**
 * The starting point for crux.  Prints a general usage statement when
 * given no arguments.  Runs one of the crux commands, including
 * printing the current version number.
 */
int main(int argc, char** argv){

  CruxApplicationList applications("crux");

  applications.add(new CreateIndex());

  applications.add(new MatchSearch());
  applications.add(new SequestSearch());

  applications.add(new ComputeQValues());
  applications.add(new Percolator());
  applications.add(new QRanker());

  applications.add(new Barista());

  applications.add(new PrintProcessedSpectra());

  applications.add(new SearchForXLinks());

  applications.add(new SpectralCounts());

  applications.add(new ExtractColumns());
  applications.add(new ExtractRows());
  applications.add(new StatColumn());
  applications.add(new SortColumn());

  applications.add(new PrintVersion());
  




  int ret = applications.main(argc, argv);
  return ret;

}// end main

