/**
 * \file print-processed-spectra.h
 *
 * AUTHOR: Barbara Frewen
 * CREATE DATE: September 18, 2009
 * DESCRIPTION: Main method for the print-processed-spectra command.
 *              For every spectrum in an ms2 file, process as for
 *              xcorr and print peaks in ms2 format to new file.
 * REVISION:
 */

#ifndef PRINT_PROCESSED_SPECTRA_H
#define PRINT_PROCESSED_SPECTRA_H

#include "crux-utils.h"
#include "carp.h"
#include "parameter.h"
#include "SpectrumCollection.h"
#include "FilteredSpectrumChargeIterator.h"
#include "scorer.h"

int print_processed_spectra_main(int argc, char** argv);


#endif //PRINT_PROCESSED_SPECTRA_H
