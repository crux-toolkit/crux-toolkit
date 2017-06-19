/**
 * \file SpectrumCollectionFactory.cpp 
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 14 June 2011
 * \brief Return a SpectrumCollection object of the appropriate
 * derived class.
 */
#include "parameter.h"
#include "SpectrumCollectionFactory.h"
#include "MSToolkitSpectrumCollection.h"
#include "PWIZSpectrumCollection.h"
#include "SpectrumRecordSpectrumCollection.h"
#include "util/FileUtils.h"
#include "util/Params.h"

/**
 * Instantiates a SpectrumCollection based on the extension of the
 * given file and the use-mstoolkit and msgf options.
 */
Crux::SpectrumCollection* SpectrumCollectionFactory::create(const string& filename) {
  if (!FileUtils::Exists(filename)) {
    carp(CARP_FATAL, "The file %s does not exist. \n", filename.c_str());
  }
  if (FileUtils::IsDir(filename)) {
    carp(CARP_FATAL, "Path %s is a directory. \n Please enter a spectrum filename\
(.ms2, .mgf, or .mzXML)", filename.c_str());
  }

  if (SpectrumRecordSpectrumCollection::IsSpectrumRecordFile(filename)) {
    return new SpectrumRecordSpectrumCollection(filename);
  }

  string parser = Params::GetString("spectrum-parser");
  if (parser == "pwiz") {
    carp(CARP_DEBUG, "Using protewizard to parse spectra");
    return new PWIZSpectrumCollection(filename);
  } else if (parser == "mstoolkit") {
    carp(CARP_DEBUG, "Using mstoolkit to parse spectra");
    return new MSToolkitSpectrumCollection(filename);
  }

  carp(CARP_FATAL, "Unknown spectrum parser type");
  return(NULL); // Avoid compiler warning.
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
