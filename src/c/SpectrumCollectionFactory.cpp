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
#include <sys/types.h>
#include <sys/stat.h>

/**
 * Instantiates a SpectrumCollection based on the extension of the
 * given file and the use-mstoolkit and msgf options.
 */
Crux::SpectrumCollection* SpectrumCollectionFactory::create(const char* filename){
  Crux::SpectrumCollection* collection = NULL;
  struct stat stat_buff ; 
  stat(filename, &stat_buff);
  if(stat(filename, &stat_buff)!=0){
    carp(CARP_FATAL, "The file %s does not exist. \n", filename);
  }
  if (S_ISDIR(stat_buff.st_mode)){
    carp(CARP_FATAL, "Path %s is a directory. \n Please enter a spectrum filename\
(.ms2,.mgf, or .mzXML)",filename);
  }


  switch (get_spectrum_parser_parameter("spectrum-parser")) {
    case PROTEOWIZARD_SPECTRUM_PARSER:
      carp(CARP_DEBUG, "Using protewizard to parse spectra");
      collection = new PWIZSpectrumCollection(filename);
      break;
    case MSTOOLKIT_SPECTRUM_PARSER:
      carp(CARP_DEBUG, "Using mstoolkit to parse spectra");
      collection = new MSToolkitSpectrumCollection(filename);
      break;
    case INVALID_SPECTRUM_PARSER:
    case NUMBER_SPECTRUM_PARSERS:
      carp(CARP_FATAL, "Unknown spectrum parser type");
  
  }

  return collection;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
