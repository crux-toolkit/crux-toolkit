/**
 * \file SpectrumCollectionFactory.h 
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 14 June 2011
 * \brief Return a SpectrumCollection object of the appropriate
 * derived class.
 */
#ifndef SPECTRUM_COLLECTION_FACTORY_H
#define SPECTRUM_COLLECTION_FACTORY_H

#include "SpectrumCollection.h"

/**
 * Instantiates a SpectrumCollection based on the extension of the
 * given file and the use-mstoolkit option.
 */
class SpectrumCollectionFactory {

 public:
  static Crux::SpectrumCollection* create(const char* filename);

};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif 
