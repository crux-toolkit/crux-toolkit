/**
 * \file SpectrumCollectionFactory.cpp 
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 14 June 2011
 * \brief Return a SpectrumCollection object of the appropriate
 * derived class.
 */
#include "parameter.h"
#include "SpectrumCollectionFactory.h"
#include "MS2SpectrumCollection.h"
#include "MSToolkitSpectrumCollection.h"
#include "MGFSpectrumCollection.h"
/**
 * Instantiates a SpectrumCollection based on the extension of the
 * given file and the use-mstoolkit and msgf options.
 */
SpectrumCollection* SpectrumCollectionFactory::create(const char* filename){
  SpectrumCollection* collection = NULL;

  if( has_extension(filename, ".mgf") || get_boolean_parameter("use-mgf") ){
    collection = new MGFSpectrumCollection(filename);
  } else if( has_extension(filename, ".mzXML") 
             || get_boolean_parameter("use-mstoolkit") ) {
    collection = new MSToolkitSpectrumCollection(filename);
  } else {
    collection = new MS2SpectrumCollection(filename);
  }
  return collection;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
