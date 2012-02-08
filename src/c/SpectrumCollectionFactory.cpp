/**
 * \file SpectrumCollectionFactory.cpp 
 * AUTHOR Barbara Frewen
 * CREATE DATE: 14 June 2011
 * \brief Return a SpectrumCollection object of the appropriate
 * derived class.
 */
#include "parameter.h"
#include "SpectrumCollectionFactory.h"
#include "MS2SpectrumCollection.h"
#include "MSToolkitSpectrumCollection.h"
#include "MGFSpectrumCollection.h"
#include <sys/types.h>
#include <sys/stat.h>

/**
 * Instantiates a SpectrumCollection based on the extension of the
 * given file and the use-mstoolkit and msgf options.
 */
SpectrumCollection* SpectrumCollectionFactory::create(const char* filename){
  SpectrumCollection* collection = NULL;
  struct stat stat_buff ; 
  stat(filename, &stat_buff);
  if(stat(filename, &stat_buff)!=0){
    carp(CARP_FATAL, "The file %s does not exist. \n", filename);
  }
  if (S_ISDIR(stat_buff.st_mode)){
    carp(CARP_FATAL, "Path %s is a directory. \n Please enter a spectrum filename\
(.ms2,.mgf, or .mzXML)",filename);
  }else if( has_extension(filename, ".mgf") || get_boolean_parameter("use-mgf") ){
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
