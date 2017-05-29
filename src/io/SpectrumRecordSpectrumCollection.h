#ifndef SPECTRUM_RECORD_SPECTRUM_COLLECTION_H
#define SPECTRUM_RECORD_SPECTRUM_COLLECTION_H

#include "SpectrumCollection.h"

class SpectrumRecordSpectrumCollection : public Crux::SpectrumCollection {
 public:
  SpectrumRecordSpectrumCollection(const std::string& filename);
  virtual ~SpectrumRecordSpectrumCollection();
  static bool IsSpectrumRecordFile(const std::string& filename);
  virtual bool parse();
  virtual Crux::Spectrum* getSpectrum(int first_scan);
  virtual bool getSpectrum(int first_scan, Crux::Spectrum* spectrum);
};

#endif

