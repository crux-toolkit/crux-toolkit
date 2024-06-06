#ifndef NEOPEPXMLSTRUCTS_H
#define NEOPEPXMLSTRUCTS_H

#include <string>
#include <stdio.h>
#include <string.h>
#include <iostream>

enum mzMLElement :int{
  mzActivation,
  mzAnalyzer,
  mzBinary,
  mzBinaryDataArray,
  mzBinaryDataArrayList,
  mzChromatogramList,
  mzComponentList,
  mzCv,
  mzCvList,
  mzCvParam,
  mzDataProcessing,
  mzDataProcessingList,
  mzDetector,
  mzFileChecksum,
  mzFileContent,
  mzFileDescription,
  mzIndex,
  mzIndexedmzML,
  mzIndexList,
  mzIndexListOffset,
  mzInstrumentConfiguration,
  mzInstrumentConfigurationList,
  mzIsolationWindow,
  mzMzML,
  mzOffset,
  mzPrecursor,
  mzPrecursorList,
  mzProcessingMethod,
  mzReferenceableParamGroup,
  mzReferenceableParamGroupList,
  mzReferenceableParamGroupRef,
  mzRun,
  mzSample,
  mzSampleList,
  mzScan,
  mzScanList,
  mzScanWindow,
  mzScanWindowList,
  mzSelectedIon,
  mzSelectedIonList,
  mzSoftware,
  mzSoftwareList,
  mzSoftwareRef,
  mzSource,
  mzSourceFile,
  mzSourceFileList,
  mzSpectrum,
  mzSpectrumList,
  mzUserParam,
  MZML_NUM_ELEMENTS
};

enum mzMLWriteState{
  mzWS_Spectrum,
  mzWS_Chromat,
  mzWS_Done
};

static void NMZprintTabs(FILE* f, int tabs){
  for (int i = 0; i<tabs; i++) fprintf(f, " ");
}

static void NMZerrMsg(std::string el, std::string attr){
  std::cerr << el << "::" << attr << " required." << std::endl;
  exit(69);
}

typedef struct nmzDate {
  int year;
  int month;
  int day;
  nmzDate() {
    year = 0;
    month = 0;
    day = 0;
  }
} nmzDate;

typedef struct nmzTime {
  int hour;
  int minute;
  int second;
  nmzTime() {
    hour = 0;
    minute = 0;
    second = 0;
  }
} nmzTime;

typedef struct nmzDateTime{
  nmzDate date;
  nmzTime time;
  void clear(){
    date.day = 0;
    date.month = 0;
    date.year = 0;
    time.hour = 0;
    time.minute = 0;
    time.second = 0;
  }
  void parseDateTime(const char* dt){
    if (strlen(dt)<2){
      clear();
      return;
    }
    int x = sscanf(dt, "%d-%d-%dT%d:%d:%d", &date.year, &date.month, &date.day, &time.hour, &time.minute, &time.second);
  }
  void parseDateTime(std::string s){
    parseDateTime(s.c_str());
  }
  std::string write(){
    std::string s;
    char str[64];
    sprintf(str, "%4d-%02d-%02dT%02d:%02d:%02d", date.year, date.month, date.day, time.hour, time.minute, time.second);
    s = str;
    return s;
  }
} nmzDateTime;

#ifdef _MSC_VER
//#define __inline__ __inline
//typedef _int64  __int64_t;
typedef unsigned _int32 uint32_t;
typedef unsigned _int64 uint64_t;
typedef __int64 f_off_nmz;
#define nmzfseek(h,p,o) _fseeki64(h,p,o)
#define nmzftell(h) _ftelli64(h)
#define nmzatoi64(h) _atoi64(h)
#endif

//For MinGW toolset, which lacks the ftello, fseeko, etc functions
#ifdef MINGW
#define nmzfseek(h,p,o) fseeko64(h,p,o)
#define nmzftell(h) ftello64(h)
#define nmzatoi64(h) _atoi64(h)
#endif

#if defined(GCC) || defined(__LINUX__) || defined(GNUC) || defined(__MINGW32__)
#include <stdint.h>
#include <stdexcept>
#ifndef _LARGEFILE_SOURCE
#error "need to define _LARGEFILE_SOURCE!!"
#endif    /* end _LARGEFILE_SOURCE */
#if _FILE_OFFSET_BITS<64
#error "need to define _FILE_OFFSET_BITS=64"
#endif

typedef off_t f_off_nmz;
#define nmzfseek(h,p,o) fseeko(h,p,o)
#define nmzftell(h) ftello(h)
#define nmzatoi64(h) atoll(h)
#endif

#endif
