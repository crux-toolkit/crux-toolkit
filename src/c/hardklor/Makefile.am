# Automake template for xlink-search library


AM_CXXFLAGS = -Wall \
              -W \
              $(USEDOUBLES) \
              $(USEDARWIN) \
              -D_LARGEFILE_SOURCE=1 \
              -D_FILE_OFFSET_BITS=64 \
              -D_GNU_SOURCE


libhardklor_la_SOURCES = \
  CAveragine.cpp \
  CAveragine.h \
  CHardklor.h \
  CHardklor.cpp \
  CHardklor2.h \
  CHardklor2.cpp \
  CHardklorFileReader.h \
  CHardklorFileReader.cpp \
  CHardklorParser.h \
  CHardklorParser.cpp \
  CHardklorProtein.h \
  CHardklorProtein.cpp \
  CHardklorSetting.h \
  CHardklorSetting.cpp \
  CHardklorVariant.h \
  CHardklorVariant.cpp \
  CMercury8.h \
  CMercury8.cpp \
  CModelLibrary.h \
  CModelLibrary.cpp \
  CNoiseReduction.h \
  CNoiseReduction.cpp \
  CPeriodicTable.h \
  CPeriodicTable.cpp \
  CSettings.h \
  CSettings.cpp \
  CSpecAnalyze.h \
  CSpecAnalyze.cpp \
  CSplitSpectrum.h \
  CSplitSpectrum.cpp \
  FFT.h \
  FFT.cpp \
  FFT-HK.h \
  FFT-HK.cpp \
  HardklorApp.cpp \
  HardklorTypes.h \
  mercury.h \
  resource1.h \
  resource.h \
  S2N.cpp \
  S2N.h \
  Smooth.h \
  Smooth.cpp \
  SpecAnalyzeSupport.cpp \
  SpecAnalyzeSupport.h \
  CruxHardklorApplication.cpp \
  CruxHardklorApplication.h

libhardklor_la_CXXFLAGS = $(AM_CXXFLAGS) \
  -I$(top_srcdir)/src/c/ \
  -I$(top_srcdir)/src/external/MSToolkit/include \
  -DGCC \
  -D_LARGEFILE_SOURCE \
  -DHAVE_EXPAT_CONFIG_H \
  -D_FILE_OFFSET_BITS=64 \
  -DCRUX


noinst_LTLIBRARIES = libhardklor.la

