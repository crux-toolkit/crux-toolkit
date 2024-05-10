cmake_minimum_required(VERSION 3.15)
cmake_policy(VERSION 3.15)

###### Set up include and link directories  ######
include_directories(
  ${CMAKE_SOURCE_DIR}/LibraPeptideParser
  ${CMAKE_SOURCE_DIR}/LibraProteinParser
  ${CMAKE_BINARY_DIR}/ext/include
)
if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-pie")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
elseif (WIN32 AND NOT CYGWIN)
  # On windows set static linking of runtime libraries
  set(
    CMAKE_MSVC_RUNTIME_LIBRARY 
    "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL"
  )
endif(UNIX)

link_directories(
  ${CMAKE_BINARY_DIR}/ext/lib
)

ExternalProject_Add(
  libarchive-3.1.2
  URL file://${CMAKE_CURRENT_SOURCE_DIR}/extern/libarchive-3.1.2.tar.gz
  DOWNLOAD_NAME libarchive-3.1.2.tar.gz
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
  CMAKE_ARGS ${TOOLCHAIN_FILE}
    -DCMAKE_VERBOSE_MAKEFILE=ON
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_BUILD_TYPE=Release
  PATCH_COMMAND ""
)

set(
  libra_lib_files
  Common/crypt.cpp
  Common/get_random_number.cpp
  Common/hooks_tpp.cpp
  Common/kwset.cpp
  Common/obstack.cpp
  Common/spectStrct.cpp
  Common/tpp_tarball.cpp
  Common/TPPVersionInfo.cpp
  Common/util.cpp
  Common/wglob.cpp
  Common/Enzyme/EnzymeDigestion/AspNEnzymeDigestion.cpp
  Common/Enzyme/EnzymeDigestion/ChymotrypticEnzymeDigestion.cpp
  Common/Enzyme/EnzymeDigestion/CNBrEnzymeDigestion.cpp
  Common/Enzyme/EnzymeDigestion/Elastase.cpp
  Common/Enzyme/EnzymeDigestion/EnzymeDigestion.cpp
  Common/Enzyme/EnzymeDigestion/EnzymeSpecificity.cpp
  Common/Enzyme/EnzymeDigestion/GluC_bicarbEnzymeDigestion.cpp
  Common/Enzyme/EnzymeDigestion/GluCEnzymeDigestion.cpp
  Common/Enzyme/EnzymeDigestion/MultipleEnzymeDigestion.cpp
  Common/Enzyme/EnzymeDigestion/NonspecificEnzymeDigestion.cpp
  Common/Enzyme/EnzymeDigestion/TrypChymAspnEnzymeDigestion.cpp
  Common/Enzyme/EnzymeDigestion/TrypticCNBrEnzymeDigestion.cpp
  Common/Enzyme/EnzymeDigestion/TrypticEnzymeDigestion.cpp
  Common/Enzyme/EnzymeDigestionParser/EnzymeDigestionParser.cpp
  Common/Enzyme/EnzymeDigestionParser/EnzymeDigestionParserMain.cpp
  Common/Enzyme/ProteolyticEnzyme/ProteolyticEnzyme.cpp
  Common/Enzyme/ProteolyticEnzyme/ProteolyticEnzymeMain.cpp
  Common/Enzyme/ProteolyticEnzyme/ProteolyticEnzymeFactory/ProteolyticEnzymeFactory.cpp
  Common/ModificationInfo/ModificationInfo.cpp
  Common/ModifiedResidueMass/ModifiedResidueMass.cpp
  Common/ResidueMass/ResidueMass.cpp
  LibraPeptideParser/Intensity.cpp
  LibraPeptideParser/LibraConditionHandler.cpp
  LibraPeptideParser/LibraPeptideParser.cpp
  LibraPeptideParser/LibraPeptideParserMain.cpp
  LibraPeptideParser/LibraResult.cpp
  LibraPeptideParser/LibraSummary.cpp
  LibraPeptideParser/LibraWrapper.cpp
  LibraPeptideParser/Quantitation.cpp
  LibraProteinRatioParser/LibraConditionHandler2.cpp
  LibraProteinRatioParser/LibraGroupPeptideParser.cpp
  LibraProteinRatioParser/LibraProteinRatioParser.cpp
  LibraProteinRatioParser/LibraProteinRatioParserMain.cpp
  LibraProteinRatioParser/QuantitationFile.cpp
  LibraProteinRatioParser/StringConvertor.cpp
)
add_library(libra SHARED ${libra_lib_files})

if (WIN32 AND NOT CYGWIN)
  set_property(
    TARGET tide-support 
    PROPERTY 
      COMPILE_DEFINITIONS 
      GFLAGS_DLL_DECL= 
      GFLAGS_DLL_DECLARE_FLAG= 
      GFLAGS_DLL_DEFINE_FLAG=
  )
elseif (APPLE)
  set_property(
    TARGET tide-support 
    PROPERTY 
      COMPILE_DEFINITIONS 
      DARWIN 
  )
endif (WIN32 AND NOT CYGWIN)
