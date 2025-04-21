include(ExternalProject)

INCLUDE_DIRECTORIES(BEFORE SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/include)

if (UNIX)
  SET(ENV{BOOSTROOT} ${CMAKE_CURRENT_BINARY_DIR})
  SET(ENV{LIBPATH} ${CMAKE_CURRENT_BINARY_DIR}/lib:$ENV{LIBPATH})
  SET(ENV{LD_LIBRARY_PATH} ${CMAKE_CURRENT_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH})
  SET(ENV{PATH} ${CMAKE_CURRENT_BINARY_DIR}/bin:$ENV{PATH})
  SET(CPP_ARGS "-I${CMAKE_CURRENT_BINARY_DIR}/include")
  SET(LD_ARGS "-L${CMAKE_CURRENT_BINARY_DIR}/lib")
endif (UNIX)
SET(ENV{CMAKE_PREFIX_PATH} ${CMAKE_CURRENT_BINARY_DIR})
SET(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR})
SET(USE_LATEST_PWIZ OFF CACHE BOOL "Use lastest version of ProteoWizard")

if (WIN32 AND NOT CYGWIN)
  SET(MSBUILD_PLATFORM "x64")
  # On Windows we should build vendor libraries by default
  if (NOT DEFINED INCLUDE_VENDOR_LIBRARIES)
    set(INCLUDE_VENDOR_LIBRARIES TRUE)
  endif (NOT DEFINED INCLUDE_VENDOR_LIBRARIES)
  if (INCLUDE_VENDOR_LIBRARIES)
    message(STATUS "Including vendor libraries.")
  else (INCLUDE_VENDOR_LIBRARIES)
    message(STATUS "Not including vendor libraries.")
  endif (INCLUDE_VENDOR_LIBRARIES)
endif (WIN32 AND NOT CYGWIN)

# Set up which source to use for ProteoWizard
if (USE_LATEST_PWIZ)
  # Use the ProteoWizard TeamCity repository.
  set(
      PWIZ_SOURCE
      DOWNLOAD_COMMAND
        ${CMAKE_COMMAND}
          -DDOWNLOAD_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/downloads
          -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard
          -P  ${CMAKE_CURRENT_SOURCE_DIR}/download_pwiz.cmake
      DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/downloads
  )
else (USE_LATEST_PWIZ)
  # Use tarball cached on our server
  set(
      PWIZ_SOURCE
      URL https://noble.gs.washington.edu/crux-downloads/pwiz-src-3_0_24044_fd6604f.tar.bz2
  )
endif (USE_LATEST_PWIZ)

# MSToolkit
# Get MSToolkit source from Github
ExternalProject_Add(
  MSToolkit
  GIT_REPOSITORY https://github.com/mhoopmann/mstoolkit.git
  GIT_TAG 20e99ce
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
  CMAKE_ARGS ${TOOLCHAIN_FILE}
    -DCMAKE_VERBOSE_MAKEFILE=ON
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_BUILD_TYPE=${CMAKE_CFG_INTDIR}
  PATCH_COMMAND
    COMMAND
      unzip -o src/zlib1211.zip -d src
    COMMAND
      tar xvf src/expat-2.2.9.tar.gz -C src
    COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/MSToolkit/CMakeLists.txt
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/MSToolkit/
    COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/MSToolkit/src/RAMPface.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/MSToolkit/src/mzParser/RAMPface.cpp
  BUILD_IN_SOURCE 1
)

# Download and build a partial distribution of the pwiz source
if (WIN32 AND NOT CYGWIN)
    ExternalProject_Add(
      ProteoWizard
      ${PWIZ_SOURCE}
      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
      CONFIGURE_COMMAND
        # We don't need to build the .Net/SWIG or Java bindings
        COMMAND ${CMAKE_COMMAND} -E rm -r ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/pwiz/utility/bindings/SWIG
        COMMAND ${CMAKE_COMMAND} -E rm -r ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/pwiz/utility/bindings/CLI
          ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/pwiz/utility/bindings
        # Use patched Jamroot.jam to turn off unit testing
        COMMAND ${CMAKE_COMMAND} -E copy
          ${CMAKE_CURRENT_SOURCE_DIR}/patches/pwiz/Jamroot.jam
          ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/
      BUILD_COMMAND
        ${CMAKE_COMMAND}
          -DMSBUILD_PLATFORM=${MSBUILD_PLATFORM}
          -DBUILD_TYPE:STRING=${CMAKE_CFG_INTDIR}
          -DPREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
          -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard
          -DINCLUDE_VENDOR_LIBRARIES:BOOLEAN=${INCLUDE_VENDOR_LIBRARIES}
          -P  ${CMAKE_CURRENT_SOURCE_DIR}/build_pwiz.cmake
      INSTALL_COMMAND
        COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/libraries/boost_1_76_0/boost/
          ${CMAKE_CURRENT_BINARY_DIR}/include/boost
        COMMAND ${CMAKE_COMMAND} -E copy
          ${CMAKE_CURRENT_SOURCE_DIR}/patches/boost/asio.hpp
          ${CMAKE_CURRENT_BINARY_DIR}/include/boost
        COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${CMAKE_CURRENT_SOURCE_DIR}/patches/boost/asio
          ${CMAKE_CURRENT_BINARY_DIR}/include/boost/asio
        COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${CMAKE_CURRENT_SOURCE_DIR}/patches/boost/wintls
          ${CMAKE_CURRENT_BINARY_DIR}/include/boost/wintls
        COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/libraries/boost_aux/boost/utility
          ${CMAKE_CURRENT_BINARY_DIR}/include/boost/utility
        COMMAND ${CMAKE_COMMAND} -E copy
          ${CMAKE_CURRENT_SOURCE_DIR}/patches/boost/cmath.hpp
          ${CMAKE_CURRENT_BINARY_DIR}/include/boost/core
        COMMAND ${CMAKE_COMMAND} -E copy
          ${CMAKE_CURRENT_SOURCE_DIR}/patches/boost/wintls.hpp
          ${CMAKE_CURRENT_BINARY_DIR}/include/boost/
      BUILD_IN_SOURCE 1
    )
  # On Windows we have to do some extra house keeping
  # To install the libraries
  ExternalProject_Add_Step(
    ProteoWizard
    windows_post_install
    COMMAND
      ${CMAKE_COMMAND}
        -DBUILD_TYPE:STRING=${CMAKE_CFG_INTDIR}
        -DPREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
        -DMSBUILD_PLATFORM:STRING=${MSBUILD_PLATFORM}
        -P  ${CMAKE_CURRENT_SOURCE_DIR}/post_install_pwiz.cmake
    DEPENDEES install
  )
elseif (APPLE)
  ExternalProject_Add(
    ProteoWizard
    ${PWIZ_SOURCE}
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    CONFIGURE_COMMAND
      COMMAND ${CMAKE_COMMAND} -E rm -r ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/pwiz/utility/bindings
      # Use patched Jamroot.jam to turn off unit testing
      COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/patches/pwiz/Jamroot.jam
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/
    BUILD_COMMAND
      ${CMAKE_COMMAND}
        -DBUILD_TYPE:STRING=${CMAKE_CFG_INTDIR}
        -DPREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
        -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard
        -P  ${CMAKE_CURRENT_SOURCE_DIR}/build_pwiz.cmake
    INSTALL_COMMAND
      ${CMAKE_CURRENT_SOURCE_DIR}/create_links.sh
         --prefix=${CMAKE_CURRENT_BINARY_DIR}
      COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/libraries/boost_1_76_0/boost/
        ${CMAKE_CURRENT_BINARY_DIR}/include/boost
      COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/patches/boost/asio.hpp
        ${CMAKE_CURRENT_BINARY_DIR}/include/boost
      COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/patches/boost/asio
        ${CMAKE_CURRENT_BINARY_DIR}/include/boost/asio
      COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/libraries/boost_aux/boost/utility
        ${CMAKE_CURRENT_BINARY_DIR}/include/boost/utility
    BUILD_IN_SOURCE 1
  )
else (UNIX)
  ExternalProject_Add(
    ProteoWizard
    ${PWIZ_SOURCE}
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    CONFIGURE_COMMAND 
      COMMAND ${CMAKE_COMMAND} -E rm -r 
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/pwiz/utility/bindings
      # Use patched Jamroot.jam to turn off unit testing
      COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/patches/pwiz/Jamroot.jam
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/
    BUILD_COMMAND
      ${CMAKE_COMMAND}
        -DBUILD_TYPE:STRING=${CMAKE_CFG_INTDIR}
        -DPREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
        -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard
        -P ${CMAKE_CURRENT_SOURCE_DIR}/build_pwiz.cmake ${CMAKE_CURRENT_BINARY_DIR}
    INSTALL_COMMAND
      ${CMAKE_CURRENT_SOURCE_DIR}/create_links.sh
         --prefix=${CMAKE_CURRENT_BINARY_DIR}
      COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/libraries/boost_1_76_0/boost/
        ${CMAKE_CURRENT_BINARY_DIR}/include/boost
      COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/patches/boost/asio.hpp
        ${CMAKE_CURRENT_BINARY_DIR}/include/boost
      COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/patches/boost/asio
        ${CMAKE_CURRENT_BINARY_DIR}/include/boost/asio
      COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/libraries/boost_aux/boost/utility
        ${CMAKE_CURRENT_BINARY_DIR}/include/boost/utility
    BUILD_IN_SOURCE 1
  )
endif (WIN32 AND NOT CYGWIN)

set(PERCOLATOR_GIT_VERSION 310f924)
ExternalProject_Add(
  percolator
  GIT_REPOSITORY https://github.com/percolator/percolator.git
  GIT_TAG ${PERCOLATOR_GIT_VERSION}
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
  CMAKE_ARGS ${TOOLCHAIN_FILE}
    -DCMAKE_VERBOSE_MAKEFILE=ON
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCRUX=${CMAKE_CURRENT_SOURCE_DIR}/..
    -DEXT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/include
    -DWITHOUT_GTEST=ON
  PATCH_COMMAND
    COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/CMakeLists.txt
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src
    COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/picked_protein/CMakeLists.txt
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/picked_protein
  INSTALL_COMMAND ""
)

if (UNIX)
  ExternalProject_Add(
    protobuf-3.19.4
    URL https://github.com/protocolbuffers/protobuf/releases/download/v3.19.4/protobuf-all-3.19.4.tar.gz
    DOWNLOAD_NAME protobuf-cpp-3.19.4.tar.gz
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    SOURCE_SUBDIR cmake
    CMAKE_ARGS ${TOOLCHAIN_FILE}
      -DCMAKE_VERBOSE_MAKEFILE=ON
      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_BUILD_TYPE=${CMAKE_CFG_INTDIR}
      -Dprotobuf_BUILD_TESTS=OFF
    PATCH_COMMAND
      # Patch to compile with c++14
      COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/patches/protobuf-3.19.4/CMakeLists.txt
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/protobuf-3.19.4/cmake
 )
elseif (WIN32 AND NOT CYGWIN)
  ExternalProject_Add(
    protobuf-3.19.4
    URL https://github.com/protocolbuffers/protobuf/releases/download/v3.19.4/protobuf-all-3.19.4.tar.gz
    DOWNLOAD_NAME protobuf-cpp-3.19.4.tar.gz
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    SOURCE_SUBDIR cmake
    CMAKE_ARGS
      -DCMAKE_VERBOSE_MAKEFILE=ON
      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_BUILD_TYPE=${CMAKE_CFG_INTDIR}
      -Dprotobuf_BUILD_TESTS=OFF
      -Dprotobuf_WITH_ZLIB=OFF
      -Dprotobuf_MSVC_STATIC_RUNTIME=OFF
    PATCH_COMMAND
      # There is a bug in the 3.19.4 release of protobuf that prevents it
      # from honoring the protobuf_MSVC_STATIC_RUNTIME=OFF option,
      # So for now we have to patch the CMakeLists.txt file.
      # This should be fixed in the 3.20.0 release.
      COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/patches/protobuf-3.19.4/CMakeLists.txt
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/protobuf-3.19.4/cmake
  )
endif (UNIX)

ExternalProject_Add(
  gflags-2.2.2
  URL https://codeload.github.com/gflags/gflags/tar.gz/v2.2.2
  DOWNLOAD_NAME gflags-2.2.2.tar.gz
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
  CMAKE_ARGS ${TOOLCHAIN_FILE}
    -DCMAKE_VERBOSE_MAKEFILE=ON
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_BUILD_TYPE=${CMAKE_CFG_INTDIR}
  PATCH_COMMAND ""
)

# Comet
set(comet_cxx_flags ${CMAKE_CXX_FLAGS}\ -DCRUX)
set(base_include_path ${CMAKE_BINARY_DIR}/ext/include)
set(mstoolkit_include_path ${CMAKE_BINARY_DIR}/ext/include/MSToolkit)
if (WIN32 AND NOT CYGWIN)
  set(LIBPATH "${CMAKE_CURRENT_BINARY_DIR}/lib;")
endif (WIN32 AND NOT CYGWIN)
ExternalProject_Add(
  comet
  GIT_REPOSITORY https://github.com/UWPR/Comet/
  GIT_TAG v2025.01.1
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
  CMAKE_ARGS ${TOOLCHAIN_FILE}
    -DBASE_INCLUDE_PATH=${base_include_path}
    -DMSTOOLKIT_INCLUDE_PATH=${mstoolkit_include_path}
    -DCMAKE_VERBOSE_MAKEFILE=ON
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_BUILD_TYPE=${CMAKE_CFG_INTDIR}
    PATCH_COMMAND
      COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/patches/comet/CMakeLists.txt
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/
    BUILD_IN_SOURCE 1
)

# NeoPepXMLParser 
set(base_include_path ${CMAKE_BINARY_DIR}/ext/include)
set(mstoolkit_include_path ${CMAKE_BINARY_DIR}/ext/include/MSToolkit)
ExternalProject_Add(
  neopepxmlparser
  GIT_TAG e89f166
  GIT_REPOSITORY https://github.com/mhoopmann/NeoPepXMLParser.git
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
  CMAKE_ARGS ${TOOLCHAIN_FILE}
    -DBASE_INCLUDE_PATH=${base_include_path}
    -DMSTOOLKIT_INCLUDE_PATH=${mstoolkit_include_path}
    -DCMAKE_VERBOSE_MAKEFILE=ON
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
PATCH_COMMAND
  COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/NeoPepXMLParser/CMakeLists.txt
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/neopepxmlparser/
  BUILD_IN_SOURCE 1
)

# Kojak
set(base_include_path ${CMAKE_BINARY_DIR}/ext/include)
set(mstoolkit_include_path ${CMAKE_BINARY_DIR}/ext/include/MSToolkit)
set(neopepxmlparser_include_path ${CMAKE_BINARY_DIR}/ext/build/src/neopepxmlparser)
set(hardklor_include_path ${CMAKE_SOURCE_DIR}/src/app/hardklor)
ExternalProject_Add(
  kojak
  GIT_TAG 4307dfa
  GIT_REPOSITORY https://github.com/mhoopmann/kojak.git
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
  CMAKE_ARGS ${TOOLCHAIN_FILE}
    -DBASE_INCLUDE_PATH=${base_include_path}
    -DMSTOOLKIT_INCLUDE_PATH=${mstoolkit_include_path}
    -DHARDKLOR_INCLUDE_PATH=${hardklor_include_path}
    -DNEOPEPXMLPARSER_INCLUDE_PATH=${neopepxmlparser_include_path}
    -DCMAKE_VERBOSE_MAKEFILE=ON
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
PATCH_COMMAND
  COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/kojak/CMakeLists.txt
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/kojak/
  COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/kojak/KLog.h
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/kojak/
  COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/kojak/KParams.h
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/kojak/
  BUILD_IN_SOURCE 1
)

# On Windows we have to include the vendor libraries in the package
if (WIN32 AND NOT CYGWIN AND INCLUDE_VENDOR_LIBRARIES)
  set(
    vendor_binaries
    ${CMAKE_CURRENT_BINARY_DIR}/lib/agtsampleinforw.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/baf2sql_c.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BaseCommon.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BaseDataAccess.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BaseError.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BaseTof.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BDal.BCO.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BDal.BCO.Interfaces.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BDal.CCO.Calibration.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BDal.CCO.Interfaces.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BDal.CCO.TemperatureCompensation.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BDal.CCO.Transformation.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BDal.CCO.Utilities.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BDal.CXt.Lc.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BDal.CXt.Lc.Factory.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BDal.CXt.Lc.Interfaces.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/BDal.CXt.Lc.UntU2.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/CLFIO32.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/CompassXtractMS.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/CRHAKEI2.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/DualProbeInterfaceParametersCS.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Google.Protobuf.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/HSReadWrite.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/IdentityModel.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Interop.EDAL.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Interop.HSREADWRITELib.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/IOModuleQTFL.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/MassCalcWrapObject.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/MassSpecDataReader.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/mfc90.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/MIDAC.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/mkl_sequential.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/MSMSDBCntl.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/msvcm90.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/msvcp110.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/msvcp120.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/msvcp140.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/msvcp90.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/msvcr110.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/msvcr120.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/msvcr90.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Newtonsoft.Json.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/ParallelExtensionsExtras.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/PeakItgLSS.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/QTFLDebugLog.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/QTFLProtoCS.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Shimadzu.LabSolutions.IO.IoModule.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/System.Data.SQLite.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/System.Runtime.Caching.Generic.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/ThermoFisher.CommonCore.Data.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/ThermoFisher.CommonCore.RawFileReader.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/timsdata.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/ucrtbase.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/UIMFLibrary.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/unifi-protobuf-net.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/vcomp110.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/vcomp140.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/vcomp90.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/vcruntime140.dll
  )
  if (CMAKE_GENERATOR MATCHES "^.*Win64$")
    # The 64-bit Windows buid needs these additional DLL
    set(
      vendor_binaries
      ${vendor_binaries}
      ${CMAKE_CURRENT_BINARY_DIR}/lib/ThermoFisher.CommonCore.Data.dll
      ${CMAKE_CURRENT_BINARY_DIR}/lib/ThermoFisher.CommonCore.RawFileReader.dll
    )
  endif (CMAKE_GENERATOR MATCHES "^.*Win64$")
  install (FILES ${vendor_binaries} DESTINATION bin)
endif (WIN32 AND NOT CYGWIN AND INCLUDE_VENDOR_LIBRARIES)
