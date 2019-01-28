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
SET(USE_LATEST_PWIZ ON CACHE BOOL "Use lastest version of ProteoWizard")

if (UNIX AND BUILD_32)
  SET(TOOLCHAIN_FILE "-DCMAKE_TOOLCHAIN_FILE=${CMAKE_SOURCE_DIR}/linux_i686.toolchain.cmake")
else (UNIX AND BUILD_32)
  SET(TOOLCHAIN "")
endif (UNIX AND BUILD_32)

if (WIN32 AND NOT CYGWIN)
  # If windows build we need to select 32/64 bit platform
  if (CMAKE_GENERATOR MATCHES "^.*Win64$")
    SET(MSBUILD_PLATFORM "x64")
  else (CMAKE_GENERATOR MATCHES "^.*Win64$")
    SET(MSBUILD_PLATFORM "Win32")
  endif (CMAKE_GENERATOR MATCHES "^.*Win64$")

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
      URL https://noble.gs.washington.edu/crux-downloads/pwiz-src-3_0_18137_c65b6f6.tar.bz2
  )
endif (USE_LATEST_PWIZ)

# Download and build a partial distribution of the pwiz source
if (WIN32 AND NOT CYGWIN)
    ExternalProject_Add(
      ProteoWizard
      ${PWIZ_SOURCE}
      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
      CONFIGURE_COMMAND ""
      BUILD_COMMAND
        ${CMAKE_COMMAND}
          -DMSBUILD_PLATFORM=${MSBUILD_PLATFORM}
          -DBUILD_TYPE:STRING=${CMAKE_CFG_INTDIR}
          -DPREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
          -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard
          -DINCLUDE_VENDOR_LIBRARIES:BOOLEAN=${INCLUDE_VENDOR_LIBRARIES}
          -P  ${CMAKE_CURRENT_SOURCE_DIR}/build_pwiz.cmake
      INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/patches/boost ${CMAKE_CURRENT_BINARY_DIR}/include/boost
      BUILD_IN_SOURCE 1
    )
  # On Windows we have to do some extra house keeping
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
    PATCH_COMMAND
      ${CMAKE_COMMAND}  -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/patches/pwiz/tar.jam
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/libraries &&
      ${CMAKE_COMMAND}  -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/patches/pwiz/WhittakerSmoother.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/pwiz/analysis/common
    CONFIGURE_COMMAND ""
    BUILD_COMMAND
      ${CMAKE_COMMAND}
        -DBUILD_TYPE:STRING=${CMAKE_CFG_INTDIR}
        -DPREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
        -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard
        -P  ${CMAKE_CURRENT_SOURCE_DIR}/build_pwiz.cmake
    INSTALL_COMMAND
      ${CMAKE_CURRENT_SOURCE_DIR}/create_links.sh
         --prefix=${CMAKE_CURRENT_BINARY_DIR}
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
    CONFIGURE_COMMAND ""
    BUILD_COMMAND
      ${CMAKE_COMMAND}
        -DBUILD_32:BOOL=${BUILD_32}
        -DBUILD_TYPE:STRING=${CMAKE_CFG_INTDIR}
        -DPREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
        -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard
        -P  ${CMAKE_CURRENT_SOURCE_DIR}/build_pwiz.cmake
    INSTALL_COMMAND
      ${CMAKE_CURRENT_SOURCE_DIR}/create_links.sh
         --prefix=${CMAKE_CURRENT_BINARY_DIR}
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

if (UNIX)
  ExternalProject_Add(
    percolator
    GIT_REPOSITORY https://github.com/percolator/percolator.git
    GIT_TAG ${PERCOLATOR_GIT_VERSION}
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    CMAKE_ARGS ${TOOLCHAIN_FILE}
      -DCMAKE_VERBOSE_MAKEFILE=ON
      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_BUILD_TYPE=${CMAKE_CFG_INTDIR}
      -DCRUX=${CMAKE_CURRENT_SOURCE_DIR}/..
      -DEXT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/include
    PATCH_COMMAND ${CMAKE_COMMAND}  -E copy
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/src/CMakeLists.txt
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/CMakeLists.txt
    INSTALL_COMMAND ""
  )
else (UNIX)
  ExternalProject_Add(
    percolator
    GIT_REPOSITORY https://github.com/percolator/percolator.git
    GIT_TAG ${PERCOLATOR_GIT_VERSION}
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    CMAKE_ARGS -DCMAKE_VERBOSE_MAKEFILE=ON
      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_BUILD_TYPE=${CMAKE_CFG_INTDIR}
      -DCRUX=${CMAKE_CURRENT_SOURCE_DIR}/..
      -DEXT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/include
    INSTALL_COMMAND ""
  )
endif (UNIX)

if (UNIX)
  ExternalProject_Add(
    protobuf-3.0.0
    URL https://github.com/protocolbuffers/protobuf/releases/download/v3.0.0/protobuf-cpp-3.0.0.tar.gz
    DOWNLOAD_NAME protobuf-cpp-3.0.0.tar.gz
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    SOURCE_SUBDIR cmake
    CMAKE_ARGS ${TOOLCHAIN_FILE}
      -DCMAKE_VERBOSE_MAKEFILE=ON
      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_BUILD_TYPE=${CMAKE_CFG_INTDIR}
      -Dprotobuf_BUILD_TESTS=OFF
 )
elseif (WIN32 AND NOT CYGWIN)
  ExternalProject_Add(
    protobuf-3.0.0
    URL https://github.com/protocolbuffers/protobuf/releases/download/v3.0.0/protobuf-cpp-3.0.0.tar.gz
    DOWNLOAD_NAME protobuf-cpp-3.0.0.tar.gz
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    SOURCE_SUBDIR cmake
    CMAKE_ARGS
      -DCMAKE_VERBOSE_MAKEFILE=ON
      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_BUILD_TYPE=${CMAKE_CFG_INTDIR}
      -Dprotobuf_MSVC_STATIC_RUNTIME=OFF
    PATCH_COMMAND ""
  )
endif (UNIX)

ExternalProject_Add(
  gflags-2.2.1
  URL https://codeload.github.com/gflags/gflags/tar.gz/v2.2.1
  DOWNLOAD_NAME gflags-2.2.1.tar.gz
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
  CMAKE_ARGS ${TOOLCHAIN_FILE}
    -DCMAKE_VERBOSE_MAKEFILE=ON
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_BUILD_TYPE=${CMAKE_CFG_INTDIR}
  PATCH_COMMAND ""
)

if (UNIX)
  # Get Comet source from SVN and buld using it's own makefile
  set(comet_cxx_flags ${CMAKE_CXX_FLAGS}\ -DCRUX)
  set(mstoolkit_c_flags ${CMAKE_C_FLAGS}\ -DCRUX)
  ExternalProject_Add(
    comet
    SVN_REPOSITORY svn://svn.code.sf.net/p/comet-ms/code/tags/release_2018012/
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    PATCH_COMMAND
      cp ${CMAKE_CURRENT_SOURCE_DIR}/patches/comet/CometSearch/CometSearchManager.cpp
         ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/ &&
      cp ${CMAKE_CURRENT_SOURCE_DIR}/patches/comet/MSToolkit/Makefile
         ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/MSToolkit/ &&
      cp ${CMAKE_CURRENT_SOURCE_DIR}/patches/comet/MSToolkit/include/MSReader.h
         ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/MSToolkit/include/ &&
      cp ${CMAKE_CURRENT_SOURCE_DIR}/patches/comet/MSToolkit/src/MSToolkit/MSReader.cpp
         ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/MSToolkit/src/MSToolkit/
    BUILD_COMMAND 
      COMMAND make -C MSToolkit CXXFLAGS=${comet_cxx_flags} CFLAGS=${mstoolkit_c_flags} 
                                LDFLAGS=${CMAKE_LD_FLAGS}
      COMMAND make -C CometSearch CXXFLAGS=${comet_cxx_flags} LDFLAGS=${CMAKE_LD_FLAGS}
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
  )
  # Copy the Comet and MSToolkit include files and libraries to
  # the desired directories.
  add_custom_command(
    TARGET comet
    POST_BUILD
    COMMAND
      [ -d ${CMAKE_CURRENT_BINARY_DIR}/include/MSToolkit ] ||
        mkdir ${CMAKE_CURRENT_BINARY_DIR}/include/MSToolkit
    COMMAND
      cp -p
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/MSToolkit/include/*
      ${CMAKE_CURRENT_BINARY_DIR}/include/MSToolkit
    COMMAND
      cp -p
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/MSToolkit/libmstoolkitlite.a
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/libcometsearch.a
      ${CMAKE_CURRENT_BINARY_DIR}/lib
    COMMAND
      [ -d ${CMAKE_CURRENT_BINARY_DIR}/include/CometSearch ] ||
        mkdir ${CMAKE_CURRENT_BINARY_DIR}/include/CometSearch
    COMMAND
      cp -p
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/Common.h
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/CometData.h
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/CometSearch.h
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/CometSearchManager.h
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/CometPreprocess.h
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/ThreadPool.h
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/Threading.h
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/OSSpecificThreading.h
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/CometInterfaces.h
      ${CMAKE_CURRENT_BINARY_DIR}/include/CometSearch
  )
elseif (WIN32 AND NOT CYGWIN)
  set(LIBPATH "${CMAKE_CURRENT_BINARY_DIR}/lib;")
  ExternalProject_Add(
    comet
    SVN_REPOSITORY http://svn.code.sf.net/p/comet-ms/code/tags/release_2018012/
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    UPDATE_COMMAND ""
    # Copy the Comet and MSToolkit include files and libraries to
    # the desired directories.
    CONFIGURE_COMMAND
      ${CMAKE_COMMAND} -E copy_directory
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/MSToolkit/include
      ${CMAKE_CURRENT_BINARY_DIR}/include/MSToolkit/
    PATCH_COMMAND
      ${CMAKE_COMMAND}
        -DSOURCE_DIR:STRING=${CMAKE_CURRENT_SOURCE_DIR}
        -DBINARY_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/patch-comet.cmake
    BUILD_COMMAND
      ${CMAKE_COMMAND}
        -DLIBPATH:STRING="${LIBPATH}"
        -DMSBUILD_ARGS:STRING=/p:Configuration=${CMAKE_CFG_INTDIR}
        -DMSBUILD_PLATFORM=${MSBUILD_PLATFORM}
        -DPROJ_FILE:STRING=comet.sln
        -P ${CMAKE_CURRENT_SOURCE_DIR}/build-comet.cmake
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
  ) 
  ExternalProject_Add_Step(
    comet
    copy_comet_library
    COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/${MSBUILD_PLATFORM}/${CMAKE_CFG_INTDIR}/CometSearch.lib
      ${CMAKE_CURRENT_BINARY_DIR}/lib/CometSearch.lib
    COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/${MSBUILD_PLATFORM}/${CMAKE_CFG_INTDIR}/MSToolkitLite.lib
      ${CMAKE_CURRENT_BINARY_DIR}/lib/MSToolkitLite.lib
    DEPENDEES install
  )
endif (UNIX)


# On Windows we have to include the vendor libraries in the package
if (WIN32 AND NOT CYGWIN AND INCLUDE_VENDOR_LIBRARIES)
  set(
    vendor_binaries
    ${CMAKE_CURRENT_BINARY_DIR}/lib/agtsampleinforw.dll
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
    ${CMAKE_CURRENT_BINARY_DIR}/lib/concrt140.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/CRHAKEI2.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/DataReader.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/DualProbeInterfaceParametersCS.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Google.Protobuf.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/HSReadWrite.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/IdentityModel.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Interop.EDAL.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Interop.HSREADWRITELib.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/IOModuleQTFL.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/MassCalcWrapObject.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/MassSpecDataReader.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/mfc140.dll
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
    ${CMAKE_CURRENT_BINARY_DIR}/lib/pwiz_bindings_cli.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/pwiz_data_cli.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/QTFLDebugLog.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/QTFLProtoCS.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Shimadzu.LabSolutions.IO.IoModule.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/System.Data.SQLite.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/System.Runtime.Caching.Generic.dll
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
  else()
    set(
      # The 32-bit Windows buid needs these additional DLL
      vendor_binaries
      ${vendor_binaries}
      ${CMAKE_CURRENT_BINARY_DIR}/lib/Interop.DataExplorer.dll
      ${CMAKE_CURRENT_BINARY_DIR}/lib/MSFileReader.XRawfile2.dll
    )
  endif (CMAKE_GENERATOR MATCHES "^.*Win64$")
  install (FILES ${vendor_binaries} DESTINATION bin)
endif (WIN32 AND NOT CYGWIN AND INCLUDE_VENDOR_LIBRARIES)
