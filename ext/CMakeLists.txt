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

# Download and build a partial distribution of the pwiz source
# from the ProteoWizard TeamCity repository.
if (WIN32 AND NOT CYGWIN)
  if (USE_LATEST_PWIZ)
    ExternalProject_Add(
      ProteoWizard
      DOWNLOAD_COMMAND 
        ${CMAKE_COMMAND} 
          -DDOWNLOAD_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/downloads
          -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard
          -P  ${CMAKE_CURRENT_SOURCE_DIR}/download_pwiz.cmake
      DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/downloads
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
  else (USE_LATEST_PWIZ)
    ExternalProject_Add(
      ProteoWizard
      URL http://cruxtoolkit.sourceforge.net/download/pwiz-src-3_0_9870.tar.bz2
      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
      CONFIGURE_COMMAND ""
      BUILD_COMMAND 
        ${CMAKE_COMMAND} 
          -DBUILD_TYPE:STRING=${CMAKE_CFG_INTDIR}
          -DPREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
          -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard
          -DINCLUDE_VENDOR_LIBRARIES:BOOLEAN=${INCLUDE_VENDOR_LIBRARIES}
          -P  ${CMAKE_CURRENT_SOURCE_DIR}/build_pwiz.cmake
      INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/patches/boost ${CMAKE_CURRENT_BINARY_DIR}/include/boost
      BUILD_IN_SOURCE 1
    )
  endif (USE_LATEST_PWIZ)
  # On Windows we have to do some extra house keeping
  ExternalProject_Add_Step(
    ProteoWizard
    windows_post_install
    COMMAND
      ${CMAKE_COMMAND} 
        -DBUILD_TYPE:STRING=${CMAKE_CFG_INTDIR}
        -DPREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
        -P  ${CMAKE_CURRENT_SOURCE_DIR}/post_install_pwiz.cmake
    DEPENDEES install
  )
else (WIN32 AND NOT CYGWIN)
  if (USE_LATEST_PWIZ)
    ExternalProject_Add(
      ProteoWizard
      DOWNLOAD_COMMAND 
        ${CMAKE_COMMAND} 
          -DDOWNLOAD_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/downloads
          -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard
          -P  ${CMAKE_CURRENT_SOURCE_DIR}/download_pwiz.cmake
      DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/downloads
      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
      CONFIGURE_COMMAND ""
      BUILD_COMMAND 
        ${CMAKE_COMMAND} 
          -DBUILD_TYPE:STRING=${CMAKE_CFG_INTDIR}
          -DPREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
          -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard
          -P  ${CMAKE_CURRENT_SOURCE_DIR}/build_pwiz.cmake
      INSTALL_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/create_links.sh 
        --prefix=${CMAKE_CURRENT_BINARY_DIR} &&
        ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/patches/boost ${CMAKE_CURRENT_BINARY_DIR}/include/boost
      INSTALL_COMMAND ""
      BUILD_IN_SOURCE 1
    )
  else (USE_LATEST_PWIZ)
    ExternalProject_Add(
      ProteoWizard
      URL http://cruxtoolkit.sourceforge.net/download/pwiz-src-3_0_9870.tar.bz2
      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
      CONFIGURE_COMMAND ""
      BUILD_COMMAND 
        ${CMAKE_COMMAND} 
          -DBUILD_TYPE:STRING=${CMAKE_CFG_INTDIR}
          -DPREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
          -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard
          -P  ${CMAKE_CURRENT_SOURCE_DIR}/build_pwiz.cmake
      INSTALL_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/create_links.sh 
        --prefix=${CMAKE_CURRENT_BINARY_DIR} &&
        ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/patches/boost ${CMAKE_CURRENT_BINARY_DIR}/include/boost
      INSTALL_COMMAND ""
      BUILD_IN_SOURCE 1
    )
  endif (USE_LATEST_PWIZ)
endif (WIN32 AND NOT CYGWIN)

ExternalProject_Add(
  percolator
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
  URL http://cruxtoolkit.sourceforge.net/download/percolator-rel-3-00.tar.gz
  URL_MD5 d8a8e4c3ea261ae572b0f8b860ffefb4
  DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/downloads
  PATCH_COMMAND ${CMAKE_COMMAND} -E remove -f
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fisher/Database.h
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fisher/Database.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fisher/objects.h
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fisher/PeptideConstraint.h
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fisher/PeptideConstraint.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fisher/Peptide.h
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fisher/Peptide.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fisher/PeptideSrc.h
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fisher/PeptideSrc.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fisher/Protein.h
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fisher/Protein.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fisher/ProteinPeptideIterator.h
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fisher/ProteinPeptideIterator.cpp
    && ${CMAKE_COMMAND} -E copy_directory
       ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/src/
       ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/
  CMAKE_ARGS -DCMAKE_VERBOSE_MAKEFILE=ON
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_BUILD_TYPE=${CMAKE_CFG_INTDIR}
    -DCRUX_ROOT=${PROJECT_SOURCE_DIR}
)

if (UNIX) 
  ExternalProject_Add(
    protobuf-2.5.0
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    DOWNLOAD_COMMAND
          ${CMAKE_COMMAND} 
            -DDOWNLOAD_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/downloads
            -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src
            -P  ${CMAKE_CURRENT_SOURCE_DIR}/download_pbuffers.cmake
    CONFIGURE_COMMAND
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/protobuf-2.5.0/configure
      --prefix=${CMAKE_CURRENT_BINARY_DIR}
      --enable-shared=no
    BUILD_COMMAND make
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE 1
 )
elseif (WIN32 AND NOT CYGWIN)
  # Update Microsoft VS project files and build ProtocolBuffers using MSBuild
  ExternalProject_Add(
    protobuf-2.5.0
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    DOWNLOAD_COMMAND
          ${CMAKE_COMMAND} 
            -DDOWNLOAD_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/downloads
            -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src
            -P  ${CMAKE_CURRENT_SOURCE_DIR}/download_pbuffers.cmake
    PATCH_COMMAND ${CMAKE_COMMAND} -E copy_directory 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/protobuf/vsprojects/ 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/protobuf-2.5.0/vsprojects/ &&
      ${CMAKE_COMMAND} -E copy_directory 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/protobuf/src/ 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/protobuf-2.5.0/src/
    CONFIGURE_COMMAND ""
    BUILD_COMMAND msbuild 
      /p:Configuration=${CMAKE_CFG_INTDIR}
      /p:Platform=${MSBUILD_PLATFORM} 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/protobuf-2.5.0/vsprojects/libprotobuf.vcxproj &&
      msbuild
      /p:Configuration=${CMAKE_CFG_INTDIR}
      /p:Platform=${MSBUILD_PLATFORM} 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/protobuf-2.5.0/vsprojects/libprotoc.vcxproj &&
      msbuild
      /p:Configuration=${CMAKE_CFG_INTDIR}
      /p:Platform=${MSBUILD_PLATFORM} 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/protobuf-2.5.0/vsprojects/protoc.vcxproj
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
  )
  ExternalProject_Add_Step(
    protobuf-2.5.0
    copy_binaries
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/protobuf-2.5.0/vsprojects/${MSBUILD_PLATFORM}/${CMAKE_CFG_INTDIR}/protoc.exe
      ${CMAKE_CURRENT_BINARY_DIR}/bin/protoc.exe
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/protobuf-2.5.0/vsprojects/${MSBUILD_PLATFORM}/${CMAKE_CFG_INTDIR}/libprotobuf.lib
      ${CMAKE_CURRENT_BINARY_DIR}/lib/protobuf.lib
    DEPENDEES install
  )
endif (UNIX)

ExternalProject_Add(
  gflags-2.1.2
  DOWNLOAD_COMMAND
        ${CMAKE_COMMAND} 
          -DDOWNLOAD_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/downloads
          -DWORKING_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/build/src
          -P  ${CMAKE_CURRENT_SOURCE_DIR}/download_gflags.cmake
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
  CMAKE_ARGS -DCMAKE_VERBOSE_MAKEFILE=ON
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_BUILD_TYPE=${CMAKE_CFG_INTDIR}
  PATCH_COMMAND ""
)

if (UNIX) 
  # Get Comet source from SVN and buld using it's own makefile
  ExternalProject_Add(
    comet
    SVN_REPOSITORY svn://svn.code.sf.net/p/comet-ms/code/tags/release_2016011/
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    PATCH_COMMAND cp
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/comet/CometSearch/CometSearchManager.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/
    BUILD_COMMAND make -C MSToolkit
      COMMAND make -C CometSearch CXXFLAGS=-DCRUX 
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
    SVN_REPOSITORY http://svn.code.sf.net/p/comet-ms/code/tags/release_2016011/
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

add_dependencies(percolator ProteoWizard comet)

# On Windows we have to include the vendor libraries in the package
if (WIN32 AND NOT CYGWIN)
  if (INCLUDE_VENDOR_LIBRARIES)
    if (CMAKE_GENERATOR MATCHES "^.*Win64$")
      set(
        vendor_binaries
        ${CMAKE_CURRENT_BINARY_DIR}/lib/agtsampleinforw.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/baf2sql_c.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/BaseCommon.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/BaseDataAccess.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/BaseError.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/BaseTof.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Compression.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Data.AnalystDataProvider.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Data.CommonInterfaces.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Data.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Data.WiffReader.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.InternalRawXYProcessing.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Muni.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.RawXYProcessing.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.StructuredStorage.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Utility.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/fileio_x64.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/fregistry_x64.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/MassSpecDataReader.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/MIDAC.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/MSFileReader.XRawfile2.dll
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/pwiz_aux/msrc/utility/vendor_api/Waters/vc12_x64/cdt.dll
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/pwiz_aux/msrc/utility/vendor_api/Waters/vc12_x64/MassLynxRaw.dll
      )
    else (CMAKE_GENERATOR MATCHES "^.*Win64$")
      set(
        vendor_binaries
        ${CMAKE_CURRENT_BINARY_DIR}/lib/agtsampleinforw.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/baf2sql_c.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/BaseCommon.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/BaseDataAccess.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/BaseError.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/BaseTof.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Compression.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Data.AnalystDataProvider.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Data.CommonInterfaces.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Data.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Data.WiffReader.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.InternalRawXYProcessing.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Muni.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.RawXYProcessing.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.StructuredStorage.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Utility.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/fileio.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/fregistry.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/MassSpecDataReader.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/MIDAC.dll
        ${CMAKE_CURRENT_BINARY_DIR}/lib/MSFileReader.XRawfile2.dll
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/pwiz_aux/msrc/utility/vendor_api/Waters/vc12_x86/cdt.dll
        ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard/pwiz_aux/msrc/utility/vendor_api/Waters/vc12_x86/MassLynxRaw.dll
      )
    endif (CMAKE_GENERATOR MATCHES "^.*Win64$")
    install (FILES ${vendor_binaries} DESTINATION bin)
  endif (INCLUDE_VENDOR_LIBRARIES)
endif (WIN32 AND NOT CYGWIN)
