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

# Download and build a partial distribution of the pwiz source
# from the ProteoWizard TeamCity repository.
if (WIN32 AND NOT CYGWIN)
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
    INSTALL_COMMAND ""
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
        -P  ${CMAKE_CURRENT_SOURCE_DIR}/post_install_pwiz.cmake
    DEPENDEES install
  )
else (WIN32 AND NOT CYGWIN)
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
      --prefix=${CMAKE_CURRENT_BINARY_DIR}
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
  )
endif (WIN32 AND NOT CYGWIN)

if (WIN32 AND NOT CYGWIN)
  ExternalProject_Add(
    percolator
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    URL http://noble.gs.washington.edu/proj/crux/download/percolator-2.07.tar.gz
    URL_MD5 4a02ecc2a05d3294eae4aa44c8cb2c0e
    CMAKE_ARGS -DCMAKE_VERBOSE_MAKEFILE=ON
      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_BUILD_TYPE=${CMAKE_CFG_INTDIR}
    PATCH_COMMAND ""
  )
  # On Windows we use this to patch percolator
  ExternalProject_Add_Step(
    percolator
    windows_patch
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/windows/FidoInterface.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/FidoInterface.cpp
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/windows/ProteinFDRestimator.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/ProteinFDRestimator.cpp
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/windows/ProteinProbEstimator.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/ProteinProbEstimator.cpp
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/fix_find_boost/CMakeLists.txt 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/windows/fido/Numerical.h
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fido/Numerical.h
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/windows/fido/Numerical.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fido/Numerical.cpp
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/windows/fido/GroupPowerBigraph.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src/fido/GroupPowerBigraph.cpp
    DEPENDEES patch
    DEPENDERS configure
  )
elseif (APPLE)
  ExternalProject_Add(
    percolator
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    URL http://noble.gs.washington.edu/proj/crux/download/percolator-2.07.tar.gz
    URL_MD5 4a02ecc2a05d3294eae4aa44c8cb2c0e
    CMAKE_ARGS -DCMAKE_VERBOSE_MAKEFILE=ON 
      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR} 
      -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    INSTALL_COMMAND ""
  )
  # On Apple we use this to patch percolator
  ExternalProject_Add_Step(
    percolator
    apple_patch
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/fix_common/CommonCMake.txt 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/fix_find_boost/CMakeLists.txt 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src
    DEPENDEES patch
    DEPENDERS configure
  )
else(APPLE)
  ExternalProject_Add(
    percolator
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    URL http://noble.gs.washington.edu/proj/crux/download/percolator-2.07.tar.gz
    URL_MD5 4a02ecc2a05d3294eae4aa44c8cb2c0e
    CMAKE_ARGS -DCMAKE_VERBOSE_MAKEFILE=ON
      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    PATCH_COMMAND cp 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/fix_find_boost/CMakeLists.txt 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src
  )
endif (WIN32 AND NOT CYGWIN)

if (UNIX) 
  ExternalProject_Add(
    ProtocolBuffers
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    URL http://protobuf.googlecode.com/files/protobuf-2.5.0.tar.bz2
    URL_MD5 a72001a9067a4c2c4e0e836d0f92ece4
    DOWNLOAD_DIR downloads
    CONFIGURE_COMMAND
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProtocolBuffers/configure
      --prefix=${CMAKE_CURRENT_BINARY_DIR}
      --enable-shared=no
    BUILD_COMMAND make
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE 1
 )
elseif (WIN32 AND NOT CYGWIN)
  # Update Microsoft VS project files and build ProtocolBuffers using MSBuild
  ExternalProject_Add(
    ProtocolBuffers
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    URL http://protobuf.googlecode.com/files/protobuf-2.5.0.tar.bz2
    URL_MD5 a72001a9067a4c2c4e0e836d0f92ece4
    DOWNLOAD_DIR downloads
    PATCH_COMMAND ${CMAKE_COMMAND} -E copy_directory 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/protobuf/vs2010/ 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProtocolBuffers/vsprojects/
    CONFIGURE_COMMAND ""
    BUILD_COMMAND msbuild 
      /p:Configuration=${CMAKE_CFG_INTDIR}
      /p:Platform=Win32 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProtocolBuffers/vsprojects/protobuf.sln
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
  )
  ExternalProject_Add_Step(
    ProtocolBuffers
    copy_binaries
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProtocolBuffers/vsprojects/${CMAKE_CFG_INTDIR}/protoc.exe
      ${CMAKE_CURRENT_BINARY_DIR}/bin/protoc.exe
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProtocolBuffers/vsprojects/${CMAKE_CFG_INTDIR}/libprotobuf.lib
      ${CMAKE_CURRENT_BINARY_DIR}/lib/protobuf.lib
    DEPENDEES install
  )
endif (UNIX)

if (UNIX) 
ExternalProject_Add(
  google-gflags
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
  URL http://gflags.googlecode.com/files/gflags-2.0.tar.gz
  URL_MD5 e02268c1e90876c809ec0ffa095da864
  DOWNLOAD_DIR downloads
  CONFIGURE_COMMAND
    ${CMAKE_CURRENT_BINARY_DIR}/build/src/google-gflags/configure
    --prefix=${CMAKE_CURRENT_BINARY_DIR}
    --enable-shared=no
  BUILD_COMMAND make
  INSTALL_COMMAND make install
  BUILD_IN_SOURCE 1
)
elseif (WIN32 AND NOT CYGWIN)
  # Update Microsoft VS project files and build ProtocolBuffers using MSBuild
  ExternalProject_Add(
    google-gflags
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    URL http://gflags.googlecode.com/files/gflags-2.0.tar.gz
    URL_MD5 e02268c1e90876c809ec0ffa095da864
    DOWNLOAD_DIR downloads
    PATCH_COMMAND ${CMAKE_COMMAND} -E copy_directory 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/google-gflags/vs2010/ 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/google-gflags/vsprojects/libgflags/
    CONFIGURE_COMMAND ""
    BUILD_COMMAND msbuild 
      /p:Configuration=${CMAKE_CFG_INTDIR} 
      /p:Platform=Win32 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/google-gflags/vsprojects/libgflags/libgflags.vcxproj
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
  )
  ExternalProject_Add_Step(
    google-gflags
    copy_gflags_library
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/google-gflags/vsprojects/libgflags/${CMAKE_CFG_INTDIR}/libgflags.lib
      ${CMAKE_CURRENT_BINARY_DIR}/lib/gflags.lib
    DEPENDEES install
  )
endif (UNIX)

if (UNIX) 
  # Get Comet source from SVN and buld using it's own makefile
  ExternalProject_Add(
    comet
    SVN_REPOSITORY http://svn.code.sf.net/p/comet-ms/code/tags/release_2014011
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    PATCH_COMMAND cp
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/comet/MSToolkit/src/MSToolkit/MSReader.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/MSToolkit/src/MSToolkit
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
    SVN_REPOSITORY http://svn.code.sf.net/p/comet-ms/code/tags/release_2014011
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
        -DPROJ_FILE:STRING=comet.sln
        -P ${CMAKE_CURRENT_SOURCE_DIR}/build-comet.cmake
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
  )  
  ExternalProject_Add_Step(
    comet
    copy_comet_library
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/Win32/${CMAKE_CFG_INTDIR}/CometSearch.lib
      ${CMAKE_CURRENT_BINARY_DIR}/lib/CometSearch.lib
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/Win32/${CMAKE_CFG_INTDIR}/MSToolkitLite.lib
      ${CMAKE_CURRENT_BINARY_DIR}/lib/MSToolkitLite.lib
    DEPENDEES install
  )
endif (UNIX)

add_dependencies(percolator Proteowizard)

# On Windows we have to include the vendor libraries in the package
if (WIN32 AND NOT CYGWIN)
  set(
    vendor_binaries
    ${CMAKE_CURRENT_BINARY_DIR}/lib/agtsampleinforw.dll
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
    ${CMAKE_CURRENT_BINARY_DIR}/lib/cdt.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Compression.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Data.AnalystDataProvider.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Data.CommonInterfaces.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Data.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Data.WiffReader.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.InternalRawXYProcessing.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Muni.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.ProjectUtilities.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.RawXYProcessing.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.StructuredStorage.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Clearcore2.Utility.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/CompassXtractMS.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Compressor_12451AC8-3BF1-48e1-8E66-DA05BF3852A0.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Compressor_382D5134-9727-4be6-B6F8-754F577426D6.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Compressor_4F1EA0CB-5A6A-43ED-9562-D6D873034577.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Compressor_6BB2E64A-27A0-4575-A66A-4E312C8B9AD7.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Compressor_6EEE649E-09AA-4332-8E82-8188652E8AB5.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Compressor_7F88E97F-2038-40ba-8F3A-FCA9A9719569.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/fileio.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/fregistry.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/HSReadWrite.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Interop.DataExplorer.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Interop.EDAL.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/Interop.HSREADWRITELib.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/MassLynxRaw.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/MassSpecDataReader.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/MIDAC.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/mkl_sequential.dll
    ${CMAKE_CURRENT_BINARY_DIR}/lib/MSFileReader.XRawfile2.dll
  )
  install (FILES ${vendor_binaries} DESTINATION bin)
endif (WIN32 AND NOT CYGWIN)
