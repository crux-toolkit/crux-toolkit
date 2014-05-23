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

# Set up quickbuild arguments used by pwiz build.
# quickbuild.sh is the script provided by proteowizard 
# to bootstrap bjam and to run bjam for building the external
# libraries and Proteowizard components..
if (WIN32 AND NOT CYGWIN)
  set(protwiz_args ./quickbuild.bat)
  set(protwiz_args ${protwiz_args} --toolset=msvc-10.0)
  set(protwiz_args ${protwiz_args} --i-agree-to-the-vendor-licenses)
  set(protwiz_args ${protwiz_args} --without-mz5)
  set(protwiz_args ${protwiz_args} --prefix=${CMAKE_CURRENT_BINARY_DIR})
  set(protwiz_args ${protwiz_args} libraries)
else(WIN32 AND NOT CYGWIN)
  set(protwiz_args ./quickbuild.sh)
  set(protwiz_args ${protwiz_args} --without-binary-msdata)
  set(protwiz_args ${protwiz_args} --layout=system)
  set(protwiz_args ${protwiz_args} --prefix=${CMAKE_CURRENT_BINARY_DIR})
  set(protwiz_args ${protwiz_args} runtime-link=shared)
  set(protwiz_args ${protwiz_args} libraries)
endif (WIN32 AND NOT CYGWIN)

# Download and build a partial distribution of the pwiz source
# from the ProteoWizard TeamCity repository.
if (WIN32 AND NOT CYGWIN)
  ExternalProject_Add(
    ProteoWizard
    DOWNLOAD_COMMAND 
      ${CMAKE_COMMAND} 
        -DDOWNLOAD_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/downloads
        -DWORKING_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard
        -P  ${CMAKE_CURRENT_SOURCE_DIR}/download_pwiz.cmake
    DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/downloads
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${protwiz_args}
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
  )
else (WIN32 AND NOT CYGWIN)
  ExternalProject_Add(
    ProteoWizard
    DOWNLOAD_COMMAND 
      ${CMAKE_COMMAND} 
        -DDOWNLOAD_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/downloads
        -DWORKING_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/build/src/ProteoWizard
        -P  ${CMAKE_CURRENT_SOURCE_DIR}/download_pwiz.cmake
    DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/downloads
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${protwiz_args}
    INSTALL_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/create_links.sh 
      --prefix=${CMAKE_CURRENT_BINARY_DIR}
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
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    PATCH_COMMAND ""
  )
  # On Windows we use this to patch percolator
  ExternalProject_Add_Step(
    percolator
    windows_patch
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/fix_find_boost/CMakeLists.txt 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/percolator/src
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
      /p:Configuration=Debug 
      /p:Platform=Win32 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/ProtocolBuffers/vsprojects/protobuf.sln
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
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
      /p:Configuration=Debug 
      /p:Platform=Win32 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/google-gflags/vsprojects/libgflags/libgflags.vcxproj
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
  )
endif (UNIX)

if (UNIX) 
  # Get Comet source from SVN and buld using it's own makefile
  ExternalProject_Add(
    comet
    SVN_REPOSITORY http://svn.code.sf.net/p/comet-ms/code/branches/release_2014010
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    UPDATE_COMMAND ""
    PATCH_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/patch_comet.sh 
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_BINARY_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make CXXFLAGS=-DCRUX 
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
  )  
  # Copy the Comet and MSToolkit include files and libraries to
  # the desired directories.
  add_custom_command(
    TARGET comet
    POST_BUILD
    COMMAND 
      ${CMAKE_COMMAND} -E copy_directory 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/MSToolkit/include
      ${CMAKE_CURRENT_BINARY_DIR}/include/MSToolkit/
    COMMAND 
      ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/MSToolkit/libmstoolkitlite.a
      ${CMAKE_CURRENT_BINARY_DIR}/lib
    COMMAND 
      ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/Common.h
      ${CMAKE_CURRENT_BINARY_DIR}/include/CometSearch/Common.h
    COMMAND 
      ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/CometData.h
      ${CMAKE_CURRENT_BINARY_DIR}/include/CometSearch/CometData.h
    COMMAND 
      ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/CometSearch.h
      ${CMAKE_CURRENT_BINARY_DIR}/include/CometSearch/CometSearch.h
    COMMAND 
      ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/CometSearchManager.h
      ${CMAKE_CURRENT_BINARY_DIR}/include/CometSearch/CometSearchManager.h
    COMMAND 
      ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/CometPreprocess.h
      ${CMAKE_CURRENT_BINARY_DIR}/include/CometSearch/CometPreprocess.h
    COMMAND 
      ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/ThreadPool.h
      ${CMAKE_CURRENT_BINARY_DIR}/include/CometSearch/ThreadPool.h
    COMMAND 
      ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/Threading.h
      ${CMAKE_CURRENT_BINARY_DIR}/include/CometSearch/Threading.h
    COMMAND 
      ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/OSSpecificThreading.h
      ${CMAKE_CURRENT_BINARY_DIR}/include/CometSearch/OSSpecificThreading.h
    COMMAND 
      ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/CometInterfaces.h
      ${CMAKE_CURRENT_BINARY_DIR}/include/CometSearch/CometInterfaces.h
    COMMAND 
      ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/CometSearch/libcometsearch.a
      ${CMAKE_CURRENT_BINARY_DIR}/lib
  )
elseif (WIN32 AND NOT CYGWIN)
  set(LIBPATH "${CMAKE_CURRENT_BINARY_DIR}/lib;")
  ExternalProject_Add(
    comet
    SVN_REPOSITORY http://svn.code.sf.net/p/comet-ms/code/tags/release_2013022@232
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/build
    UPDATE_COMMAND ""
    # Copy the Comet and MSToolkit include files and libraries to
    # the desired directories.
    CONFIGURE_COMMAND 
      ${CMAKE_COMMAND} -E copy_directory 
      ${CMAKE_CURRENT_BINARY_DIR}/build/src/comet/MSToolkit/include
      ${CMAKE_CURRENT_BINARY_DIR}/include/MSToolkit/
    BUILD_COMMAND 
      ${CMAKE_COMMAND} 
        -DLIBPATH:STRING="${LIBPATH}"
        -DMSBUILD_ARGS:STRING=/p:Configuration=Release
        -DPROJ_FILE:STRING=comet.sln
        -P ${CMAKE_CURRENT_SOURCE_DIR}/build-comet.cmake
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
  )  
endif (UNIX)

add_dependencies(percolator Proteowizard)
