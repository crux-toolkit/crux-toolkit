include(ExternalProject)

INCLUDE_DIRECTORIES(BEFORE SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/include)

SET(ENV{BOOSTROOT} ${CMAKE_CURRENT_SOURCE_DIR})
SET(ENV{LIBPATH} ${CMAKE_CURRENT_SOURCE_DIR}/lib:$ENV{LIBPATH})
SET(ENV{LD_LIBRARY_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/lib:$ENV{LD_LIBRARY_PATH})
SET(ENV{PATH} ${CMAKE_CURRENT_SOURCE_DIR}/bin:$ENV{PATH})
SET(ENV{CMAKE_PREFIX_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
SET(CPP_ARGS "-I${CMAKE_CURRENT_SOURCE_DIR}/include")
SET(LD_ARGS "-L${CMAKE_CURRENT_SOURCE_DIR}/lib")
SET(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})

# Set up quickbuild arguments
# quickbuild.sh is the script provided by proteowizard 
# to bootstrap bjam and to run bjam for building the external
# libraries and Proteowizard components..
if (WIN32 AND NOT CYGWIN)
  set(protwiz_args ./quickbuild.bat)
  set(protwiz_args ${protwiz_args} --without-binary-msdata)
  set(protwiz_args ${protwiz_args} --prefix=${CMAKE_CURRENT_SOURCE_DIR})
  set(protwiz_args ${protwiz_args} libraries)
else(WIN32 AND NOT CYGWIN)
  set(protwiz_args ./quickbuild.sh)
  set(protwiz_args ${protwiz_args} --without-binary-msdata)
  set(protwiz_args ${protwiz_args} --layout=system)
  set(protwiz_args ${protwiz_args} --prefix=${CMAKE_CURRENT_SOURCE_DIR})
  set(protwiz_args ${protwiz_args} runtime-link=shared)
  set(protwiz_args ${protwiz_args} libraries)
endif (WIN32 AND NOT CYGWIN)

if (UNIX)
add_custom_target(
  deep-clean-external
  COMMAND rm -rf build lib/* include/* bin/*
  COMMENT "Deleting all files for external tools."
)
endif (UNIX)

# Download a partial distribution of pwiz r4837 stored
# as a tar ball on the Noble Lab server.
if (WIN32 AND NOT CYGWIN)
ExternalProject_Add(
  ProteoWizard
  URL http://noble.gs.washington.edu/proj/crux/download/pwiz-4837.tar.gz
  URL_MD5 b59e972e193d8d78e89f2a74d2fc94ba
  DOWNLOAD_DIR downloads
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ${protwiz_args}
  INSTALL_COMMAND ""
  BUILD_IN_SOURCE 1
)
else (WIN32 AND NOT CYGWIN)
ExternalProject_Add(
  ProteoWizard
  URL http://noble.gs.washington.edu/proj/crux/download/pwiz-4837.tar.gz
  URL_MD5 b59e972e193d8d78e89f2a74d2fc94ba
  DOWNLOAD_DIR downloads
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ${protwiz_args}
  INSTALL_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/create_links.sh --prefix=${CMAKE_CURRENT_SOURCE_DIR}
  BUILD_IN_SOURCE 1
)
endif (WIN32 AND NOT CYGWIN)

ExternalProject_Add(
  XercesC
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  URL http://www.motorlogy.com/apache//xerces/c/3/sources/xerces-c-3.1.1.tar.gz
  URL_MD5 6a8ec45d83c8cfb1584c5a5345cb51ae
  DOWNLOAD_DIR downloads
  CONFIGURE_COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/build/src/XercesC/configure 
    --prefix=${CMAKE_CURRENT_SOURCE_DIR}
    --enable-transcoder-iconv
    --disable-network 
    --disable-threads
    --disable-shared
  BUILD_COMMAND make LD_ARGS="-L${CMAKE_CURRENT_SOURCE_DIR}/lib -framework CoreServices"
  BUILD_IN_SOURCE 1
)

if (APPLE OR CYGWIN)
  ExternalProject_Add(
    XSD
    PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
    URL http://www.codesynthesis.com/download/xsd/3.3/xsd-3.3.0-2+dep.tar.bz2
    URL_MD5 7ecee5029c119ef678bbe5ca4ef8a288
    DOWNLOAD_DIR downloads
    PATCH_COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/patches/XSD/mac_xdr_fix/xdr-stream-insertion.hxx ${CMAKE_CURRENT_SOURCE_DIR}/patches/XSD/mac_xdr_fix/xdr-stream-extraction.hxx ${CMAKE_CURRENT_SOURCE_DIR}/build/src/XSD/xsd/libxsd/xsd/cxx/tree/
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make CPPFLAGS=${CPP_ARGS} LDFLAGS=${LD_ARGS} verbose=1 .LIBPATTERNS=lib%.a 
    INSTALL_COMMAND make CPPFLAGS="${CPP_ARGS}" LDFLAGS="${LD_ARGS}" install_prefix=${CMAKE_CURRENT_SOURCE_DIR} install
    BUILD_IN_SOURCE 1
  )
else(APPLE OR CYGWIN)
  ExternalProject_Add(
    XSD
    PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
    URL http://www.codesynthesis.com/download/xsd/3.3/xsd-3.3.0-2+dep.tar.bz2
    URL_MD5 7ecee5029c119ef678bbe5ca4ef8a288
    DOWNLOAD_DIR downloads 
    CONFIGURE_COMMAND "" 
    BUILD_COMMAND make CPPFLAGS=${CPP_ARGS} LDFLAGS=${LD_ARGS} verbose=1 .LIBPATTERNS=lib%.a 
    INSTALL_COMMAND make CPPFLAGS="${CPP_ARGS}" LDFLAGS="${LD_ARGS}" install_prefix=${CMAKE_CURRENT_SOURCE_DIR} install
    BUILD_IN_SOURCE 1
  )
endif(APPLE OR CYGWIN)


if (APPLE)
ExternalProject_Add(
  percolator
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  URL http://noble.gs.washington.edu/proj/crux/download/percolator_2_04_full_src.tar.gz
  URL_MD5 5c328688493d57c5b41d413338dfcb66
  CMAKE_ARGS -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR} -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_SOURCE_DIR}
  PATCH_COMMAND cp  ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/fix_xml_write/Caller.h ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/fix_var_mods/DataSet.cpp ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/fix_m_pi_def/EludeModel.cpp ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/fix_find_boost/CMakeLists.txt ${CMAKE_CURRENT_SOURCE_DIR}/build/src/percolator/src
)
else(APPLE)
ExternalProject_Add(
  percolator
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  URL http://noble.gs.washington.edu/proj/crux/download/percolator_2_04_full_src.tar.gz
  URL_MD5 5c328688493d57c5b41d413338dfcb66
  CMAKE_ARGS -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR} -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_SOURCE_DIR}
  PATCH_COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/fix_xml_write/Caller.h ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/fix_var_mods/DataSet.cpp ${CMAKE_CURRENT_SOURCE_DIR}/patches/percolator/fix_find_boost/CMakeLists.txt ${CMAKE_CURRENT_SOURCE_DIR}/build/src/percolator/src
)
endif(APPLE)

if (UNIX) 
ExternalProject_Add(
  ProtocolBuffers
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  URL http://protobuf.googlecode.com/files/protobuf-2.5.0.tar.bz2
  URL_MD5 a72001a9067a4c2c4e0e836d0f92ece4
  DOWNLOAD_DIR downloads
  CONFIGURE_COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/build/src/ProtocolBuffers/configure
    --prefix=${CMAKE_CURRENT_SOURCE_DIR}
    --enable-shared=no
  BUILD_COMMAND make
  INSTALL_COMMAND make install
  BUILD_IN_SOURCE 1
)
endif (UNIX)

if (UNIX) 
ExternalProject_Add(
  google-gflags
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  URL http://gflags.googlecode.com/files/gflags-2.0.tar.gz
  URL_MD5 e02268c1e90876c809ec0ffa095da864
  DOWNLOAD_DIR downloads
  CONFIGURE_COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/build/src/google-gflags/configure
    --prefix=${CMAKE_CURRENT_SOURCE_DIR}
   --enable-shared=no
  BUILD_COMMAND make
  INSTALL_COMMAND make install
  BUILD_IN_SOURCE 1
)
endif (UNIX)

ExternalProject_Add(
  comet
  SVN_REPOSITORY svn://svn.code.sf.net/p/comet-ms/code/trunk/comet-ms
  SVN_REVISION -r 147
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build/src/comet
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make 
  INSTALL_COMMAND ""
  PATCH_COMMAND patch -d ${CMAKE_CURRENT_SOURCE_DIR}/build/src/comet/src/comet -p0 -N -f -i ${CMAKE_CURRENT_SOURCE_DIR}/patches/comet/compile_patch.diff
  BUILD_IN_SOURCE 1
)  


add_dependencies(XSD Proteowizard)
add_dependencies(XSD XercesC)
add_dependencies(percolator Proteowizard)
add_dependencies(percolator XercesC)
add_dependencies(percolator XSD)
