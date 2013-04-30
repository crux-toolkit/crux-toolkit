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
set(protwiz_args ./quickbuild.sh)
set(protwiz_args ${protwiz_args} --without-binary-msdata)
set(protwiz_args ${protwiz_args} --layout=system)
set(protwiz_args ${protwiz_args} --prefix=${CMAKE_CURRENT_SOURCE_DIR})
set(protwiz_args ${protwiz_args} runtime-link=shared)
set(protwiz_args ${protwiz_args} libraries)
set(protwiz_args ${protwiz_args} pwiz/data/common//pwiz_data_common)
set(protwiz_args ${protwiz_args} pwiz/data/identdata//pwiz_data_identdata)
set(protwiz_args ${protwiz_args} pwiz/data/identdata//pwiz_data_identdata_version)
set(protwiz_args ${protwiz_args} pwiz/data/msdata//pwiz_data_msdata)
set(protwiz_args ${protwiz_args} pwiz/data/msdata//pwiz_data_msdata_version)
set(protwiz_args ${protwiz_args} pwiz/data/proteome//pwiz_data_proteome)
set(protwiz_args ${protwiz_args} pwiz/utility/chemistry//pwiz_utility_chemistry)
set(protwiz_args ${protwiz_args} pwiz/utility/minimxml//pwiz_utility_minimxml)
set(protwiz_args ${protwiz_args} pwiz/utility/misc//pwiz_utility_misc)

add_custom_target(
  deep-clean-external
  COMMAND rm -rf build lib/* include/* bin/*
  COMMENT "Deleting all files for external tools."
)

# We want to avoid checking out the entire Proteowizard repository
# because of its massive size. We check out just the top level directory, 
# and them use 'svn update' to obtain the needed files and directories.
ExternalProject_Add(
  ProteoWizard
  SVN_REPOSITORY http://proteowizard.svn.sourceforge.net/svnroot/proteowizard/trunk/pwiz/
  SVN_REVISION -r 3950 --depth empty
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ${protwiz_args}
  INSTALL_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/create_links.sh --prefix=${CMAKE_CURRENT_SOURCE_DIR}
  BUILD_IN_SOURCE 1
)

# Since we're not checking out the entire Proteowizard tree,
# We need a separate step to use 'svn update' to pick up the files
# and directories we need.
ExternalProject_Add_Step(
  ProteoWizard
  update_selected_directories
  DEPENDEES download
  DEPENDERS configure
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/src/ProteoWizard
  COMMAND svn update --force -r 3950 libraries pwiz Jamroot.jam quickbuild.sh clean.sh
)

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

add_dependencies(XSD Proteowizard)
add_dependencies(XSD XercesC)
add_dependencies(percolator Proteowizard)
add_dependencies(percolator XercesC)
add_dependencies(percolator XSD)
add_subdirectory(MSToolkit)

