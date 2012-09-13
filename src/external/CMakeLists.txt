include(ExternalProject)

INCLUDE_DIRECTORIES(BEFORE SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/include)

SET(ENV{BOOSTROOT} ${CMAKE_CURRENT_SOURCE_DIR})
SET(ENV{LIBPATH} ${CMAKE_CURRENT_SOURCE_DIR}/lib:$ENV{LIBPATH})
SET(ENV{LD_LIBRARY_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/lib:$ENV{LD_LIBRARY_PATH})
SET(ENV{PATH} ${CMAKE_CURRENT_SOURCE_DIR}/bin:$ENV{PATH})
SET(ENV{CMAKE_PREFIX_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
SET(CPP_ARGS "-I${CMAKE_CURRENT_SOURCE_DIR}/include")
SET(LD_ARGS "-L${CMAKE_CURRENT_SOURCE_DIR}/lib")

ExternalProject_Add(
  XercesC
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  URL http://www.motorlogy.com/apache//xerces/c/3/sources/xerces-c-3.1.1.tar.gz
  URL_MD5 6a8ec45d83c8cfb1584c5a5345cb51ae
  DOWNLOAD_DIR downloads
  CONFIGURE_COMMAND 
    ${CMAKE_CURRENT_SOURCE_DIR}/build/src/XercesC/configure 
    --disable-network 
    --prefix=${CMAKE_CURRENT_SOURCE_DIR}
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
)

ExternalProject_Add(
  Boost
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  URL http://sourceforge.net/projects/boost/files/boost/1.44.0/boost_1_44_0.tar.bz2/download 
  DOWNLOAD_DIR downloads
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ./bjam --stage
  INSTALL_COMMAND ./bjam --prefix=${CMAKE_CURRENT_SOURCE_DIR} install
  BUILD_IN_SOURCE 1
)

ExternalProject_Add_Step(
  Boost
  boostrap
  DEPENDEES download
  DEPENDERS configure
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/src/Boost
  COMMAND ./bootstrap.sh
)

ExternalProject_Add(
  XSD
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  URL http://www.codesynthesis.com/download/xsd/3.3/xsd-3.3.0-2+dep.tar.bz2
  DOWNLOAD_DIR downloads
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make 
    CPPFLAGS=${CPP_ARGS}
    LDFLAGS=${LD_ARGS}
    verbose=1
  INSTALL_COMMAND make 
    CPPFLAGS="${CPP_ARGS}"
    LDFLAGS="${LD_ARGS}"
    install_prefix=${CMAKE_CURRENT_SOURCE_DIR}
    install
  BUILD_IN_SOURCE 1
)

ExternalProject_Add(
  percolator
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/percolator
  CMAKE_ARGS 
    -DBoost_NO_SYSTEM_PATHS=1
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_SOURCE_DIR}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}
)

add_dependencies(XSD XercesC)
add_dependencies(XSD Boost)
add_dependencies(percolator XercesC)
add_dependencies(percolator Boost)
add_dependencies(percolator XSD)
add_subdirectory(MSToolkit)

