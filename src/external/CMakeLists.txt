INCLUDE_DIRECTORIES(BEFORE SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/include)
SET(ENV{LIBPATH} ${CMAKE_CURRENT_SOURCE_DIR}/lib:$ENV{LIBPATH})
SET(ENV{PATH} ${CMAKE_CURRENT_SOURCE_DIR}/bin:$ENV{PATH})
include(ExternalProject)
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR})
ExternalProject_Add(
  XercesC
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  URL http://www.trieuvan.com/apache//xerces/c/3/sources/xerces-c-3.1.1.tar.gz 
  DOWNLOAD_DIR downloads
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/build/src/XercesC/configure --disable-net --prefix=${CMAKE_CURRENT_SOURCE_DIR}
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
)
ExternalProject_Add(
  BJam
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  URL http://sourceforge.net/projects/boost/files/boost-jam/3.1.13/boost-jam-3.1.13.tgz/download  
  DOWNLOAD_DIR downloads
  CONFIGURE_COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/bin
  BUILD_COMMAND ./build.sh
  INSTALL_COMMAND cp bin.linux/bjam ${CMAKE_CURRENT_SOURCE_DIR}/bin
  BUILD_IN_SOURCE 1
)
ExternalProject_Add(
  Boost
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  URL http://sourceforge.net/projects/boost/files/boost/1.37.0/boost_1_37_0.tar.bz2/download 
  DOWNLOAD_DIR downloads
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bin/bjam --without-python stage
  INSTALL_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bin/bjam install --prefix=${CMAKE_CURRENT_SOURCE_DIR}
  BUILD_IN_SOURCE 1
)
#ExternalProject_Add(
#  Boost
#  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
#  URL http://sourceforge.net/projects/boost/files/boost/1.44.0/boost_1_44_0.tar.bz2/download 
#  DOWNLOAD_DIR downloads
#  CONFIGURE_COMMAND ""
#  BUILD_COMMAND bjam
#  INSTALL_COMMAND bjam install --prefix=${CMAKE_CURRENT_SOURCE_DIR}
#  BUILD_IN_SOURCE 1
#)
#ExternalProject_Add_Step(
#  Boost
#  boostrap
#  DEPENDEES download
#  DEPENDERS configure
#  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/src/Boost
#  COMMAND ./bootstrap.sh
#)
ExternalProject_Add(
  XSD
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  URL http://www.codesynthesis.com/download/xsd/3.3/xsd-3.3.0+dep.tar.bz2
  DOWNLOAD_DIR downloads
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make CPPFLAGS=-I${CMAKE_CURRENT_SOURCE_DIR}/include LDFLAGS=-L${CMAKE_CURRENT_SOURCE_DIR}/lib
  BUILD_IN_SOURCE 1
)
ExternalProject_Add(
  percolator
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/percolator
  CMAKE_ARGS -DBoost_NO_BOOST_CMAKE=1 -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_SOURCE_DIR}
)
add_dependencies(Boost BJam)
add_dependencies(XSD XercesC)
add_dependencies(XSD Boost)
add_dependencies(percolator XercesC)
add_dependencies(percolator Boost)
add_dependencies(percolator XSD)
add_subdirectory(MSToolkit)

