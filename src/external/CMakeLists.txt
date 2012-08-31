include(ExternalProject)
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
  Boost
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  URL http://sourceforge.net/projects/boost/files/boost/1.51.0/boost_1_51_0.tar.bz2/download 
  DOWNLOAD_DIR downloads
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ./b2
  INSTALL_COMMAND ./b2 install --prefix=${CMAKE_CURRENT_SOURCE_DIR}
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
  URL http://www.codesynthesis.com/download/xsd/3.3/xsd-3.3.0+dep.tar.bz2
  DOWNLOAD_DIR downloads
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
)
ExternalProject_Add(
  percolator
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/percolator
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/percolator
  CMAKE_ARGS -DBoost_NO_BOOST_CMAKE=1 -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_SOURCE_DIR}/percolator
)
add_dependencies(percolator XercesC)
add_dependencies(percolator Boost)
add_dependencies(percolator XSD)
add_subdirectory(MSToolkit)

