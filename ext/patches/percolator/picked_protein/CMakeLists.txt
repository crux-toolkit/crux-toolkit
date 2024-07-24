include_directories(${PERCOLATOR_SOURCE_DIR}/src)
link_directories(${PERCOLATOR_SOURCE_DIR}/src)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
file(GLOB PICKED_PROTEIN_SOURCES PickedProteinCaller.cpp Database.cpp Protein.cpp ProteinPeptideIterator.cpp Peptide.cpp PeptideSrc.cpp PeptideConstraint.cpp ../Option.cpp ../Globals.cpp ../MyException.cpp ../Logger.cpp)
add_library(picked_protein STATIC ${PICKED_PROTEIN_SOURCES})
