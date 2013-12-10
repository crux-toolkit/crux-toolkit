include_directories(${PERCOLATOR_SOURCE_DIR}/src)
link_directories(${PERCOLATOR_SOURCE_DIR}/src)
if (WIN32 AND NOT CYGWIN)
  add_definitions(-DNOMINMAX)
endif (WIN32 AND NOT CYGWIN)

file(GLOB FIDO_SOURCES Set.cpp Vector.cpp PackedVector.cpp Numerical.cpp Matrix.cpp PackedMatrix.cpp GramSchmidt.cpp LPSolver.cpp ForwardSampler.cpp Random.cpp ProteinIdentifier.cpp BasicBigraph.cpp BasicGroupBigraph.cpp GroupPowerBigraph.cpp PivdoSplitter.cpp)
add_library(fido ${FIDO_SOURCES})

