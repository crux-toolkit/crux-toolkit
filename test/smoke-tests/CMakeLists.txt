cmake_minimum_required(VERSION 3.15)
cmake_policy(SET CMP0091 NEW)
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.24")
  cmake_policy(SET CMP0135 NEW)
endif()

# Copy files for smoke tests in case we are building out of source.
file(COPY runall DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY crux-test.pl DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY crux-test.cmds DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY clean.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY compare-by-field.pl DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY test.fasta DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY for-sequest-comparison.fasta DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY small-yeast.fasta DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY duplicates.fasta DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY demo.ms2 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY test.ms2 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY test-no-z1.ms2 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY small.raw2xml.mzXML DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY test.mgf DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY hardklor.test.ms1 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY bullseye.test.ms2 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY params DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY good_results DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# LFQ smoke-test inputs
file(COPY lfq-psms-no-norm.tsv DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY lfq-psms-with-norm.tsv DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY lfq-specfile-replicates.tsv DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# mzML files shared with the unit-test suite
file(COPY
  ${CMAKE_SOURCE_DIR}/test/unit-tests/test-data/sliced-mzml.mzML
  ${CMAKE_SOURCE_DIR}/test/unit-tests/test-data/20100614_Velos1_TaGe_SA_K562_3.mzML
  ${CMAKE_SOURCE_DIR}/test/unit-tests/test-data/20100614_Velos1_TaGe_SA_K562_4.mzML
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(smoke-tests COMMAND ./runall)
