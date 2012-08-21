project(MSToolkit)

cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_property(
      DIRECTORY
      PROPERTY
      COMPILE_DEFINITIONS
      _CRT_SECURE_NO_WARNINGS
      _SCL_SECURE_NO_WARNINGS
      COMPILED_FROM_DSP
      XML_BUILD_EXPAT
      XML_STATIC
    )
else (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_property(
      DIRECTORY
      PROPERTY
      COMPILE_DEFINITIONS
      _FILE_OFFSET_BITS=64
      _LARGEFILE_SOURCE
      GCC
      HAVE_EXPAT_CONFIG_H
    )
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

include_directories(include)
add_library(
  MSToolkit
  STATIC
  src/zLib-1.2.5/adler32.c
  src/zLib-1.2.5/compress.c
  src/zLib-1.2.5/crc32.c
  src/zLib-1.2.5/deflate.c
  src/zLib-1.2.5/inffast.c
  src/zLib-1.2.5/inflate.c
  src/zLib-1.2.5/infback.c
  src/zLib-1.2.5/inftrees.c
  src/zLib-1.2.5/trees.c
  src/zLib-1.2.5/uncompr.c
  src/zLib-1.2.5/zutil.c
  src/mzParser/mzp_base64.cpp
  src/mzParser/BasicSpectrum.cpp
  src/mzParser/mzParser.cpp
  src/mzParser/RAMPface.cpp
  src/mzParser/saxhandler.cpp
  src/mzParser/saxmzmlhandler.cpp
  src/mzParser/saxmzxmlhandler.cpp
  src/mzParser/Czran.cpp
  src/expat-2.0.1/xmlparse.c
  src/expat-2.0.1/xmlrole.c
  src/expat-2.0.1/xmltok.c
  src/expat-2.0.1/xmltok_impl.c
  src/expat-2.0.1/xmltok_ns.c
  src/sqlite-3.7.7.1/sqlite3.c
  src/MSToolkit/Spectrum.cpp
  src/MSToolkit/MSReader.cpp
  src/MSToolkit/MSObject.cpp
)
