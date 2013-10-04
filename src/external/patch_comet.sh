#!/bin/sh
#
# Script for patching comet for building with CRUX.
#

echo "PATCHING COMET SOURCE FOR BUILDING WITH CRUX"

SRC_PREFIX=$1
BIN_PREFIX=$2

echo "Copying from " $1 "to" $2

cp  $SRC_PREFIX/patches/comet/Makefile  $BIN_PREFIX/build/src/comet/Makefile
cp  $SRC_PREFIX/patches/comet/CometSearch/CometPreprocess.h $BIN_PREFIX/build/src/comet/CometSearch/CometPreprocess.h
cp  $SRC_PREFIX/patches/comet/CometSearch/Common.h $BIN_PREFIX/build/src/comet/CometSearch/Common.h
cp  $SRC_PREFIX/patches/comet/CometSearch/Makefile $BIN_PREFIX/build/src/comet/CometSearch/Makefile
cp $SRC_PREFIX/patches/comet/MSToolkit/Makefile $BIN_PREFIX/build/src/comet/MSToolkit/Makefile
cp $SRC_PREFIX/patches/comet/MSToolkit/include/H5public.h  $BIN_PREFIX/build/src/comet/MSToolkit/include/H5public.h
cp $SRC_PREFIX/patches/comet/MSToolkit/include/mzParser.h $BIN_PREFIX/build/src/comet/MSToolkit/include/mzParser.h
