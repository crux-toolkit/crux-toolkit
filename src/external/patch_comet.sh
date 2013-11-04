#!/bin/sh
#
# Script for patching comet for building with CRUX.
#

echo "PATCHING COMET SOURCE FOR BUILDING WITH CRUX"

SRC_PREFIX=$1
BIN_PREFIX=$2

echo "Copying from " $1 "to" $2

cp  $SRC_PREFIX/patches/comet/Makefile  $BIN_PREFIX/build/src/comet/Makefile
cp  $SRC_PREFIX/patches/comet/CometSearch/Makefile $BIN_PREFIX/build/src/comet/CometSearch/Makefile
cp  $SRC_PREFIX/patches/comet/CometSearch/CometWriteTxt.cpp $BIN_PREFIX/build/src/comet/CometSearch/CometWriteTxt.cpp

