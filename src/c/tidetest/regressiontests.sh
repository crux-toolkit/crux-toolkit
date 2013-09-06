#!/bin/bash

scriptdir=$(dirname "$BASH_SOURCE")
testbinary="$scriptdir/tide-regression"
binarycpp="$scriptdir/TideRegressionTestDriver.cpp $scriptdir/TideRegressionTest.cpp"
crux="$scriptdir/../crux"

if ! [ -f "$testbinary" ]; then
  g++ -o "$testbinary" $binarycpp
fi

"$testbinary" "$crux"

