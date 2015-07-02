#!/bin/bash

# assume that you have been running runall and that all the preliminaries are
# in place

# Get the name of the tests to run
testnames="$*"

echo "Running test $testnames"

# pull the tests out of the command file
rm -f select-tests
for t in $testnames
do
 grep " $t " crux-test.cmds >> select-tests

 if [ $? != 0 ] ; then
    echo "Could not find test $t."
    exit
 fi
done

# run the test
./crux-test.pl -p ../../ select-tests 1>out-select 2>error-select

# print success
if [ $? == 25 ] ; then
  echo "Test aborted"
else
  tail -n6 out-select
fi
