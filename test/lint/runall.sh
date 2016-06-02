#!/bin/bash
hostname
date
echo PID=$$

# These three options make it harder to for some part of a script to
# fail without you recognizing it. nounset means that references to an
# unset variable result in an error. This means that you can no longer
# do stuff like this if VAR is potentially unset, because "$VAR" returns
# an error rather than "":
#
# if [ "$VAR" ]; then
#
# fi
#
# To explicitly indicate that you are OK with the variable potentially
# being empty you can instead use ${VAR:-}.

# To avoid errexit for a single command, use "|| true", e.g.,
#    diff foo foobar || true

set -o nounset
set -o pipefail
set -o errexit 
#set -o xtrace

prevDirectory=""

# Initialize the output file.
results=lint-results.html
echo "<html><title>Crux lint test results</title>" > $results
echo "<body><h1>Crux lint test results</h1>" >> $results
echo "<p>" >> $results
date >> $results
echo "</p>" >> $results

for filename in `find ../../src \( -name \*.cpp -o -name \*.h \)`; do

  truncatedName=`echo $filename | sed 's|../../src/||g'`
  echo $truncatedName

  # If we're in a new directory, make a header.
  directory=`dirname $truncatedName`
  if [[ $directory != $prevDirectory ]]; then
    prevDirectory=$directory
    echo \<h2\>$directory\</h2\> >> $results
    mkdir -p $directory
  fi
  outFile=$truncatedName.txt
  
  
  # For some reason, cruxlint outputs to stderr rather than stdout.
  ./cruxlint.sh $filename 2> $outFile || true

  errorCount=`grep "Total errors found" $outFile | awk '{print $4}'`

  if [[ $errorCount == "0" ]]; then
    echo $filename >> $results
  else
    echo \<a href=\"$outFile\"\> $filename \</a\> >> $results
  fi
  echo \($errorCount errors found\) \<br\> >> $results
done

echo "</ol></body></html>" >> $results
