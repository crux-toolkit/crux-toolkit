#!/bin/bash

# testruntime.sh [tide|crux] [index|search]

timecmd="/usr/bin/time"

scriptdir=$(dirname "$BASH_SOURCE")
cruxdir="$scriptdir/.."
tidedir="$scriptdir/../tide"

fasta="$scriptdir/worm.fasta"
spectra="$scriptdir/worm-06-10000.spectrumrecords"
cruxindexdir="$scriptdir/wormcrux"
tideindexpep="$scriptdir/wormtide.pepix"
tideindexprot="$scriptdir/wormtide.protix"
tideindexaux="$scriptdir/wormtide.auxlocs"
tideindexfiles="$tideindexpep $tideindexprot $tideindexaux $tideindexpep.nopeaks.tmp"

crux="$cruxdir/crux"
cruxoutput="$scriptdir/crux-output"
tideindex="$tidedir/tide-index"
tidesearch="$tidedir/tide-search"

if ! [ -f "$crux" ]; then
  echo "$crux not found"
  exit 1
elif ! [ -f "$tideindex" ]; then
  echo "$tideindex not found"
  exit 1
elif ! [ -f "$tidesearch" ]; then
  echo "$tidesearch not found"
  exit 1
fi

runtideindex=0
runtidesearch=0
runcruxindex=0
runcruxsearch=0

if [ -z "$1" ]; then
  runtideindex=1
  runtidesearch=1
  runcruxindex=1
  runcruxsearch=1
elif [ "$1" == "tide" ]; then
  if [ -z "$2" ]; then
    runtideindex=1
    runtidesearch=1
  elif [ "$2" == "index" ]; then
    runtideindex=1
  elif [ "$2" == "search" ]; then
    runtidesearch=1
  else
    echo "invalid command $2, should be index or search"
    exit 1
  fi
elif [ "$1" == "crux" ]; then
  if [ -z "$2" ]; then
    runcruxindex=1
    runcruxsearch=1
  elif [ "$2" == "index" ]; then
    runcruxindex=1
  elif [ "$2" == "search" ]; then
    runcruxsearch=1
  else
    echo "invalid command $2, should be index or search"
    exit 1
  fi
else
  echo "invalid program $1, should be crux or tide"
  exit 1
fi

if [ $runtideindex -eq 1 ]; then
  for file in $tideindexfiles; do
    if [ -f "$file" ]; then
      rm "$file"
    fi
  done
  echo -e "\e[1;31mRunning tide-index (Tide)...\e[0m"
  $timecmd -f "[%U/%S]" "$tideindex" --peptides="$tideindexpep" --proteins="$tideindexprot" --aux_locations="$tideindexaux" --enzyme=trypsin --mods_spec=C+57.0214637206 --fasta="$fasta"
fi

if [ $runtidesearch -eq 1 ]; then
  echo -e "\e[1;31mRunning tide-search (Tide)...\e[0m"
  $timecmd -f "[%U/%S]" "$tidesearch" --peptides="$tideindexpep" --proteins="$tideindexprot" --spectra="$spectra" 1>/dev/null
fi

if { [ $runcruxindex -eq 1 ] || [ $runcruxsearch -eq 1 ]; } && [ -d "$cruxoutput" ]; then
  rm -rf "$cruxoutput"
fi

if [ $runcruxindex -eq 1 ]; then
  if [ -d "$cruxindexdir" ]; then
    rm -rf "$cruxindexdir"
  fi
  echo -e "\e[1;31mRunning tide-index (Crux)...\e[0m"
  $timecmd -f "[%U/%S]" "$crux" tide-index --max-length 50 --min-length 6 --max-mass 7200 --min-mass 200 --monoisotopic-precursor F --use-flanking-peaks T --mods-spec C+57.0214637206 --decoy-format none --output-dir "$cruxoutput" "$fasta" "$cruxindexdir"
fi

if [ $runcruxsearch -eq 1 ]; then
  echo -e "\e[1;31mRunning tide-search (Crux)...\e[0m"
  $timecmd -f "[%U/%S]" "$crux" tide-search --precursor-window 3 --precursor-window-type mass --compute-sp F --min-peaks 0 --top-match 5 --txt-output F --output-dir "$cruxoutput" "$spectra" "$cruxindexdir"
fi

