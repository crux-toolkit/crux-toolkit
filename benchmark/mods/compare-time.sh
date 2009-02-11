#!/bin/bash

# Script to benchmark the time to search with modifications
# as compared to a search without them.  Get a n-fold increase
# measurement for crux and for sequest

# Get input files
ln -sf ../../doc/user/data/demo.ms2 .
ln -sf ../../src/c/crux .
ln -sf ../../results/paper-figure/score/sequest27 .
chmod +x ../../results/paper-figure/score/sequest27
ln -sf ../../results/paper-figure/score/ms22dta.pl .
chmod +x ../../results/paper-figure/score/ms22dta.pl
ln -sf ../../doc/user/data/small-yeast.fasta .

### For SEQUEST ###

# write dtas
if [ ! -e dtas ]
then
 mkdir dtas
 ./ms22dta.pl demo.ms2
 mv *dta dtas/
fi 

# search without mods
echo "Sequest no mods"
rm -f sequest.params
ln -s sequest-nomod.params sequest.params
time ./sequest27 dtas/*.dta 1>/dev/null

# search with mods
echo "Sequest one mod"
rm -f sequest.params
ln -s sequest-mod.params sequest.params
time ./sequest27 dtas/*dta 1>/dev/null

### For crux ###

# create index
./crux create-index --parameter-file crux-nomod.params small-yeast.fasta index

# search without mods
echo "Crux no mods"
time ./crux search-for-matches --parameter-file crux-nomod.params demo.ms2 index/
#time ./crux search-for-matches --parameter-file crux-nomod.params demo.ms2 small-yeast.fasta

# search with mods
echo "crux one mods"
time ./crux search-for-matches --parameter-file crux-mod.params demo.ms2 index/
#time ./crux search-for-matches --parameter-file crux-mod.params demo.ms2 small-yeast.fasta
