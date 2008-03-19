# ../../crux-create-index yeast-shuffled.fasta yeast-shuffled.index --parameter-file crux.params

# ../../crux-create-index yeast.fasta yeast.index --parameter-file crux.params


export MS2_FILE="60cm.ms2"

rm -f *.csm
rm -f target.sqt
./crux-search-for-matches $MS2_FILE yeast.index --overwrite T \
  --use-index T --score-type xcorr-logp --number-decoy-set 0 --verbosity 40 \
  --parameter-file crux.params 
mkdir -p target
mv *.csm target
./crux-analyze-matches target yeast.index --use-index T --overwrite T \
  --sqt-output-file target.sqt \
  --algorithm none --parameter-file crux.params > target.proteins


rm -f decoy.sqt
./crux-search-for-matches $MS2_FILE yeast-shuffled.index --overwrite T \
  --use-index T --sqt-output-file decoy.sqt \
  --score-type xcorr-logp --number-decoy-set 0 \
  --parameter-file crux.params 
mkdir -p decoy
mv *.csm decoy
./crux-analyze-matches decoy yeast-shuffled.index --use-index T  \
  --overwrite T --sqt-output-file decoy.sqt \
  --algorithm none --parameter-file crux.params > shuffle.proteins

