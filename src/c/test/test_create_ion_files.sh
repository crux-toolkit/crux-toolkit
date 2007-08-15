# cd ../
# make clean
# rm -f ../../bin/create_ion_files
# make
# cd test
 rm -fr output
../../../bin/create_ion_files peptides.txt 13002 test-3.ms2 output --verbosity 100
