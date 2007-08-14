 cd ../
 make clean
 rm -f ../../bin/create_ion_files
 make 
 cd test
 rm -fr output
../../../bin/create_ion_files LLRKLEAMAPK 18 test.ms2 output --verbosity 0
