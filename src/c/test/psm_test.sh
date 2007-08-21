# cd ../
# make clean
# rm -f ../../bin/create_ion_files
# make
# cd test
cd ..
make
cd -
rm -fr output
../../../bin/create_psm_files peptides.txt 13002 test-3.ms2 output --verbosity 100
cd output/
# pfile_create -o 0.pfile -i 0.prepfile -b -f 1 -l 1 
# pfile_create -o 1.pfile -i 1.prepfile -b -f 1 -l 1 
pfile_create -o 0.pfile -i 0.prepfile -b -f 3 -l 9 
# pfile_create -o 1.pfile -i 1.prepfile -b -f 3 -l 9
pfile_print -i 0.pfile
od -x 0.prepfile
