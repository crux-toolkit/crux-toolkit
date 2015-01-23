cd ../
# make clean
# rm -f ../../bin/create_ion_files
make
cd test
rm -fr output
../../../../bin/create_psm_files peptides.txt 16 test-3.ms2 output --verbosity 100 --starting-sentence-idx 0
cd output/
pfile_create -o 0.pfile -i 0.prepfile -b -f 3 -l 9 
pfile_create -o 1.pfile -i 1.prepfile -b -f 3 -l 9 
echo
pfile_print -i 0.pfile | head
echo
pfile_print -i 1.pfile | head
