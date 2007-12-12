cd ..
rm -f ../../bin/create_psm_files
make
cd test
 rm -fr output
# gdb --args ~/crux/bin/create_psm_files peptides.txt 13002 test-3.ms2 output --verbosity 100
~/crux/bin/create_psm_files peptides.txt 13002 test-3.ms2 output --verbosity 100
cd output
pfile_create -i 0.prepfile -o 0.pfile -l 3 -f 9 -b
pfile_create -i 1.prepfile -o 1.pfile -l 3 -f 9 -b
pfile_print -i 0.pfile | head
pfile_print -i 1.pfile | head
