cd ..
# rm -f ../../bin/create_psm_files
make
cd test
rm -fr output
rm -fr paired single
# gdb --args ~/crux/bin/create_psm_files peptides.txt 13002 test-3.ms2 output single --verbosity 100
~/crux/bin/create_psm_files \
  peptides.txt 13002 test-3.ms2 single single --verbosity 100 

cd single
pfile_create -i 0.prepfile -o 0.pfile -l 3 -f 9 -b
pfile_print -i 0.pfile | head
cd ..

# gdb --args ~/crux/bin/create_psm_files peptides.txt 13002 test-3.ms2 output paired --verbosity 100

~/crux/bin/create_psm_files \
   peptides.txt 13002 test-3.ms2 paired paired --verbosity 100 

cd paired
for idx in `echo 0 1 2 3 4 5 6 7 8 9 10 11 12`;{
  echo $idx;
  pfile_create -i $idx.prepfile -o $idx.pfile -l 12 -f 6 -b
  pfile_print -i $idx.pfile | head
}

rm -fr paired single
