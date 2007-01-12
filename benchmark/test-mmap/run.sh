echo   test-seek human.fasta 4000 10 
time ./test-seek human.fasta 4000 10
echo
echo   test-seek human.fasta 4000 1000
time ./test-seek human.fasta 4000 1000
echo
echo   test-seek human.fasta 4000 100000 
time ./test-seek human.fasta 4000 100000
echo
echo   test-seek human.fasta 4000 10000000 
time ./test-seek human.fasta 4000 10000000
echo
echo

echo        human.fasta 4000 10 
time ./test human.fasta 4000 10
echo
echo        human.fasta 4000 1000
time ./test human.fasta 4000 1000
echo
echo        human.fasta 4000 100000 
time ./test human.fasta 4000 100000
echo
echo        human.fasta 4000 10000000 
time ./test human.fasta 4000 10000000
echo

echo

echo   test-seek nrdb.fasta 400000 10 
time ./test-seek nrdb.fasta 400000 10
echo
echo   test-seek nrdb.fasta 400000 1000
time ./test-seek nrdb.fasta 400000 1000
echo
echo   test-seek nrdb.fasta 400000 100000 
time ./test-seek nrdb.fasta 400000 100000
echo
echo   test-seek nrdb.fasta 400000 10000000 
time ./test-seek nrdb.fasta 400000 10000000
echo
echo

echo        nrdb.fasta 400000 10 
time ./test nrdb.fasta 400000 10
echo
echo        nrdb.fasta 400000 1000
time ./test nrdb.fasta 400000 1000
echo
echo        nrdb.fasta 400000 100000 
time ./test nrdb.fasta 400000 100000
echo
echo        nrdb.fasta 400000 10000000 
time ./test nrdb.fasta 400000 10000000
echo


