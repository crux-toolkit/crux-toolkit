# name files
srcms2=090306-20cm-yeast-2h-01.ms2; 
listfile=file-1.exact; 
rm -f select-spectra.ms2; 

# list scans, those with one charge state and those with two
cut -f3 $listfile | sort | uniq -d > two-charge-scan-list
cut -f3 $listfile | sort | uniq -u > one-charge-scan-list

# get spec with multiple charge state matches
for s in `cat two-charge-scan-list`
do rm -f scan.sqt
   ./crux-get-ms2-spectrum --verbosity 0 $s $srcms2 scan.sqt; 
   sed -e '/^[0-9]/ s/\(\.[0-9]\)[0-9]*$/\1/' scan.sqt >> select-spectra.ms2
done

# get spec with only one charge state (remove z 2 or 3 before adding to ms2)
for s in `cat one-charge-scan-list`; 
do rm -f scan.sqt; 
   ./crux-get-ms2-spectrum --verbosity 0 $s $srcms2 scan.sqt; 
   z=`grep "	$s	" $listfile | cut -f 4`; 
   sed -n -e '1 p' -ne "/^Z	$z/ p" -ne '/^[0-9]/ s/\(\.[0-9]\)[0-9]*$/\1/ p' scan.sqt >> select-spectra.ms2; 
done
