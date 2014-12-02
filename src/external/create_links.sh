#!/bin/bash
#
# script for creating symbolic links to the boost libraries (for LINUX)
#

echo "CREATING GENERICALLY NAMED SYMBOLIC LINKS TO LIBRARIES"

install_path=.
for arg in $*; do
    echo "arg:$arg"
    if [[ $arg = --prefix=* ]]; then
        echo "parsing $arg"
        install_path=`echo "$arg" | cut -d"=" -f2`
    fi

done

echo "path:$install_path"



cd $install_path/lib

ln -s -f libboost_chrono*-mt.a libboost_chrono.a
ln -s -f libboost_date_time-*-mt.a libboost_date_time-mt-s.a
ln -s -f libboost_date_time-*-mt.a libboost_date_time.a
ln -s -f libboost_iostreams-*-mt.a libboost_iostreams-mt-s.a 
ln -s -f libboost_iostreams-*-mt.a libboost_iostreams.a 
ln -s -f libboost_serialization-*-mt.a libboost_serialization-mt-s.a
ln -s -f libboost_serialization-*-mt.a libboost_serialization.a
ln -s -f libboost_thread-*-mt.a libboost_thread-mt-s.a 
ln -s -f libboost_thread-*-mt.a libboost_thread.a 
ln -s -f libgd-*-mt-2_1.a libgd-mt-s-2_1.a   
ln -s -f libgd-*-mt-2_1.a libgd.a   
[ -f libzlib.a ] &&
  { ln -s -f libzlib.a libz-mt-s-1_2.a; ln -s -f libzlib.a libz.a; } ||
  { ln -s -f libz-*-mt-1_2.a libz-mt-s-1_2.a; ln -s -f libz-*-mt-1_2.a libz.a; }
ln -s -f libboost_filesystem-*-mt.a libboost_filesystem-mt-s.a 
ln -s -f libboost_filesystem-*-mt.a libboost_filesystem.a 
ln -s -f libboost_system-*-mt.a libboost_system-mt-s.a 
ln -s -f libboost_system-*-mt.a libboost_system.a 
ln -s -f libfreetype-*-mt-2_4.a  libfreetype-mt-s-2_4.a
ln -s -f libfreetype-*-mt-2_4.a  libfreetype.a

# Copy the file needed for Gregorian date time

cd ..

cp $install_path/build/src/ProteoWizard/libraries/boost_1_54_0/boost/date_time/gregorian_calendar.ipp include/boost/date_time
