#!/bin/sh
#
# script for creating symbolic links to the boost libraries (for LINUX)
#

echo "CREATING GENERICALLY NAMED SYMBOLIC LINKS TO LIBRARIES"

install_path=.
STARING_DIR=$(pwd)
for arg in $*; do
    echo "arg:$arg"
    if [[ $arg = --prefix=* ]]; then
        echo "parsing $arg"
        install_path=`echo "$arg" | cut -d"=" -f2`
    fi

done

echo "path:$install_path"



cd $install_path/lib



ln -s libboost_date_time-*-mt-s.a libboost_data_time-mt-s.a
ln -s libboost_iostreams-*-mt-s.a libboost_iostreams-mt-s.a 
ln -s libboost_serialization-*-mt-s.a libboost_serialization-mt-s.a
ln -s libboost_thread-*-mt-s.a libboost_thread-mt-s.a 
ln -s libgd-*-mt-s-2_1.a libgd-mt-s-2_1.a   
ln -s libz-*-mt-s-1_2.a libz-mt-s-1_2.a
ln -s libboost_filesystem-*-mt-s.a libboost_filesystem-mt-s.a 
ln -s libboost_regex-*-mt-s.a  libboost_regex-mt-s.a     
ln -s libboost_system-*-mt-s.a libboost_system-mt-s.a 
ln -s libfreetype-*-mt-s-2_4.a  libfreetype-mt-s-2_4.a
ln -s libpng-*-mt-s-1_5.a libpng-mt-s-1_5.a

cd $STARTING_DIR

