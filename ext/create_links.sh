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

ln -s -f libboost_chrono-*.a libboost_chrono.a
ln -s -f libboost_date_time-**.a libboost_date_time.a
ln -s -f libboost_filesystem-*.a libboost_filesystem.a 
ln -s -f libboost_iostreams-*.a libboost_iostreams.a 
ln -s -f libboost_serialization-*.a libboost_serialization.a
ln -s -f libboost_system-*.a libboost_system.a 
ln -s -f libboost_thread-*.a libboost_thread.a 
ln -s -f libfreetype-*.a  libfreetype.a

