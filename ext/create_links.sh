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


mkdir -p $install_path/lib
cd $install_path/lib

find $install_path/build/src/ProteoWizard -name '*\.a' -exec cp {} . \;
ln -s -f libboost_chrono-*gcc*a libboost_chrono.a
ln -s -f libboost_filesystem-*gcc*.a libboost_filesystem.a 
ln -s -f libboost_iostreams-*gcc*.a libboost_iostreams.a 
ln -s -f libboost_serialization-*gcc*.a libboost_serialization.a
ln -s -f libboost_system-*gcc*.a libboost_system.a 
ln -s -f libboost_thread-*gcc*.a libboost_thread.a 
ln -s -f libfreetype-*gcc*.a  libfreetype.a

