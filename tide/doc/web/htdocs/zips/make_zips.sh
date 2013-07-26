#!/bin/sh
VERSION=1.0
cd $HOME/tide/doc/web/htdocs/zips || exit
zip -j tide-$VERSION-i686-binaries.zip ../i686/*
zip -j tide-$VERSION-x86_64-binaries.zip ../x86_64/*
zip -j tide-$VERSION-sample-data.zip ../data/* ../data/scripts/*
