AM_CXXFLAGS=$(DEBUG_OPTIONS) \
 -I$(top_srcdir)/src/c \
 -I$(top_srcdir)/src/c/ArrayLibrary \
 -I$(top_srcdir)/src/external/percolator/src \
 -I$(top_srcdir)/src/external/percolator/src/fido \
 -I$(top_srcdir)/src/external/xerces/src \
 -I/net/gs/vol3/software/modules-sw/xsd/3.3.0/Linux/RHEL5/x86_64/libxsd/ \
 -I/net/gs/vol3/software/modules-sw/boost/1.47.0/Linux/RHEL5/x86_64/include \
 -D_LARGEFILE_SOURCE

libpercolator_la_SOURCES = \
  PercolatorApplication.cpp\
  PercolatorApplication.h 

noinst_LTLIBRARIES = libpercolator.la

ignore-files:
	svn update .svnignore
	svn propset svn:ignore -F .svnignore .

