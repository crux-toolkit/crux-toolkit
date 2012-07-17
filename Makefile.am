# To set up debugging:
#
#   make clean
#   make CXXFLAGS=" -ggdb -O0 -DDEBUG " CFLAGS=" -ggdb -O0 -DDEBUG "
# 
# CXXFLAGS contains the options for the C++ compiler. The CFLAGS
# argument is for the C compiler, and is only needed if you are going
# to debug down into the MSTOOLKIT routines.  The -DDEBUG will cause
# the program to abort and dump core, rather than simply exit, when an
# error is encountered.

SUBDIRS = src doc

if SEARCH
  SEARCH_ENABLED = true
else
  SEARCH_ENABLED =
endif

distdir = $(PACKAGE)_$(VERSION)
distclean-local: clean-local

clean-local: clean-ac
clean-ac:
	rm -rf autom4te.cache

# makes a binary only distribution, using zip.
# make sure to 'make' before you make a 'bdist'
ar=`arch`
name=`uname`
bdist:
	rm -fr bdist bdist-data
	mkdir bdist
	mkdir bdist-data
	mkdir bdist/doc
	mkdir bdist/bin
	cp -R \
		doc/user/data/demo.ms2 \
		doc/user/data/small-yeast.fasta \
		doc/user/data/search.target.txt \
		doc/user/data/target.sqt \
		doc/user/broken-link.html \
		doc/user/compute-q-values.html \
		doc/user/crux-analyze-matches.html \
		doc/user/crux-create-index.html \
		doc/user/crux-generate-peptides.html \
		doc/user/crux-get-ms2-spectrum.html \
		doc/user/crux.html \
		doc/user/crux-predict-peptide-ions.html \
		doc/user/crux-search-for-matches.html \
		doc/user/customize-tutorial.html \
		doc/user/fasta-format.html \
		doc/user/index-tutorial.html \
		doc/user/install-tutorial.html \
		doc/user/ms2-format.html \
		doc/user/parameter-file.html \
		doc/user/percolator.html \
		doc/user/q-ranker.html \
		doc/user/release-notes.html \
		doc/user/search-tutorial.html \
		doc/user/sqt-format.html \
		doc/user/start-tutorial.html \
		doc/user/txt-format.html \
		doc/user/*.png \
		doc/user/*.gif \
		doc/user/default.params \
		doc/user/sample.fasta \
		doc/user/*.css \
		bdist/doc/
	cp 	src/c/crux \
		bdist/bin/
	if [ "$(SEARCH_ENABLED)" = 'true' -a "x"$(name) = "xLinux" ] ; then \
		mkdir bdist/tide; \
		mkdir bdist/tide/images; \
		cp tide/images/nibib.gif bdist/tide/images; \
		cp tide/images/UWlogo150p.gif bdist/tide/images; \
		cp tide/images/tide-performance.png bdist/tide/images; \
		cp tide/*.html bdist/tide; \
		cp tide/*.css bdist/tide; \
		cp tide/worm-demo*.sh bdist/tide; \
		cp tide/README bdist/tide; \
		cp tide/yeast-demo*.sh bdist/tide; \
		cp tide/$(ar)/tide-index bdist/tide; \
		cp tide/$(ar)/tide-search bdist/tide; \
		cp tide/$(ar)/msconvert bdist/tide; \
		cp tide/$(ar)/read-spectrumrecords bdist/tide; \
		scp nexus.gs.washington.edu:/net/noble/vol1/data/crux-datasets/tide-data-1/yeast-02-10000.spectrumrecords bdist-data/; \
		scp nexus.gs.washington.edu:/net/noble/vol1/data/crux-datasets/tide-data-1/yeast.fasta bdist-data/; \
		scp nexus.gs.washington.edu:/net/noble/vol1/data/crux-datasets/tide-data-1/worm-06-10000.spectrumrecords bdist-data/; \
		scp nexus.gs.washington.edu:/net/noble/vol1/data/crux-datasets/tide-data-1/worm.fasta bdist-data/;\
	fi
	zip -r $(distdir)-$(ar)-$(name).zip bdist
	rm -fr bdist
	if [ "x"$(name) = "xLinux" ] ; then \
		zip -r $(distdir)-data.zip bdist-data; \
		rm -rf bdist-data; \
	fi

checksum:
	md5sum $(distdir)-$(ar)-$(name).zip > $(distdir)-$(ar)-$(name).zip.md5sum.txt
	md5sum $(distdir).tar.gz > $(distdir).tar.gz.md5sum.txt

ignore-files:
	svn update .svnignore
	svn propset svn:ignore -F .svnignore .

# This is needed to keep 'make distcheck' from complaining about
# the file version.h, which doesn't get cleaned up (and doesn't
# need to be  cleaned up) in 'make distclean'
distcleancheck_listfiles = $(top_srcdir)/src/c/version.h

# Copy the appropriate copyright statement, COPYING for 
# limited license distribution, COPYRIGHT for open source distribution.
dist-hook:
	if test -z "$(SEARCH_ENABLED)"; \
	then cp COPYRIGHT $(distdir)/COPYING; fi
