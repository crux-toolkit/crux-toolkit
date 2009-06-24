SUBDIRS = src doc

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
	rm -fr bdist
	mkdir bdist
	mkdir bdist/doc
	mkdir bdist/bin
	cp -R \
		doc/user/data/demo.ms2 \
		doc/user/data/small-yeast.fasta \
		doc/user/*.html \
		doc/user/*.png \
		doc/user/*.gif \
		doc/user/default.params \
		doc/user/*.css \
		bdist/doc/
	cp 	src/c/crux \
		src/c/crux-get-ms2-spectrum \
		src/c/crux-generate-peptides \
		src/c/crux-predict-peptide-ions \
		bin/spectacle.pl \
		bdist/bin/
	zip -r $(distdir)-$(ar)-$(name).zip bdist
	rm -fr bdist/

checksum:
	md5sum $(distdir)-$(ar)-$(name).zip > $(distdir)-$(ar)-$(name).zip.md5sum.txt
	md5sum $(distdir).tar.gz > $(distdir).tar.gz.md5sum.txt
	
ignore-files:
	svn update .svnignore
	svn propset svn:ignore -F .svnignore .
