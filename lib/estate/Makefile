# ------------------------------------------------------------------
# Makefile for ESTate.
# ------------------------------------------------------------------

all:	build documentation

clean:
	touch estate.tar.gz; rm estate.tar.gz;
	cd src; make clean; \
	cd ../doc; make clean;

vclean: clean
	rm -rf bin

documentation:
	cd doc; make;

build:
	cd src; make install;


tarfile:	vclean
	cd ..; tar cvf estate.tar estate; mv estate.tar estate;

dist:	tarfile
	gzip estate.tar

backup:	dist
	mv estate.tar.gz estate.tar.gz_`date +"%d_%m_%Y"`_backup

# --
