#
# $Id$
# $Revision$
# $LastChangedBy$
# $LastChangedDate$
#
# This makefile provides main rules for documentation

doc: flow123d.doc 

clean:
	rm -f flow123d/*

%.doc: %.doxy
	cat flow123d.doxy short.doxy | doxygen - >doxygen.log
