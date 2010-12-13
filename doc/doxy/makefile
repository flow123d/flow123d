#
# $Id$
# $Revision$
# $LastChangedBy$
# $LastChangedDate$
#
# This makefile provides main rules for documentation

include ../../makefile.in

doc: flow123d.doc 

clean:
	rm -f flow123d/*

%.doc: %.doxy
	doxygen $<
